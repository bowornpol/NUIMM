#' Internal MLN Integration & Interactive Visualization Helper
#' @keywords internal
utils::globalVariables(c("type", "weight", "name", "edge_score", "layer"))

con_mln_int <- function(
  gsea_file, mpn_file, ppn_file, pmn_file, output_dir,
  visualize, layout_method, node_colors, node_shapes, base_node_size, plot_width, plot_height, plot_dpi
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  mpn <- read_input_file(mpn_file, file_type = "csv", stringsAsFactors = FALSE)
  ppn <- if (!is.na(ppn_file) && file.exists(ppn_file)) read_input_file(ppn_file, file_type = "csv") else NULL
  pmn <- if (!is.null(pmn_file) && file.exists(pmn_file)) read_input_file(pmn_file, file_type = "csv", stringsAsFactors = FALSE) else NULL
  gsea <- read_input_file(gsea_file, file_type = "csv")

  valid_paths <- gsea$ID
  edges <- data.frame()

  mpn <- mpn[mpn$FunctionID %in% valid_paths, ]
  if (nrow(mpn) > 0) edges <- rbind(edges, data.frame(Feature1 = mpn$TaxonID, Feature2 = mpn$FunctionID, edge_score = mpn$relative_contribution, edge_type = "Microbe-Pathway"))

  if (!is.null(ppn)) {
    ppn <- ppn[ppn$FunctionID_1 %in% valid_paths & ppn$FunctionID_2 %in% valid_paths, ]
    if (nrow(ppn) > 0) edges <- rbind(edges, data.frame(Feature1 = ppn$FunctionID_1, Feature2 = ppn$FunctionID_2, edge_score = ppn$jaccard_index, edge_type = "Pathway-Pathway"))
  }

  if (!is.null(pmn)) {
    pmn <- pmn[pmn$FunctionID %in% valid_paths, ]
    if (nrow(pmn) > 0) edges <- rbind(edges, data.frame(Feature1 = pmn$FunctionID, Feature2 = pmn$MetaboliteID, edge_score = abs(pmn$correlation), edge_type = "Pathway-Metabolite"))
  }

  out_path <- file.path(output_dir, paste0("final_mln_", basename(gsea_file)))
  write.csv(edges, out_path, row.names = FALSE)

  if (nrow(edges) > 0) {
    library(igraph)
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)

    igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% edges$Feature1[edges$edge_type == "Microbe-Pathway"], "Microbe",
      ifelse(igraph::V(g)$name %in% edges$Feature2[edges$edge_type == "Pathway-Metabolite"], "Metabolite", "Pathway")
    )
    igraph::V(g)$layer <- ifelse(igraph::V(g)$type == "Microbe", 1, ifelse(igraph::V(g)$type == "Pathway", 2, 3))

    # --- ADJUSTMENT: Always export GraphML for large-scale external visualization (Cytoscape/Gephi)
    igraph::write_graph(g, file.path(output_dir, paste0("final_mln_", tools::file_path_sans_ext(basename(gsea_file)), ".graphml")), format = "graphml")

    if (visualize) {
      tryCatch(
        {
          if (!requireNamespace("visNetwork", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) {
            warning("Install 'visNetwork' and 'htmlwidgets' for interactive HTML plotting.")
          } else {
            # --- ADJUSTMENT: Build Interactive Plot instead of static ggraph
            nodes_df <- data.frame(
              id = igraph::V(g)$name,
              label = igraph::V(g)$name,
              group = igraph::V(g)$type,
              level = igraph::V(g)$layer, # Hierarchical layering
              title = paste0("<p><b>Type:</b> ", igraph::V(g)$type, "<br><b>ID:</b> ", igraph::V(g)$name, "</p>"), # Hover tooltip
              stringsAsFactors = FALSE
            )

            shape_map <- c("Microbe" = "triangle", "Pathway" = "dot", "Metabolite" = "square")
            nodes_df$shape <- shape_map[nodes_df$group]

            # Map colors using the user's provided palette
            color_map <- c("Microbe" = "#D55E00", "Pathway" = "#0072B2", "Metabolite" = "#009E73")
            nodes_df$color <- color_map[nodes_df$group]

            edges_df <- data.frame(
              from = edges$Feature1,
              to = edges$Feature2,
              value = edges$edge_score, # Scales edge thickness
              title = paste0("<p><b>Score:</b> ", round(edges$edge_score, 4), "</p>"), # Hover tooltip
              color = "rgba(200, 200, 200, 0.4)" # Semi-transparent light gray
            )

            # Build Plot (Dark Academic Presentation Theme)
            vis_plot <- visNetwork::visNetwork(
              nodes_df, edges_df,
              width = "100%", height = "900px",
              main = list(text = "Multi-Layered Multi-Omics Network", style = "color: white; font-family: sans-serif;")
            ) |>
              visNetwork::visNodes(
                font = list(color = "white", size = 20, face = "sans-serif"),
                borderWidth = 1.5,
                borderWidthSelected = 4,
                scaling = list(min = 10, max = 30)
              ) |>
              visNetwork::visEdges(smooth = list(enabled = TRUE, type = "continuous")) |>
              visNetwork::visOptions(
                highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), # Highlights connected nodes on click
                nodesIdSelection = list(enabled = TRUE, style = "width: 200px; background: #333; color: white;"), # Dropdown search
                selectedBy = list(variable = "group", style = "width: 200px; background: #333; color: white;") # Filter by omic layer
              ) |>
              visNetwork::visInteraction(navigationButtons = TRUE, dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) |>
              visNetwork::visPhysics(stabilization = TRUE, barnesHut = list(gravitationalConstant = -6000, centralGravity = 0.3)) |>
              visNetwork::visExport(type = "png", name = paste0("network_", tools::file_path_sans_ext(basename(gsea_file))), label = "Save Current View as PNG")

            if (layout_method == "sugiyama") {
              vis_plot <- vis_plot |> visNetwork::visHierarchicalLayout(direction = "UD", levelSeparation = 300, nodeSpacing = 150)
            }

            vis_plot$x$background <- "#121212"

            # Save Self-Contained HTML
            html_path <- file.path(output_dir, paste0("interactive_mln_", tools::file_path_sans_ext(basename(gsea_file)), ".html"))
            htmlwidgets::saveWidget(vis_plot, file = html_path, selfcontained = TRUE, background = "#121212")
            message("  Saved Interactive HTML visualization to: ", html_path)
          }
        },
        error = function(e) warning("Interactive visualization failed: ", e$message)
      )
    }
  }
  return(out_path)
}
