#' Internal MLN Integration & Interactive Visualization Helper
#' @keywords internal
utils::globalVariables(c("type", "weight", "name", "edge_score", "layer"))

# --- Helper Function to Clean Microbial Taxonomy ---
clean_taxonomy <- function(tax_string) {
  if (grepl("g__", tax_string)) {
    parts <- unlist(strsplit(tax_string, ";\\s*|;"))
    g_part <- parts[grepl("^g__", parts)][1]
    s_part <- parts[grepl("^s__", parts)][1]

    if (!is.na(s_part) && nchar(s_part) > 3) {
      return(paste(g_part, s_part, sep = " "))
    } else if (!is.na(g_part)) {
      return(g_part)
    }
  }
  return(tax_string)
}

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

    igraph::write_graph(g, file.path(output_dir, paste0("final_mln_", tools::file_path_sans_ext(basename(gsea_file)), ".graphml")), format = "graphml")

    if (visualize) {
      tryCatch(
        {
          if (!requireNamespace("visNetwork", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) {
            warning("Install 'visNetwork' and 'htmlwidgets' for interactive HTML plotting.")
          } else {
            # 1. Prepare Nodes & Apply Taxonomy Cleaner
            node_names_clean <- sapply(igraph::V(g)$name, function(x) {
              type <- igraph::V(g)$type[igraph::V(g)$name == x]
              if (type == "Microbe") {
                return(clean_taxonomy(x))
              } else {
                return(x)
              }
            })

            nodes_df <- data.frame(
              id = igraph::V(g)$name,
              label = unname(node_names_clean),
              group = igraph::V(g)$type,
              title = paste0("<div style='padding:8px; font-family:sans-serif; background:white; color:black; border-radius:5px;'><b>Type:</b> ", igraph::V(g)$type, "<br><b>ID:</b> ", igraph::V(g)$name, "</div>"),
              stringsAsFactors = FALSE
            )

            # EXPLICIT HIERARCHY: Microbe (1) -> Pathway (2) -> Metabolite (3)
            nodes_df$level <- ifelse(nodes_df$group == "Microbe", 1,
              ifelse(nodes_df$group == "Pathway", 2, 3)
            )

            # Colors & Shapes
            shape_map <- c("Microbe" = "hexagon", "Pathway" = "dot", "Metabolite" = "diamond")
            color_map <- c("Microbe" = "#0ea5e9", "Pathway" = "#8b5cf6", "Metabolite" = "#f59e0b")

            nodes_df$shape <- shape_map[nodes_df$group]
            nodes_df$color <- color_map[nodes_df$group]

            edges_df <- data.frame(
              from = edges$Feature1,
              to = edges$Feature2,
              value = edges$edge_score,
              title = paste0("<div style='padding:5px; font-family:sans-serif;'><b>Edge Score:</b> ", round(edges$edge_score, 4), "</div>")
            )

            # 3. Build Interactive Plot
            # Using 100vh height ensures the canvas perfectly fits the browser window (no cut-off menus)
            vis_plot <- visNetwork::visNetwork(
              nodes_df, edges_df,
              width = "100%", height = "100vh"
            ) |>
              visNetwork::visNodes(
                font = list(color = "#1e293b", size = 16, face = "sans-serif", strokeWidth = 2, strokeColor = "#ffffff"),
                borderWidth = 1.5,
                borderWidthSelected = 5,
                scaling = list(min = 15, max = 35),
                shadow = list(enabled = TRUE, color = "rgba(0,0,0,0.15)", size = 8, x = 2, y = 2)
              ) |>
              visNetwork::visEdges(
                smooth = list(enabled = TRUE, type = "cubicBezier", roundness = 0.5),
                color = list(color = "rgba(148, 163, 184, 0.6)", highlight = "#e11d48", hover = "#e11d48"),
                selectionWidth = 3
              ) |>
              visNetwork::visHierarchicalLayout(
                direction = "UD", # Up-Down
                levelSeparation = 250, # Vertical gap between the 3 layers
                nodeSpacing = 80, # Horizontal gap between nodes in the same layer
                sortMethod = "directed"
              ) |>
              visNetwork::visOptions(
                highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE, hideColor = "rgba(255, 255, 255, 0.8)"),
                selectedBy = list(variable = "group", multiple = TRUE, main = "🔍 Highlight by Group")
              ) |>
              visNetwork::visInteraction(
                navigationButtons = TRUE, # Shows standard zoom/pan controls
                dragNodes = TRUE,
                dragView = TRUE,
                zoomView = TRUE,
                hover = TRUE,
                multiselect = TRUE # ENABLES CTRL + DRAG TO SELECT MULTIPLE NODES
              ) |>
              visNetwork::visPhysics(enabled = FALSE) |> # KILLS THE JIGGLE. IT WILL STAY STILL.
              visNetwork::visExport(
                type = "png",
                name = paste0("network_", tools::file_path_sans_ext(basename(gsea_file))),
                label = "📸 Save Network Image",
                # Custom CSS to pin the button to the bottom left corner
                style = "position: absolute; left: 20px; bottom: 30px; background: #0f172a; color: white; border: none; padding: 12px 24px; border-radius: 8px; cursor: pointer; font-weight: bold; font-family: sans-serif; box-shadow: 0 4px 6px rgba(0,0,0,0.2); z-index: 1000;"
              )

            vis_plot$x$background <- "#ffffff"

            # 4. Save Self-Contained HTML
            html_path <- file.path(output_dir, paste0("interactive_mln_", tools::file_path_sans_ext(basename(gsea_file)), ".html"))
            htmlwidgets::saveWidget(vis_plot, file = html_path, selfcontained = TRUE, background = "#ffffff")
            message("  Saved Beautiful Interactive HTML visualization to: ", html_path)
          }
        },
        error = function(e) warning("Interactive visualization failed: ", e$message)
      )
    }
  }
  return(out_path)
}
