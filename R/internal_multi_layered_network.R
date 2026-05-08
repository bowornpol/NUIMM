#' Internal MLN Integration & Interactive Visualization Helper
#' @keywords internal
utils::globalVariables(c("type", "weight", "name", "edge_score", "layer"))

# --- Helper Function to Clean Microbial Taxonomy ---
clean_taxonomy <- function(tax_string) {
  if (grepl("g__", tax_string)) {
    # Split by semicolon (with or without spaces)
    parts <- unlist(strsplit(tax_string, ";\\s*|;"))

    g_part <- parts[grepl("^g__", parts)][1]
    s_part <- parts[grepl("^s__", parts)][1]

    # Check if species exists and is actually named (longer than just "s__")
    if (!is.na(s_part) && nchar(s_part) > 3) {
      return(paste(g_part, s_part, sep = " "))
    } else if (!is.na(g_part)) {
      return(g_part)
    }
  }
  return(tax_string) # Return original if it doesn't match expected taxonomy pattern
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

    # Strict Biological Layering
    igraph::V(g)$layer <- ifelse(igraph::V(g)$type == "Microbe", 1, ifelse(igraph::V(g)$type == "Pathway", 2, 3))

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
              label = unname(node_names_clean), # Uses the cleaned short names
              group = igraph::V(g)$type,
              level = igraph::V(g)$layer,
              title = paste0("<div style='padding:5px; font-family:sans-serif;'><b>Type:</b> ", igraph::V(g)$type, "<br><b>ID:</b> ", igraph::V(g)$name, "</div>"),
              stringsAsFactors = FALSE
            )

            # Custom Modern Aesthetic Map
            shape_map <- c("Microbe" = "hexagon", "Pathway" = "dot", "Metabolite" = "diamond")
            color_map <- c("Microbe" = "#00cec9", "Pathway" = "#6c5ce7", "Metabolite" = "#fdcb6e")

            nodes_df$shape <- shape_map[nodes_df$group]
            nodes_df$color <- color_map[nodes_df$group]

            # 2. Prepare Edges
            edges_df <- data.frame(
              from = edges$Feature1,
              to = edges$Feature2,
              value = edges$edge_score,
              title = paste0("<div style='padding:5px; font-family:sans-serif;'><b>Edge Score:</b> ", round(edges$edge_score, 4), "</div>")
            )

            # 3. Build Interactive Plot (Frozen & Beautiful)
            vis_plot <- visNetwork::visNetwork(
              nodes_df, edges_df,
              width = "100%", height = "900px",
              main = list(text = "Multi-Omics Mechanistic Network", style = "color: #f8fafc; font-family: sans-serif; font-weight: 400;")
            ) |>
              visNetwork::visNodes(
                font = list(color = "#cbd5e1", size = 18, face = "sans-serif"),
                borderWidth = 2,
                borderWidthSelected = 6,
                scaling = list(min = 15, max = 35),
                shadow = list(enabled = TRUE, color = "rgba(0,0,0,0.6)", size = 10, x = 3, y = 3)
              ) |>
              visNetwork::visEdges(
                smooth = list(enabled = TRUE, type = "cubicBezier", roundness = 0.5), # Elegant curved lines
                color = list(color = "rgba(148, 163, 184, 0.2)", highlight = "#ff7675", hover = "#ff7675"), # Lights up pink on click
                selectionWidth = 3
              ) |>
              visNetwork::visOptions(
                highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE, hideColor = "rgba(15, 23, 42, 0.3)"),
                nodesIdSelection = list(enabled = TRUE, style = "width: 250px; background: #334155; color: white; border: none; border-radius: 5px; padding: 5px;"),
                selectedBy = list(variable = "group", style = "width: 250px; background: #334155; color: white; border: none; border-radius: 5px; padding: 5px;")
              ) |>
              visNetwork::visInteraction(
                navigationButtons = TRUE,
                dragNodes = TRUE,
                dragView = TRUE,
                zoomView = TRUE,
                hover = TRUE
              ) |>
              visNetwork::visHierarchicalLayout(
                direction = "UD", # Up-Down orientation
                levelSeparation = 250, # Space between Omics layers
                nodeSpacing = 80, # Space between nodes in the same layer
                edgeMinimization = TRUE,
                parentCentralization = TRUE
              ) |>
              # THIS COMMAND KILLS THE PHYSICS JIGGLE:
              visNetwork::visPhysics(enabled = FALSE) |>
              visNetwork::visExport(
                type = "png",
                name = paste0("network_", tools::file_path_sans_ext(basename(gsea_file))),
                label = "📥 Download High-Res PNG",
                style = "background: #00cec9; color: #0f172a; border: none; padding: 8px 15px; border-radius: 5px; cursor: pointer; font-weight: bold;"
              )

            vis_plot$x$background <- "#0f172a" # Deep Slate Blue background

            # 4. Save Self-Contained HTML
            html_path <- file.path(output_dir, paste0("interactive_mln_", tools::file_path_sans_ext(basename(gsea_file)), ".html"))
            htmlwidgets::saveWidget(vis_plot, file = html_path, selfcontained = TRUE, background = "#0f172a")
            message("  Saved Beautiful Interactive HTML visualization to: ", html_path)
          }
        },
        error = function(e) warning("Interactive visualization failed: ", e$message)
      )
    }
  }
  return(out_path)
}
