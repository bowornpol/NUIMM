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

    if (visualize) {
      tryCatch(
        {
          # 1. Node data preparation
          nodes_df <- data.frame(
            id = igraph::V(g)$name,
            label = unname(sapply(igraph::V(g)$name, function(x) {
              if (igraph::V(g)$type[igraph::V(g)$name == x] == "Microbe") {
                return(clean_taxonomy(x))
              } else {
                return(x)
              }
            })),
            group = igraph::V(g)$type,
            stringsAsFactors = FALSE
          )

          # 2. CALCULATION FOR TRI-CIRCLE LAYOUT
          # Manually setting offsets: Left (-600), Center (0), Right (600)
          types <- c("Microbe", "Pathway", "Metabolite")
          x_offsets <- c(-600, 0, 600)

          nodes_df$x <- 0
          nodes_df$y <- 0

          for (i in 1:3) {
            idx <- which(nodes_df$group == types[i])
            n <- length(idx)
            if (n > 0) {
              angle <- seq(0, 2 * pi, length.out = n + 1)[1:n]
              radius <- 180 + (n * 1.5)
              nodes_df$x[idx] <- x_offsets[i] + radius * cos(angle)
              nodes_df$y[idx] <- radius * sin(angle)
            }
          }

          # Styling for White Background
          nodes_df$shape <- c("Microbe" = "hexagon", "Pathway" = "dot", "Metabolite" = "diamond")[nodes_df$group]
          nodes_df$color <- c("Microbe" = "#0ea5e9", "Pathway" = "#8b5cf6", "Metabolite" = "#f59e0b")[nodes_df$group]

          edges_df <- data.frame(from = edges$Feature1, to = edges$Feature2, value = edges$edge_score)

          # 3. Build Simplified Plot with Right-Side Legend
          vis_plot <- visNetwork::visNetwork(nodes_df, edges_df, width = "100%", height = "100vh") |>
            visNetwork::visNodes(
              font = list(color = "#1e293b", size = 18, face = "sans-serif"),
              shadow = list(enabled = TRUE, color = "rgba(0,0,0,0.1)", size = 10)
            ) |>
            visNetwork::visEdges(
              smooth = list(enabled = TRUE, type = "curvedCW", roundness = 0.2),
              color = list(color = "rgba(180, 180, 180, 0.4)", highlight = "#e11d48")
            ) |>
            # Define Groups for the Legend
            visNetwork::visGroups(groupname = "Microbe", color = "#0ea5e9", shape = "hexagon") |>
            visNetwork::visGroups(groupname = "Pathway", color = "#8b5cf6", shape = "dot") |>
            visNetwork::visGroups(groupname = "Metabolite", color = "#f59e0b", shape = "diamond") |>
            visNetwork::visLegend(main = "Node Key", position = "right", width = 0.1) |>
            visNetwork::visInteraction(
              dragNodes = TRUE, dragView = TRUE, zoomView = TRUE, multiselect = TRUE, navigationButtons = TRUE
            ) |>
            visNetwork::visPhysics(enabled = FALSE) |>
            visNetwork::visExport(
              type = "png", label = "📸 Save Layout",
              style = "position:absolute; left:30px; bottom:30px; background:#1e293b; color:white; padding:15px 25px; border-radius:12px; border:none; cursor:pointer; font-weight:bold; font-family:sans-serif; box-shadow: 0 4px 10px rgba(0,0,0,0.2); z-index:1001;"
            )

          vis_plot$x$background <- "#ffffff"

          html_path <- file.path(output_dir, paste0("interactive_mln_", tools::file_path_sans_ext(basename(gsea_file)), ".html"))
          htmlwidgets::saveWidget(vis_plot, file = html_path, selfcontained = TRUE, background = "#ffffff")
          message("Successfully generated clean Tri-Cluster network at: ", html_path)
        },
        error = function(e) warning("Visualization failed: ", e$message)
      )
    }
  }
  return(out_path)
}
