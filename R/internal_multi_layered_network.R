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

  # Ensure path is defined early to prevent 'out_path not found' error
  out_path <- file.path(output_dir, paste0("final_mln_data_", basename(gsea_file)))

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

  if (nrow(edges) > 0) {
    write.csv(edges, out_path, row.names = FALSE)
    library(igraph)
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)

    igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% edges$Feature1[edges$edge_type == "Microbe-Pathway"], "Microbe",
      ifelse(igraph::V(g)$name %in% edges$Feature2[edges$edge_type == "Pathway-Metabolite"], "Metabolite", "Pathway")
    )

    if (visualize) {
      tryCatch(
        {
          # 1. Prepare Nodes
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
            title = paste0("<div style='padding:8px; font-family:sans-serif;'><b>ID:</b> ", igraph::V(g)$name, "</div>"),
            stringsAsFactors = FALSE
          )

          # Initialize Coordinates for 3-Zone Cluster
          types <- c("Microbe", "Pathway", "Metabolite")
          x_offsets <- c(-700, 0, 700)
          nodes_df$x <- 0
          nodes_df$y <- 0
          for (i in 1:3) {
            idx <- which(nodes_df$group == types[i])
            if (length(idx) > 0) {
              nodes_df$x[idx] <- x_offsets[i] + runif(length(idx), -50, 50)
              nodes_df$y[idx] <- seq(-400, 400, length.out = length(idx))
            }
          }

          nodes_df$shape <- c("Microbe" = "hexagon", "Pathway" = "dot", "Metabolite" = "diamond")[nodes_df$group]
          nodes_df$color <- c("Microbe" = "#0ea5e9", "Pathway" = "#8b5cf6", "Metabolite" = "#f59e0b")[nodes_df$group]

          edges_df <- data.frame(from = edges$Feature1, to = edges$Feature2, value = edges$edge_score)

          # 3. Build Plot
          vis_plot <- visNetwork::visNetwork(nodes_df, edges_df, width = "100%", height = "100vh") |>
            visNetwork::visNodes(
              font = list(color = "#1e293b", size = 18, background = "rgba(255,255,255,0.7)"),
              borderWidth = 1.5, shadow = list(enabled = TRUE)
            ) |>
            visNetwork::visEdges(
              smooth = list(enabled = TRUE, type = "continuous"),
              color = list(color = "rgba(148, 163, 184, 0.4)", highlight = "#e11d48")
            ) |>
            visNetwork::visInteraction(
              navigationButtons = TRUE, dragNodes = TRUE, multiselect = TRUE, hover = TRUE
            ) |>
            # ADVANCED PHYSICS FOR PRACTICAL USE (ANTI-OVERLAP)
            visNetwork::visPhysics(
              enabled = TRUE,
              solver = "forceAtlas2Based",
              forceAtlas2Based = list(
                gravitationalConstant = -150,
                centralGravity = 0.01,
                springLength = 100,
                springConstant = 0.08,
                avoidOverlap = 1 # THIS STOPS NODES FROM OVERLAPPING
              ),
              stabilization = list(enabled = TRUE, iterations = 200)
            ) |>
            visNetwork::visExport(
              type = "png", label = "💾 SAVE NETWORK",
              style = "position:absolute; left:20px; bottom:40px; background:#1e293b; color:white; padding:12px 20px; border-radius:8px; border:none; cursor:pointer; font-weight:bold; z-index:1000;"
            ) |>
            visNetwork::visConfigure(enabled = TRUE, filter = c("physics", "layout"))

          vis_plot$x$background <- "#ffffff"

          html_path <- file.path(output_dir, paste0("interactive_mln_", tools::file_path_sans_ext(basename(gsea_file)), ".html"))
          htmlwidgets::saveWidget(vis_plot, file = html_path, selfcontained = TRUE, background = "#ffffff")
        },
        error = function(e) warning("Visualization failed: ", e$message)
      )
    }
  }
  return(out_path)
}
