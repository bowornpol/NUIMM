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

  # --- FIX: Define out_path immediately so it is never 'not found' ---
  out_path <- file.path(output_dir, paste0("final_mln_data_", basename(gsea_file)))

  mpn <- read_input_file(mpn_file, file_type = "csv", stringsAsFactors = FALSE)
  ppn <- if (!is.na(ppn_file) && file.exists(ppn_file)) read_input_file(ppn_file, file_type = "csv") else NULL
  pmn <- if (!is.null(pmn_file) && file.exists(pmn_file)) read_input_file(pmn_file, file_type = "csv", stringsAsFactors = FALSE) else NULL
  gsea <- read_input_file(gsea_file, file_type = "csv")

  valid_paths <- gsea$ID
  edges <- data.frame()

  # Microbe-Pathway
  mpn <- mpn[mpn$FunctionID %in% valid_paths, ]
  if (nrow(mpn) > 0) edges <- rbind(edges, data.frame(Feature1 = mpn$TaxonID, Feature2 = mpn$FunctionID, edge_score = mpn$relative_contribution, edge_type = "Microbe-Pathway"))

  # Pathway-Pathway
  if (!is.null(ppn)) {
    ppn <- ppn[ppn$FunctionID_1 %in% valid_paths & ppn$FunctionID_2 %in% valid_paths, ]
    if (nrow(ppn) > 0) edges <- rbind(edges, data.frame(Feature1 = ppn$FunctionID_1, Feature2 = ppn$FunctionID_2, edge_score = ppn$jaccard_index, edge_type = "Pathway-Pathway"))
  }

  # Pathway-Metabolite
  if (!is.null(pmn)) {
    pmn <- pmn[pmn$FunctionID %in% valid_paths, ]
    if (nrow(pmn) > 0) edges <- rbind(edges, data.frame(Feature1 = pmn$FunctionID, Feature2 = pmn$MetaboliteID, edge_score = abs(pmn$correlation), edge_type = "Pathway-Metabolite"))
  }

  # Save the raw edge data
  if (nrow(edges) > 0) {
    write.csv(edges, out_path, row.names = FALSE)

    library(igraph)
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)

    igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% edges$Feature1[edges$edge_type == "Microbe-Pathway"], "Microbe",
      ifelse(igraph::V(g)$name %in% edges$Feature2[edges$edge_type == "Pathway-Metabolite"], "Metabolite", "Pathway")
    )

    # Save GraphML for Cytoscape
    igraph::write_graph(g, file.path(output_dir, paste0("network_architecture_", tools::file_path_sans_ext(basename(gsea_file)), ".graphml")), format = "graphml")

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
            # Strict level mapping for column layout
            level = ifelse(igraph::V(g)$type == "Microbe", 1, ifelse(igraph::V(g)$type == "Pathway", 2, 3)),
            title = paste0("<div style='padding:8px; font-family:sans-serif;'><b>ID:</b> ", igraph::V(g)$name, "</div>"),
            stringsAsFactors = FALSE
          )

          # Advanced Colors & Shapes
          nodes_df$shape <- c("Microbe" = "hexagon", "Pathway" = "dot", "Metabolite" = "diamond")[nodes_df$group]
          nodes_df$color <- c("Microbe" = "#0ea5e9", "Pathway" = "#8b5cf6", "Metabolite" = "#f59e0b")[nodes_df$group]

          edges_df <- data.frame(from = edges$Feature1, to = edges$Feature2, value = edges$edge_score)

          # 3. Build Futuristic Plot
          vis_plot <- visNetwork::visNetwork(nodes_df, edges_df, width = "100%", height = "100vh") |>
            visNetwork::visNodes(
              font = list(color = "#1e293b", size = 20, face = "sans-serif", background = "rgba(255,255,255,0.7)"),
              borderWidth = 1.5,
              shadow = list(enabled = TRUE, color = "rgba(0,0,0,0.1)")
            ) |>
            visNetwork::visEdges(
              smooth = list(enabled = TRUE, type = "diagonalCross"),
              color = list(color = "rgba(148, 163, 184, 0.4)", highlight = "#e11d48")
            ) |>
            visNetwork::visHierarchicalLayout(
              direction = "LR", # Left to Right for your columns
              levelSeparation = 400, # Wide space between Microbe, Pathway, Metabolite
              nodeSpacing = 100, # Space to prevent overlap within columns
              sortMethod = "directed"
            ) |>
            visNetwork::visInteraction(
              navigationButtons = TRUE,
              dragNodes = TRUE,
              multiselect = TRUE,
              hover = TRUE
            ) |>
            visNetwork::visPhysics(enabled = FALSE) |>
            visNetwork::visExport(
              type = "png", label = "💾 SAVE NETWORK",
              style = "position:absolute; left:20px; bottom:40px; background:#0f172a; color:#0ea5e9; padding:15px 25px; border-radius:8px; border:1px solid #0ea5e9; cursor:pointer; font-weight:bold; font-family:sans-serif; z-index:1000;"
            )

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
