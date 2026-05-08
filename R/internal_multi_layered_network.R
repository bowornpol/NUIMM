#' High-Tech Multi-Omics Network Visualization
#' @keywords internal
utils::globalVariables(c("type", "weight", "name", "edge_score", "layer"))

# --- Futuristic Taxonomy Cleaner ---
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

  # Data Ingestion
  mpn <- read_input_file(mpn_file, file_type = "csv", stringsAsFactors = FALSE)
  ppn <- if (!is.na(ppn_file) && file.exists(ppn_file)) read_input_file(ppn_file, file_type = "csv") else NULL
  pmn <- if (!is.null(pmn_file) && file.exists(pmn_file)) read_input_file(pmn_file, file_type = "csv", stringsAsFactors = FALSE) else NULL
  gsea <- read_input_file(gsea_file, file_type = "csv")

  valid_paths <- gsea$ID
  edges <- data.frame()

  # Build Edge Logic
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
    library(igraph)
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% edges$Feature1[edges$edge_type == "Microbe-Pathway"], "Microbe",
      ifelse(igraph::V(g)$name %in% edges$Feature2[edges$edge_type == "Pathway-Metabolite"], "Metabolite", "Pathway")
    )

    if (visualize) {
      tryCatch(
        {
          # 1. Advanced Node Prep
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

          # 2. Smart Spacing & Alignment Logic
          types <- c("Microbe", "Pathway", "Metabolite")
          x_pos <- c(-700, 0, 700)
          nodes_df$x <- 0
          nodes_df$y <- 0

          for (i in 1:3) {
            idx <- which(nodes_df$group == types[i])
            n <- length(idx)
            if (n > 0) {
              # Distribute nodes vertically with random jitter to prevent perfect lines and overlaps
              nodes_df$x[idx] <- x_pos[i] + runif(n, -50, 50)
              nodes_df$y[idx] <- seq(-400, 400, length.out = n)
            }
          }

          # Cinematic Aesthetic Mapping
          nodes_df$shape <- c("Microbe" = "hexagon", "Pathway" = "dot", "Metabolite" = "diamond")[nodes_df$group]
          nodes_df$color <- c("Microbe" = "#00d2ff", "Pathway" = "#9d50bb", "Metabolite" = "#f2994a")[nodes_df$group]
          nodes_df$shadow <- TRUE

          edges_df <- data.frame(from = edges$Feature1, to = edges$Feature2, value = edges$edge_score)

          # 3. High-Tech UI Implementation
          vis_plot <- visNetwork::visNetwork(nodes_df, edges_df, width = "100%", height = "100vh") |>
            visNetwork::visNodes(
              font = list(color = "#2c3e50", size = 20, face = "Orbitron, sans-serif", background = "rgba(255,255,255,0.7)"),
              borderWidth = 2,
              borderWidthSelected = 6,
              scaling = list(min = 20, max = 50)
            ) |>
            visNetwork::visEdges(
              smooth = list(enabled = TRUE, type = "diagonalCross"),
              color = list(color = "rgba(44, 62, 80, 0.15)", highlight = "#ff0000", hover = "#ff0000"),
              width = 2
            ) |>
            visNetwork::visInteraction(
              dragNodes = TRUE, dragView = TRUE, zoomView = TRUE,
              multiselect = TRUE, navigationButtons = TRUE, hover = TRUE
            ) |>
            # SMART COLLISION AVOIDANCE ENGINE
            visNetwork::visPhysics(
              enabled = TRUE,
              solver = "barnesHut",
              barnesHut = list(
                gravitationalConstant = -8000,
                centralGravity = 0.3,
                springLength = 150,
                springConstant = 0.05,
                avoidOverlap = 1 # CRITICAL: Anti-overlap setting
              ),
              stabilization = list(iterations = 200)
            ) |>
            # FULL CUSTOMIZATION SYSTEM (Left Panel)
            visNetwork::visConfigure(
              enabled = TRUE,
              filter = c("nodes", "edges", "physics"),
              showButton = FALSE
            ) |>
            visNetwork::visExport(
              type = "png", label = "💾 SAVE NETWORK",
              style = "position:absolute; left:20px; bottom:80px; background:linear-gradient(45deg, #2c3e50, #000000); color:#00d2ff; padding:15px 25px; border-radius:5px; border:1px solid #00d2ff; cursor:pointer; font-family:monospace; font-weight:bold; letter-spacing:2px; box-shadow: 0 0 15px rgba(0,210,255,0.4);"
            )

          vis_plot$x$background <- "#ffffff" # Default clean futuristic white

          html_path <- file.path(output_dir, paste0("interactive_mln_", tools::file_path_sans_ext(basename(gsea_file)), ".html"))
          htmlwidgets::saveWidget(vis_plot, file = html_path, selfcontained = TRUE, background = "#ffffff")
        },
        error = function(e) warning("Visualization failed: ", e$message)
      )
    }
  }
  return(out_path)
}
