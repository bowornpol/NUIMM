#' High-Tech Multi-Omics Network Visualization Helper
#' @keywords internal
utils::globalVariables(c("type", "weight", "name", "edge_score", "layer"))

# --- Helper to Clean Name and Extract Family (f__) ---
clean_and_get_family <- function(tax_string) {
  result <- list(clean = tax_string, family = "Unassigned")
  if (grepl("g__", tax_string)) {
    parts <- unlist(strsplit(tax_string, ";\\s*|;"))
    g_part <- parts[grepl("^g__", parts)][1]
    s_part <- parts[grepl("^s__", parts)][1]
    if (!is.na(s_part) && nchar(s_part) > 3) {
      result$clean <- paste(g_part, s_part, sep = " ")
    } else if (!is.na(g_part)) {
      result$clean <- g_part
    }
    f_part <- parts[grepl("^f__", parts)][1]
    if (!is.na(f_part)) result$family <- f_part
  }
  return(result)
}

con_mln_int <- function(
  gsea_file, mpn_file, ppn_file, pmn_file, output_dir,
  visualize, layout_method, node_colors, node_shapes, base_node_size, plot_width, plot_height, plot_dpi
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Define paths immediately
  base_name <- tools::file_path_sans_ext(basename(gsea_file))
  out_csv <- file.path(output_dir, paste0("final_mln_data_", base_name, ".csv"))
  out_html <- file.path(output_dir, paste0("interactive_mln_", base_name, ".html"))

  mpn <- read_input_file(mpn_file, file_type = "csv", stringsAsFactors = FALSE)
  ppn <- if (!is.na(ppn_file) && file.exists(ppn_file)) read_input_file(ppn_file, file_type = "csv") else NULL
  pmn <- if (!is.null(pmn_file) && file.exists(pmn_file)) read_input_file(pmn_file, file_type = "csv", stringsAsFactors = FALSE) else NULL
  gsea <- read_input_file(gsea_file, file_type = "csv")

  valid_paths <- gsea$ID
  edges <- data.frame()

  # Building Edges
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
    write.csv(edges, out_csv, row.names = FALSE)
    library(igraph)
    library(visNetwork)
    library(htmlwidgets)
    g <- graph_from_data_frame(edges, directed = FALSE)
    V(g)$type <- ifelse(V(g)$name %in% edges$Feature1[edges$edge_type == "Microbe-Pathway"], "Microbe",
      ifelse(V(g)$name %in% edges$Feature2[edges$edge_type == "Pathway-Metabolite"], "Metabolite", "Pathway")
    )

    if (visualize) {
      tryCatch(
        {
          # Node Prep with Family extraction
          tax_info <- lapply(V(g)$name, function(x) if (V(g)$type[V(g)$name == x] == "Microbe") clean_and_get_family(x) else list(clean = x, family = "NA"))
          nodes_df <- data.frame(id = V(g)$name, label = sapply(tax_info, function(x) x$clean), family = sapply(tax_info, function(x) x$family), group = V(g)$type, stringsAsFactors = FALSE)

          # 3-Zone Layout Logic
          types <- c("Microbe", "Pathway", "Metabolite")
          x_zones <- c(-850, 0, 850)
          nodes_df$x <- 0
          nodes_df$y <- 0
          for (i in 1:3) {
            idx <- which(nodes_df$group == types[i])
            if (length(idx) > 0) {
              if (types[i] == "Microbe") {
                fams <- unique(nodes_df$family[idx])
                y_offs <- seq(-500, 500, length.out = length(fams))
                for (j in seq_along(fams)) {
                  f_idx <- which(nodes_df$id %in% nodes_df$id[idx][nodes_df$family[idx] == fams[j]])
                  ang <- seq(0, 2 * pi, length.out = length(f_idx) + 1)[1:length(f_idx)]
                  rad <- 90 + (length(f_idx) * 4)
                  nodes_df$x[f_idx] <- x_zones[i] + rad * cos(ang)
                  nodes_df$y[f_idx] <- y_offs[j] + rad * sin(ang)
                }
              } else {
                ang <- seq(0, 2 * pi, length.out = length(idx) + 1)[1:length(idx)]
                rad <- 250 + (length(idx) * 2.5)
                nodes_df$x[idx] <- x_zones[i] + rad * cos(ang)
                nodes_df$y[idx] <- rad * sin(ang)
              }
            }
          }

          # Styling with User Colors
          nodes_df$shape <- c("Microbe" = "hexagon", "Pathway" = "dot", "Metabolite" = "diamond")[nodes_df$group]
          nodes_df$color <- c("Microbe" = "#9AA374", "Pathway" = "#C1ABAD", "Metabolite" = "#4E7286")[nodes_df$group]
          edges_df <- data.frame(from = edges$Feature1, to = edges$Feature2, value = edges$edge_score)

          # Plot Construction
          vis_plot <- visNetwork(nodes_df, edges_df,
            width = "100%", height = "100vh",
            main = list(
              text = "💡 Tip: Scroll to zoom. Drag nodes to perfect layout.",
              style = "color:#475569; font-size:14px; font-family:sans-serif;"
            )
          ) |>
            visNodes(font = list(color = "#1e293b", size = 22, background = "rgba(255,255,255,0.75)"), borderWidth = 1.5) |>
            visEdges(smooth = list(enabled = TRUE, type = "continuous"), color = list(color = "rgba(148, 163, 184, 0.35)", highlight = "#e11d48")) |>
            visGroups(groupname = "Microbe", color = "#9AA374", shape = "hexagon") |>
            visGroups(groupname = "Pathway", color = "#C1ABAD", shape = "dot") |>
            visGroups(groupname = "Metabolite", color = "#4E7286", shape = "diamond") |>
            visLegend(useGroups = TRUE, position = "right", width = 0.1, main = NULL) |>
            visInteraction(navigationButtons = FALSE, dragNodes = TRUE, multiselect = TRUE, hover = TRUE) |>
            visPhysics(
              enabled = TRUE, solver = "forceAtlas2Based",
              forceAtlas2Based = list(avoidOverlap = 1, springConstant = 0.08, centralGravity = 0.005)
            ) |>
            visExport(
              type = "png", label = "💾 SAVE NETWORK",
              style = "position:absolute; right:30px; top:30px; background:#f1f5f9; color:#475569; padding:12px 24px; border-radius:10px; border:1px solid #cbd5e1; cursor:pointer; font-weight:bold; font-family:sans-serif; z-index:1000;"
            )

          # Save HTML
          saveWidget(vis_plot, file = out_html, selfcontained = TRUE)
          message("✅ HTML saved to: ", out_html)
        },
        error = function(e) warning("Visualization failed: ", e$message)
      )
    }
  }
  return(ifelse(file.exists(out_html), out_html, out_csv))
}
