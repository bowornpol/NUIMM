#' Internal MLN Integration & Interactive Visualization Helper
#' @keywords internal
utils::globalVariables(c("type", "weight", "name", "edge_score", "layer"))

# --- Futuristic Taxonomy Cleaner & Family Extractor ---
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
          tax_info <- lapply(igraph::V(g)$name, function(x) {
            if (igraph::V(g)$type[igraph::V(g)$name == x] == "Microbe") {
              return(clean_and_get_family(x))
            } else {
              return(list(clean = x, family = "NA"))
            }
          })

          nodes_df <- data.frame(
            id = igraph::V(g)$name,
            label = unname(sapply(tax_info, function(x) x$clean)),
            family = unname(sapply(tax_info, function(x) x$family)),
            group = igraph::V(g)$type,
            stringsAsFactors = FALSE
          )

          # Zone positions (Left/Middle/Right)
          types <- c("Microbe", "Pathway", "Metabolite")
          x_zones <- c(-850, 0, 850)
          nodes_df$x <- 0
          nodes_df$y <- 0

          # Logic to prevent 'ซ้อนทับกัน' and distribute by family circles
          for (i in 1:3) {
            idx <- which(nodes_df$group == types[i])
            if (length(idx) > 0) {
              if (types[i] == "Microbe") {
                families <- unique(nodes_df$family[idx])
                fam_y_offsets <- seq(-500, 500, length.out = length(families))
                for (j in seq_along(families)) {
                  fam_idx <- which(nodes_df$id %in% nodes_df$id[idx][nodes_df$family[idx] == families[j]])
                  n_fam <- length(fam_idx)
                  angle <- seq(0, 2 * pi, length.out = n_fam + 1)[1:n_fam]
                  radius <- 90 + (n_fam * 4)
                  nodes_df$x[fam_idx] <- x_zones[i] + radius * cos(angle)
                  nodes_df$y[fam_idx] <- fam_y_offsets[j] + radius * sin(angle)
                }
              } else {
                n_grp <- length(idx)
                angle <- seq(0, 2 * pi, length.out = n_grp + 1)[1:n_grp]
                radius <- 250 + (n_grp * 2.5)
                nodes_df$x[idx] <- x_zones[i] + radius * cos(angle)
                nodes_df$y[idx] <- radius * sin(angle)
              }
            }
          }

          # NEW USER DEFINED COLOR PALETTE
          nodes_df$shape <- c("Microbe" = "hexagon", "Pathway" = "dot", "Metabolite" = "diamond")[nodes_df$group]
          nodes_df$color <- c("Microbe" = "#9AA374", "Pathway" = "#C1ABAD", "Metabolite" = "#4E7286")[nodes_df$group]
          nodes_df$shadow <- list(enabled = TRUE, color = "rgba(0,0,0,0.1)", size = 10)

          edges_df <- data.frame(from = edges$Feature1, to = edges$Feature2, value = edges$edge_score)

          vis_plot <- visNetwork::visNetwork(
            nodes_df, edges_df,
            width = "100%", height = "100vh",
            main = list(
              text = "💡 Tip: Drag clusters to arrange. Use the config panel on the right for more adjustments.",
              style = "color: #475569; font-family: sans-serif; font-size: 14px;"
            )
          ) |>
            visNetwork::visNodes(
              font = list(color = "#1e293b", size = 22, face = "sans-serif", background = "rgba(255,255,255,0.75)"),
              borderWidth = 1.5
            ) |>
            visNetwork::visEdges(
              smooth = list(enabled = TRUE, type = "continuous"),
              color = list(color = "rgba(148, 163, 184, 0.35)", highlight = "#e11d48"),
              width = 1.5
            ) |>
            visNetwork::visGroups(groupname = "Microbe", color = "#9AA374", shape = "hexagon") |>
            visNetwork::visGroups(groupname = "Pathway", color = "#C1ABAD", shape = "dot") |>
            visNetwork::visGroups(groupname = "Metabolite", color = "#4E7286", shape = "diamond") |>
            visNetwork::visLegend(useGroups = TRUE, position = "right", width = 0.1, main = NULL) |>
            visNetwork::visInteraction(
              navigationButtons = FALSE, dragNodes = TRUE, dragView = TRUE, zoomView = TRUE, multiselect = TRUE, hover = TRUE
            ) |>
            visNetwork::visPhysics(
              enabled = TRUE,
              solver = "forceAtlas2Based",
              forceAtlas2Based = list(avoidOverlap = 1, springConstant = 0.08, centralGravity = 0.005),
              stabilization = list(enabled = TRUE, iterations = 200)
            ) |>
            visNetwork::visExport(
              type = "png", label = "💾 SAVE NETWORK",
              # STYLED TOP-RIGHT LIGHT GREY BUTTON
              style = "position:absolute; right:30px; top:30px; background:#f1f5f9; color:#475569; padding:12px 24px; border-radius:10px; border:1px solid #cbd5e1; cursor:pointer; font-weight:bold; font-family:sans-serif; z-index:1000; box-shadow: 0 4px 6px rgba(0,0,0,0.05);"
            ) |>
            visNetwork::visConfigure(enabled = TRUE, filter = c("nodes", "physics"), showButton = FALSE)

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
