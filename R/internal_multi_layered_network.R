#' High-Tech Multi-Omics Network Visualization Helper (Top-Down Clustered Layout)
#' @keywords internal
utils::globalVariables(c("type", "weight", "name", "edge_score", "layer"))

# --- Robust Taxonomy Extraction ---
clean_and_get_family <- function(tax_string) {
  result <- list(clean = tax_string, family = "Unassigned")
  if (is.na(tax_string) || tax_string == "") {
    return(result)
  }

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

  base_name <- tools::file_path_sans_ext(basename(gsea_file))
  out_csv <- file.path(output_dir, paste0("final_mln_data_", base_name, ".csv"))
  out_html <- file.path(output_dir, paste0("interactive_mln_", base_name, ".html"))

  # Read Data
  mpn <- read_input_file(mpn_file, file_type = "csv", stringsAsFactors = FALSE)
  ppn <- if (!is.na(ppn_file) && file.exists(ppn_file)) read_input_file(ppn_file, file_type = "csv") else NULL
  pmn <- if (!is.null(pmn_file) && file.exists(pmn_file)) read_input_file(pmn_file, file_type = "csv", stringsAsFactors = FALSE) else NULL
  gsea <- read_input_file(gsea_file, file_type = "csv")

  valid_paths <- gsea$ID
  edges <- data.frame()

  if (nrow(mpn) > 0) {
    mpn_valid <- mpn[mpn$FunctionID %in% valid_paths, ]
    if (nrow(mpn_valid) > 0) edges <- rbind(edges, data.frame(from = mpn_valid$TaxonID, to = mpn_valid$FunctionID, value = mpn_valid$relative_contribution, type = "Microbe-Pathway"))
  }
  if (!is.null(ppn)) {
    ppn_valid <- ppn[ppn$FunctionID_1 %in% valid_paths & ppn$FunctionID_2 %in% valid_paths, ]
    if (nrow(ppn_valid) > 0) edges <- rbind(edges, data.frame(from = ppn_valid$FunctionID_1, to = ppn_valid$FunctionID_2, value = ppn_valid$jaccard_index, type = "Pathway-Pathway"))
  }
  if (!is.null(pmn)) {
    pmn_valid <- pmn[pmn$FunctionID %in% valid_paths, ]
    if (nrow(pmn_valid) > 0) edges <- rbind(edges, data.frame(from = pmn_valid$FunctionID, to = pmn_valid$MetaboliteID, value = abs(pmn_valid$correlation), type = "Pathway-Metabolite"))
  }

  if (nrow(edges) > 0) {
    write.csv(edges, out_csv, row.names = FALSE)

    if (visualize) {
      tryCatch(
        {
          library(igraph)
          library(visNetwork)
          library(htmlwidgets)

          g <- graph_from_data_frame(edges, directed = FALSE)
          V(g)$group <- ifelse(V(g)$name %in% edges$from[edges$type == "Microbe-Pathway"], "Microbe",
            ifelse(V(g)$name %in% edges$to[edges$type == "Pathway-Metabolite"], "Metabolite", "Pathway")
          )

          tax_info <- lapply(V(g)$name, function(x) if (V(g)$group[V(g)$name == x] == "Microbe") clean_and_get_family(x) else list(clean = x, family = "NA"))

          nodes_df <- data.frame(
            id = V(g)$name,
            label = sapply(tax_info, function(x) x$clean),
            family = sapply(tax_info, function(x) x$family),
            group = V(g)$group,
            # Scaling sizes: Microbes smaller, Metabolites larger (like the image)
            size = c("Microbe" = 20, "Pathway" = 30, "Metabolite" = 40)[V(g)$group],
            title = paste0("<div style='padding:10px; font-family:sans-serif;'><b>ID:</b> ", V(g)$name, "</div>"),
            stringsAsFactors = FALSE
          )

          # STRICT MATHEMATICAL PLACEMENT: TOP-DOWN
          y_zones <- c("Microbe" = -800, "Pathway" = 0, "Metabolite" = 600)
          nodes_df$x <- 0
          nodes_df$y <- 0

          for (grp in c("Microbe", "Pathway", "Metabolite")) {
            idx <- which(nodes_df$group == grp)
            if (length(idx) == 0) next

            if (grp == "Microbe") {
              # Top Layer: Horizontal Row of Family Circles
              fams <- unique(nodes_df$family[idx])
              x_offs <- seq(-1400, 1400, length.out = max(1, length(fams)))
              for (j in seq_along(fams)) {
                f_idx <- which(nodes_df$id %in% nodes_df$id[idx][nodes_df$family[idx] == fams[j]])
                n_fam <- length(f_idx)
                ang <- seq(0, 2 * pi, length.out = n_fam + 1)[1:n_fam]
                rad <- 80 + (n_fam * 8) # Circle size depends on node count
                nodes_df$x[f_idx] <- x_offs[j] + rad * cos(ang)
                nodes_df$y[f_idx] <- y_zones[grp] + rad * sin(ang)
              }
            } else if (grp == "Pathway") {
              # Middle Layer: Horizontal Spread
              n_grp <- length(idx)
              x_offs <- seq(-1000, 1000, length.out = max(1, n_grp))
              nodes_df$x[idx] <- x_offs
              # Small random Y jitter so labels don't perfectly overlap
              nodes_df$y[idx] <- y_zones[grp] + runif(n_grp, -50, 50)
            } else if (grp == "Metabolite") {
              # Bottom Layer: Horizontal Spread (tightly grouped like the image)
              n_grp <- length(idx)
              x_offs <- seq(-400, 400, length.out = max(1, n_grp))
              nodes_df$x[idx] <- x_offs
              nodes_df$y[idx] <- y_zones[grp]
            }
          }

          # MATCHING SHAPES TO REFERENCE IMAGE
          nodes_df$shape <- c("Microbe" = "triangle", "Pathway" = "dot", "Metabolite" = "square")[nodes_df$group]
          nodes_df$color <- c("Microbe" = "#9AA374", "Pathway" = "#C1ABAD", "Metabolite" = "#4E7286")[nodes_df$group]

          # Render
          vis_plot <- visNetwork(nodes_df, edges,
            width = "100%", height = "100vh",
            main = list(text = "Hierarchical Multi-Omics Architecture", style = "font-family:sans-serif; color:#1e293b; font-weight:bold; font-size:20px;")
          ) |>
            # Nodes styling (labels below node)
            visNodes(font = list(color = "#0f172a", size = 18, face = "sans-serif", background = "rgba(255,255,255,0.8)", vadjust = -40), borderWidth = 1.5, shadow = TRUE) |>
            # EDGES: Straight lines (smooth = FALSE) to match paper aesthetic
            visEdges(smooth = FALSE, color = list(color = "rgba(180, 180, 180, 0.3)", highlight = "#e11d48"), width = 1) |>
            visGroups(groupname = "Microbe", color = list(background = "#9AA374", border = "#7A825C", highlight = "#B4BE89"), shape = "triangle") |>
            visGroups(groupname = "Pathway", color = list(background = "#C1ABAD", border = "#9A898A", highlight = "#D8C5C7"), shape = "dot") |>
            visGroups(groupname = "Metabolite", color = list(background = "#4E7286", border = "#3A5565", highlight = "#6392AB"), shape = "square") |>
            visLegend(useGroups = TRUE, position = "right", width = 0.1) |>
            visInteraction(navigationButtons = FALSE, dragNodes = TRUE, multiselect = TRUE, hover = TRUE) |>
            # PHYSICS KILLED: Pure Mathematical Layout
            visPhysics(enabled = FALSE) |>
            visExport(
              type = "png", label = "📸 CAPTURE LAYOUT",
              style = "position: fixed; right: 30px; top: 30px; z-index: 999999; background: #f1f5f9; color: #334155; padding: 12px 24px; border-radius: 8px; border: 2px solid #cbd5e1; cursor: pointer; font-weight: bold; font-family: sans-serif;"
            )

          vis_plot$x$background <- "#ffffff"
          saveWidget(vis_plot, file = out_html, selfcontained = TRUE)
          message("✅ Top-Down HTML successfully built without physics: ", out_html)
        },
        error = function(e) {
          message("❌ Visualization failed: ", e$message)
        }
      )
    }
  }

  return(ifelse(file.exists(out_html), out_html, out_csv))
}
