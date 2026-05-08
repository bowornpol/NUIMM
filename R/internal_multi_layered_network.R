#' High-Tech Multi-Omics Network Visualization Helper (Fixed Layout)
#' @keywords internal
utils::globalVariables(c("type", "weight", "name", "edge_score", "layer"))

# --- Helper to Clean Taxonomic Name ---
clean_taxonomy <- function(tax_string) {
  if (is.na(tax_string) || tax_string == "") {
    return(tax_string)
  }

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

          nodes_df <- data.frame(
            id = V(g)$name,
            label = sapply(V(g)$name, function(x) if (V(g)$group[V(g)$name == x] == "Microbe") clean_taxonomy(x) else x),
            group = V(g)$group,
            size = c("Microbe" = 20, "Pathway" = 30, "Metabolite" = 40)[V(g)$group],
            title = paste0("<div style='padding:10px; font-family:sans-serif;'><b>ID:</b> ", V(g)$name, "</div>"),
            stringsAsFactors = FALSE
          )

          # STRICT MATHEMATICAL PLACEMENT
          nodes_df$x <- 0
          nodes_df$y <- 0

          for (grp in c("Microbe", "Pathway", "Metabolite")) {
            idx <- which(nodes_df$group == grp)
            if (length(idx) == 0) next

            if (grp == "Microbe") {
              # ALL MICROBES IN ONE GIANT CIRCLE AT THE TOP
              n_mic <- length(idx)
              ang <- seq(0, 2 * pi, length.out = n_mic + 1)[1:n_mic]
              rad <- 400 + (n_mic * 15) # Automatically scales up so nodes never overlap
              nodes_df$x[idx] <- rad * cos(ang)
              nodes_df$y[idx] <- -800 + rad * sin(ang)
            } else if (grp == "Pathway") {
              # PATHWAYS IN A ROW IN THE MIDDLE
              n_grp <- length(idx)
              x_offs <- seq(-1200, 1200, length.out = max(1, n_grp))
              nodes_df$x[idx] <- x_offs
              nodes_df$y[idx] <- 0
            } else if (grp == "Metabolite") {
              # METABOLITES IN A ROW AT THE BOTTOM
              n_grp <- length(idx)
              x_offs <- seq(-600, 600, length.out = max(1, n_grp))
              nodes_df$x[idx] <- x_offs
              nodes_df$y[idx] <- 600
            }
          }

          # Shapes and Colors
          nodes_df$shape <- c("Microbe" = "triangle", "Pathway" = "dot", "Metabolite" = "square")[nodes_df$group]
          nodes_df$color <- c("Microbe" = "#9AA374", "Pathway" = "#C1ABAD", "Metabolite" = "#4E7286")[nodes_df$group]

          # Render Network
          vis_plot <- visNetwork(nodes_df, edges, width = "100%", height = "100vh") |>
            # font background highlights text so edges don't ruin readability
            visNodes(font = list(color = "#0f172a", size = 18, face = "sans-serif", background = "rgba(255,255,255,0.85)"), borderWidth = 1.5, shadow = TRUE) |>
            # Straight lines
            visEdges(smooth = FALSE, color = list(color = "rgba(180, 180, 180, 0.4)", highlight = "#e11d48"), width = 1) |>
            visGroups(groupname = "Microbe", color = list(background = "#9AA374", border = "#7A825C", highlight = "#B4BE89"), shape = "triangle") |>
            visGroups(groupname = "Pathway", color = list(background = "#C1ABAD", border = "#9A898A", highlight = "#D8C5C7"), shape = "dot") |>
            visGroups(groupname = "Metabolite", color = list(background = "#4E7286", border = "#3A5565", highlight = "#6392AB"), shape = "square") |>
            visLegend(useGroups = TRUE, position = "right", width = 0.1) |>
            visInteraction(navigationButtons = FALSE, dragNodes = TRUE, multiselect = TRUE, hover = TRUE) |>
            # No Physics
            visPhysics(enabled = FALSE) |>
            # BULLETPROOF LIGHT GREY SAVE BUTTON
            visExport(
              type = "png", label = "Save Network",
              style = "position: absolute; right: 20px; top: 20px; background-color: #e2e8f0; color: #334155; padding: 10px 20px; border: 1px solid #cbd5e1; border-radius: 6px; cursor: pointer; font-weight: bold; font-family: sans-serif; z-index: 9999 !important; box-shadow: 0px 4px 6px rgba(0,0,0,0.1);"
            )

          vis_plot$x$background <- "#ffffff"
          saveWidget(vis_plot, file = out_html, selfcontained = TRUE)
          message("✅ Top-Down HTML successfully built: ", out_html)
        },
        error = function(e) {
          message("❌ Visualization failed: ", e$message)
        }
      )
    }
  }

  return(ifelse(file.exists(out_html), out_html, out_csv))
}
