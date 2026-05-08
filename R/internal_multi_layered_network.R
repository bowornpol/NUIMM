#' High-Tech Multi-Omics Network Visualization Helper (Advanced Edition)
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
  # 1. Validation & Initialization
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  base_name <- tools::file_path_sans_ext(basename(gsea_file))
  out_csv <- file.path(output_dir, paste0("final_mln_data_", base_name, ".csv"))
  out_html <- file.path(output_dir, paste0("interactive_mln_", base_name, ".html"))

  # Read Data Safely
  tryCatch(
    {
      mpn <- read_input_file(mpn_file, file_type = "csv", stringsAsFactors = FALSE)
      ppn <- if (!is.na(ppn_file) && file.exists(ppn_file)) read_input_file(ppn_file, file_type = "csv") else NULL
      pmn <- if (!is.null(pmn_file) && file.exists(pmn_file)) read_input_file(pmn_file, file_type = "csv", stringsAsFactors = FALSE) else NULL
      gsea <- read_input_file(gsea_file, file_type = "csv")
    },
    error = function(e) {
      stop("Data Load Error: Ensure all input files are formatted correctly. Details: ", e$message)
    }
  )

  # 2. Edge Construction
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

  if (nrow(edges) == 0) {
    message("Warning: No valid edges found for ", base_name)
    return(NULL)
  }

  write.csv(edges, out_csv, row.names = FALSE)

  if (visualize) {
    tryCatch(
      {
        library(igraph)
        library(visNetwork)
        library(htmlwidgets)

        # 3. Node Processing
        g <- graph_from_data_frame(edges, directed = FALSE)

        # Strict Group Assignment
        V(g)$group <- ifelse(V(g)$name %in% edges$from[edges$type == "Microbe-Pathway"], "Microbe",
          ifelse(V(g)$name %in% edges$to[edges$type == "Pathway-Metabolite"], "Metabolite", "Pathway")
        )

        tax_info <- lapply(V(g)$name, function(x) if (V(g)$group[V(g)$name == x] == "Microbe") clean_and_get_family(x) else list(clean = x, family = "NA"))

        nodes_df <- data.frame(
          id = V(g)$name,
          label = sapply(tax_info, function(x) x$clean),
          family = sapply(tax_info, function(x) x$family),
          group = V(g)$group,
          title = paste0("<div style='padding:10px; font-family:sans-serif; font-size:14px; background:#f8fafc; border:1px solid #cbd5e1; border-radius:6px;'><b>ID:</b> ", V(g)$name, "</div>"),
          stringsAsFactors = FALSE
        )

        # 4. Dynamic Pre-Layout Algorithm (Scalable 3-Zone Architecture)
        x_zones <- c("Microbe" = -1500, "Pathway" = 0, "Metabolite" = 1500)
        nodes_df$x <- 0
        nodes_df$y <- 0

        for (grp in c("Microbe", "Pathway", "Metabolite")) {
          idx <- which(nodes_df$group == grp)
          if (length(idx) == 0) next

          if (grp == "Microbe") {
            fams <- unique(nodes_df$family[idx])
            y_offs <- seq(-1000, 1000, length.out = max(1, length(fams)))
            for (j in seq_along(fams)) {
              f_idx <- which(nodes_df$id %in% nodes_df$id[idx][nodes_df$family[idx] == fams[j]])
              n_fam <- length(f_idx)
              ang <- seq(0, 2 * pi, length.out = n_fam + 1)[1:n_fam]
              rad <- 150 + (n_fam * 10) # Wider rings to prevent overlap
              nodes_df$x[f_idx] <- x_zones[grp] + rad * cos(ang)
              nodes_df$y[f_idx] <- y_offs[j] + rad * sin(ang)
            }
          } else {
            n_grp <- length(idx)
            ang <- seq(0, 2 * pi, length.out = n_grp + 1)[1:n_grp]
            rad <- 400 + (n_grp * 8)
            nodes_df$x[idx] <- x_zones[grp] + rad * cos(ang)
            nodes_df$y[idx] <- rad * sin(ang)
          }
        }

        # Shape Mapping
        nodes_df$shape <- c("Microbe" = "hexagon", "Pathway" = "dot", "Metabolite" = "diamond")[nodes_df$group]

        # 5. Advanced UI & Rendering
        vis_plot <- visNetwork(nodes_df, edges,
          width = "100%", height = "100vh",
          main = list(text = "Multi-Omics System Architecture", style = "font-family:sans-serif; color:#1e293b; font-weight:bold; font-size:20px;")
        ) |>
          # Node & Edge Design
          visNodes(font = list(color = "#0f172a", size = 26, face = "sans-serif", background = "rgba(255,255,255,0.85)"), borderWidth = 2, shadow = TRUE) |>
          visEdges(smooth = list(enabled = TRUE, type = "continuous"), color = list(color = "rgba(100, 116, 139, 0.25)", highlight = "#e11d48"), width = 1.5) |>
          # Color & Legend System
          visGroups(groupname = "Microbe", color = list(background = "#9AA374", border = "#7A825C", highlight = "#B4BE89"), shape = "hexagon") |>
          visGroups(groupname = "Pathway", color = list(background = "#C1ABAD", border = "#9A898A", highlight = "#D8C5C7"), shape = "dot") |>
          visGroups(groupname = "Metabolite", color = list(background = "#4E7286", border = "#3A5565", highlight = "#6392AB"), shape = "diamond") |>
          visLegend(useGroups = TRUE, position = "right", width = 0.1) |>
          # Interactive Search & Filter Panels (Top Left)
          visOptions(
            highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE, algorithm = "hierarchical"),
            selectedBy = list(variable = "group", style = "width: 200px; padding: 8px; border-radius: 5px; font-family: sans-serif; background: #f8fafc; border: 1px solid #cbd5e1;"),
            nodesIdSelection = list(enabled = TRUE, style = "width: 200px; padding: 8px; border-radius: 5px; font-family: sans-serif; background: #f8fafc; border: 1px solid #cbd5e1;")
          ) |>
          visInteraction(navigationButtons = FALSE, dragNodes = TRUE, multiselect = TRUE, hover = TRUE, tooltipDelay = 100) |>
          # High-Performance Physics (Aggressive separation, then Auto-Freeze)
          visPhysics(
            enabled = TRUE,
            solver = "forceAtlas2Based",
            forceAtlas2Based = list(
              gravitationalConstant = -15000, # Massive repulsion to clear labels
              centralGravity = 0.005,
              springLength = 150,
              springConstant = 0.05,
              avoidOverlap = 1
            ),
            stabilization = list(enabled = TRUE, iterations = 300, fit = TRUE)
          ) |>
          # Custom Javascript Injection: Freezes network after initial layout for performance
          visEvents(type = "on", stabilizationIterationsDone = "function () { this.setOptions( { physics: false } ); }") |>
          # Advanced Customization Sandbox
          visConfigure(enabled = TRUE, filter = c("nodes", "edges", "physics", "layout"), showButton = TRUE) |>
          # BULLETPROOF SAVE BUTTON (position: fixed guarantees it stays in viewport)
          visExport(
            type = "png", label = "💾 SAVE NETWORK",
            style = "position: fixed; right: 30px; top: 30px; z-index: 999999; background: #f1f5f9; color: #334155; padding: 12px 24px; border-radius: 8px; border: 2px solid #cbd5e1; cursor: pointer; font-weight: bold; font-family: sans-serif; box-shadow: 0px 4px 10px rgba(0,0,0,0.15); transition: all 0.2s ease-in-out;"
          )

        # Render & Save
        vis_plot$x$background <- "#ffffff"
        saveWidget(vis_plot, file = out_html, selfcontained = TRUE)
        message("✅ Interactive Architecture successfully built: ", out_html)
      },
      error = function(e) {
        message("❌ Visualization failed during construction: ", e$message)
      }
    )
  }

  return(ifelse(file.exists(out_html), out_html, out_csv))
}
