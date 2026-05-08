#' High-Tech Multi-Omics Network Visualization Helper
#' @keywords internal
utils::globalVariables(c("type", "weight", "name", "edge_score", "layer"))

# --- Helper to Clean Name and Extract Family (f__) ---
clean_and_get_family <- function(tax_string) {
  # Default
  result <- list(clean = tax_string, family = "Unassigned")

  if (grepl("g__", tax_string)) {
    # Split by semicolon (character is immune to scoping bugs)
    parts <- unlist(strsplit(tax_string, ";\\s*|;"))

    # 1. Clean Name (Genus Species)
    g_part <- parts[grepl("^g__", parts)][1]
    s_part <- parts[grepl("^s__", parts)][1]
    if (!is.na(s_part) && nchar(s_part) > 3) {
      result$clean <- paste(g_part, s_part, sep = " ")
    } else if (!is.na(g_part)) {
      result$clean <- g_part
    }

    # 2. Extract Family (f__)
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

  # Ensure path is defined early to prevent 'object not found' errors
  base_name <- tools::file_path_sans_ext(basename(gsea_file))
  out_path <- file.path(output_dir, paste0("final_mln_data_", base_name, ".csv"))

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

  if (nrow(edges) > 0) {
    write.csv(edges, out_path, row.names = FALSE)
    library(igraph)
    library(visNetwork)
    library(htmlwidgets)
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% edges$Feature1[edges$edge_type == "Microbe-Pathway"], "Microbe",
      ifelse(igraph::V(g)$name %in% edges$Feature2[edges$edge_type == "Pathway-Metabolite"], "Metabolite", "Pathway")
    )

    if (visualize) {
      tryCatch(
        {
          # 1. Advanced Node Prep (Clean Name + Family Grouping)
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
            title = paste0("<div style='padding:8px; font-family:sans-serif;'><b>ID:</b> ", igraph::V(g)$name, "<br><b>Family:</b> ", unname(sapply(tax_info, function(x) x$family)), "</div>"),
            stringsAsFactors = FALSE
          )

          # 2. Strict 3-Zone Hierarchy Logic (Left -> Middle -> Right)
          types <- c("Microbe", "Pathway", "Metabolite")
          x_zones <- c(-800, 0, 800) # Left Zone, Center Zone, Right Zone
          nodes_df$x <- 0
          nodes_df$y <- 0

          for (i in 1:3) {
            idx <- which(nodes_df$group == types[i])
            if (length(idx) > 0) {
              # Unique logic for Microbes: Multiple Family Rings
              if (types[i] == "Microbe") {
                families <- unique(nodes_df$family[idx])
                fam_y_offsets <- seq(-500, 500, length.out = length(families))
                for (j in seq_along(families)) {
                  fam_idx_in_nodes <- which(nodes_df$id %in% nodes_df$id[idx][nodes_df$family[idx] == families[j]])
                  n_fam <- length(fam_idx_in_nodes)
                  # Distribute in a small ring
                  angle <- seq(0, 2 * pi, length.out = n_fam + 1)[1:n_fam]
                  radius <- 90 + (n_fam * 3) # Wide rings help 'กระจาย'
                  nodes_df$x[fam_idx_in_nodes] <- x_zones[i] + radius * cos(angle)
                  nodes_df$y[fam_idx_in_nodes] <- fam_y_offsets[j] + radius * sin(angle)
                }
              } else {
                # General Circle for Pathway & Metabolite Zones
                n_grp <- length(idx)
                angle <- seq(0, 2 * pi, length.out = n_grp + 1)[1:n_grp]
                radius <- 200 + (n_grp * 2)
                nodes_df$x[idx] <- x_zones[i] + radius * cos(angle)
                nodes_df$y[idx] <- radius * sin(angle)
              }
            }
          }

          # Cinematic Aesthetic Mapping
          nodes_df$shape <- c("Microbe" = "hexagon", "Pathway" = "dot", "Metabolite" = "diamond")[nodes_df$group]
          nodes_df$color <- c("Microbe" = "#0ea5e9", "Pathway" = "#8b5cf6", "Metabolite" = "#f59e0b")[nodes_df$group]
          nodes_df$shadow <- list(enabled = TRUE, color = "rgba(0,0,0,0.1)", size = 10, x = 3, y = 3)

          edges_df <- data.frame(from = edges$Feature1, to = edges$Feature2, value = edges$edge_score)

          # 3. High-Tech UI Implementation
          vis_plot <- visNetwork::visNetwork(
            nodes_df, edges_df,
            width = "100%", height = "100vh", # Locked to full screen
            main = list(
              text = "💡 Tip: Scroll to zoom. Use the console on the right to change colors and physics.",
              style = "color: #334155; font-family: sans-serif; font-size: 16px; font-weight: bold; text-align: center;"
            )
          ) |>
            visNetwork::visNodes(
              font = list(color = "#1e293b", size = 20, face = "sans-serif", background = "rgba(255,255,255,0.7)"), # Stroke background for text
              borderWidth = 1.5,
              borderWidthSelected = 5,
              scaling = list(min = 20, max = 50)
            ) |>
            visNetwork::visEdges(
              smooth = list(enabled = TRUE, type = "continuous"),
              color = list(color = "rgba(148, 163, 184, 0.4)", highlight = "#e11d48"),
              width = 1.5,
              selectionWidth = 3
            ) |>
            # Pre-define Groups for Legend
            visNetwork::visGroups(groupname = "Microbe", color = "#0ea5e9", shape = "hexagon") |>
            visNetwork::visGroups(groupname = "Pathway", color = "#8b5cf6", shape = "dot") |>
            visNetwork::visGroups(groupname = "Metabolite", color = "#f59e0b", shape = "diamond") |>
            visNetwork::visLegend(useGroups = TRUE, position = "right", width = 0.1, main = NULL) |> # Untitled right-side legend
            visNetwork::visInteraction(
              navigationButtons = FALSE, # Removes messy bottom cursors that 'fall off'
              dragNodes = TRUE, dragView = TRUE, zoomView = TRUE, multiselect = TRUE, hover = TRUE
            ) |>
            # FORCEATLAS2 ENGINE FOR PRACTICAL NO-OVERLAP DISTRIBUTIONS
            visNetwork::visPhysics(
              enabled = TRUE,
              solver = "forceAtlas2Based",
              forceAtlas2Based = list(
                gravitationalConstant = -250,
                centralGravity = 0.005,
                springLength = 120,
                springConstant = 0.08,
                avoidOverlap = 1 # CRITICAL: Actively pushes nodes apart
              ),
              stabilization = list(enabled = TRUE, iterations = 200)
            ) |>
            visNetwork::visExport(
              type = "png", label = "💾 SAVE NETWORK",
              # Custom CSS to pin the button to the ขวาบน corner
              style = "position:absolute; right:30px; top:30px; background:#f1f5f9; color:#0f172a; padding:12px 24px; border-radius:8px; border:2px solid #cbd5e1; cursor:pointer; font-weight:bold; font-family:sans-serif; z-index:1000; box-shadow: 0 4px 6px rgba(0,0,0,0.1);"
            ) |>
            # ADDS THE USER COLOR/PHYSICS CONTROL PANEL ON THE RIGHT:
            visNetwork::visConfigure(enabled = TRUE, filter = c("nodes", "edges", "physics"), showButton = FALSE)

          vis_plot$x$background <- "#ffffff" # Pure white high-tech dash

          html_path <- file.path(output_dir, paste0("interactive_mln_", base_name, ".html"))
          htmlwidgets::saveWidget(vis_plot, file = html_path, selfcontained = TRUE, background = "#ffffff")
        },
        error = function(e) warning("Visualization failed: ", e$message)
      )
    }
  }
  return(out_path)
}
