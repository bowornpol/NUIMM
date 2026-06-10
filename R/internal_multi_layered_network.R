utils::globalVariables(c("type", "weight", "name", "edge_score", "layer"))

# Helper to clean taxonomic name
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

#' Internal Multi-Layered Network Assembly and Visualization
#'
#' @details
#' Merges the three sub-network layers (MPN, PPN, PMN) into a unified edge list,
#' writes the final CSV, and renders an interactive HTML network using visNetwork
#' with a circular layout and customization panel.
#'
#' @param gsea_file Path to GSEA results CSV.
#' @param mpn_file Path to Microbe-Pathway network CSV.
#' @param ppn_file Path to Pathway-Pathway Jaccard network CSV.
#' @param pmn_file Path to Pathway-Metabolite correlation CSV.
#' @param output_dir Path to output directory.
#' @param visualize Logical. If TRUE, generates interactive HTML.
#' @param layout_method Network layout algorithm.
#' @param node_colors Named character vector of group colors.
#' @param node_shapes Named character vector of group shapes.
#' @param base_node_size Base size for network nodes.
#' @param plot_width Figure width in inches.
#' @param plot_height Figure height in inches.
#' @param plot_dpi Output image resolution (DPI).
#' @return Path to the output HTML or CSV file.
#' @keywords internal
con_mln_int <- function(
  gsea_file, mpn_file, ppn_file, pmn_file, output_dir,
  visualize, layout_method, node_colors, node_shapes, base_node_size, plot_width, plot_height, plot_dpi,
  ppn_map_database = c("kegg", "metacyc", "custom"), map_file = NULL
) {
  ppn_map_database <- match.arg(ppn_map_database)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  base_name <- tools::file_path_sans_ext(basename(gsea_file))
  out_csv <- file.path(output_dir, paste0("final_mln_", base_name, ".csv"))
  out_html <- file.path(output_dir, paste0("interactive_mln_", base_name, ".html"))

  # Safe file check
  is_valid_file <- function(f) length(f) == 1 && !is.na(f) && file.exists(f)

  # Read data
  mpn <- if (is_valid_file(mpn_file)) read_input_file(mpn_file, file_type = "csv", stringsAsFactors = FALSE) else data.frame()
  ppn <- if (is_valid_file(ppn_file)) read_input_file(ppn_file, file_type = "csv") else NULL
  pmn <- if (is_valid_file(pmn_file)) read_input_file(pmn_file, file_type = "csv", stringsAsFactors = FALSE) else NULL
  gsea <- read_input_file(gsea_file, file_type = "csv")

  valid_paths <- gsea$ID
  edges <- data.frame()

  if (nrow(mpn) > 0) {
    mpn_valid <- mpn[mpn$FunctionID %in% valid_paths, ]
    if (nrow(mpn_valid) > 0) edges <- rbind(edges, data.frame(from = mpn_valid$TaxonID, to = mpn_valid$FunctionID, value = mpn_valid$relative_contribution, type = "Microbe-Pathway"))
  }
  if (!is.null(ppn) && nrow(ppn) > 0) {
    ppn_valid <- ppn[ppn$FunctionID_1 %in% valid_paths & ppn$FunctionID_2 %in% valid_paths, , drop = FALSE]
    if (nrow(ppn_valid) > 0) edges <- rbind(edges, data.frame(from = ppn_valid$FunctionID_1, to = ppn_valid$FunctionID_2, value = ppn_valid$jaccard_index, type = "Pathway-Pathway"))
  }
  if (!is.null(pmn)) {
    pmn_valid <- pmn[pmn$FunctionID %in% valid_paths, ]
    if (nrow(pmn_valid) > 0) edges <- rbind(edges, data.frame(from = pmn_valid$FunctionID, to = pmn_valid$MetaboliteID, value = abs(pmn_valid$correlation), type = "Pathway-Metabolite"))
  }

  # Load pathway ID to name mapping
  pwy_names <- list()
  if (ppn_map_database == "kegg") {
    name_path <- system.file("extdata", "KEGG_pwy_name.csv", package = "NUIMM")
    if (file.exists(name_path)) {
      name_df <- read_input_file(name_path, file_type = "csv", header = FALSE, stringsAsFactors = FALSE)
      pwy_names <- setNames(as.character(name_df[[2]]), as.character(name_df[[1]]))
    }
  } else if (ppn_map_database == "metacyc") {
    name_path <- system.file("extdata", "metacyc_pwy_name.csv", package = "NUIMM")
    if (file.exists(name_path)) {
      name_df <- read_input_file(name_path, file_type = "csv", header = FALSE, stringsAsFactors = FALSE)
      pwy_names <- setNames(as.character(name_df[[2]]), as.character(name_df[[1]]))
    }
  } else if (ppn_map_database == "custom") {
    if (!is.null(map_file) && file.exists(map_file)) {
      name_df <- read_input_file(map_file, file_type = "csv", header = FALSE, stringsAsFactors = FALSE)
      if (ncol(name_df) >= 2) {
        pwy_names <- setNames(as.character(name_df[[2]]), as.character(name_df[[1]]))
      }
    }
  }

  if (nrow(edges) > 0 && length(pwy_names) > 0) {
    translate_vec <- function(vec) {
      sapply(vec, function(x) {
        if (x %in% names(pwy_names)) {
          val <- pwy_names[[x]]
          if (!is.na(val) && val != "") {
            return(val)
          }
        }
        return(x)
      })
    }

    idx_mp <- edges$type == "Microbe-Pathway"
    if (any(idx_mp)) {
      edges$to[idx_mp] <- translate_vec(edges$to[idx_mp])
    }

    idx_pp <- edges$type == "Pathway-Pathway"
    if (any(idx_pp)) {
      edges$from[idx_pp] <- translate_vec(edges$from[idx_pp])
      edges$to[idx_pp] <- translate_vec(edges$to[idx_pp])
    }

    idx_pm <- edges$type == "Pathway-Metabolite"
    if (any(idx_pm)) {
      edges$from[idx_pm] <- translate_vec(edges$from[idx_pm])
    }
  }

  if (nrow(edges) > 0) {
    message(sprintf("    Multi-layered network assembled: |V|=%d, |E|=%d.", length(unique(c(edges$from, edges$to))), nrow(edges)))
    write.csv(edges, out_csv, row.names = FALSE)
    edges$title <- paste0("<div style='padding:10px; font-family:sans-serif;'><b>Value:</b> ", round(edges$value, 4), "</div>")

    if (visualize) {
      tryCatch(
        {
          g <- igraph::graph_from_data_frame(edges, directed = FALSE)

          node_groups <- as.character(determine_node_groups(igraph::V(g)$name, edges, "from", "to"))

          # Leave shape and color to be handled strictly by groups
          nodes_df <- data.frame(
            id = igraph::V(g)$name,
            label = sapply(seq_along(igraph::V(g)$name), function(i) {
              if (node_groups[i] == "Microbe") clean_taxonomy(igraph::V(g)$name[i]) else igraph::V(g)$name[i]
            }),
            group = node_groups,
            size = c("Microbe" = 20, "Pathway" = 30, "Metabolite" = 40)[node_groups],
            title = paste0("<div style='padding:10px; font-family:sans-serif;'><b>ID:</b> ", igraph::V(g)$name, "</div>"),
            stringsAsFactors = FALSE
          )

          # Initialize node coordinates
          nodes_df$x <- 0
          nodes_df$y <- 0

          idx_mic <- which(nodes_df$group == "Microbe")
          idx_path <- which(nodes_df$group == "Pathway")
          idx_met <- which(nodes_df$group == "Metabolite")

          r_mic <- 200 + (length(idx_mic) * 15)
          r_path <- 150 + (length(idx_path) * 20)
          r_met <- 100 + (length(idx_met) * 25)

          x_mic <- -(r_mic + r_path + 500)
          x_path <- 0
          x_met <- (r_path + r_met + 500)

          if (length(idx_mic) > 0) {
            ang <- seq(0, 2 * pi, length.out = length(idx_mic) + 1)[1:length(idx_mic)]
            nodes_df$x[idx_mic] <- x_mic + r_mic * cos(ang)
            nodes_df$y[idx_mic] <- r_mic * sin(ang)
          }
          if (length(idx_path) > 0) {
            ang <- seq(0, 2 * pi, length.out = length(idx_path) + 1)[1:length(idx_path)]
            nodes_df$x[idx_path] <- x_path + r_path * cos(ang)
            nodes_df$y[idx_path] <- r_path * sin(ang)
          }
          if (length(idx_met) > 0) {
            ang <- seq(0, 2 * pi, length.out = length(idx_met) + 1)[1:length(idx_met)]
            nodes_df$x[idx_met] <- x_met + r_met * cos(ang)
            nodes_df$y[idx_met] <- r_met * sin(ang)
          }

          # Canvas legend configuration
          max_y <- max(nodes_df$y, na.rm = TRUE)
          legend_y <- max_y + 400

          legend_nodes <- data.frame(
            id = c("LEG_MIC", "LEG_PATH", "LEG_MET"),
            label = c("Microbe", "Pathway", "Metabolite"),
            group = c("Microbe", "Pathway", "Metabolite"),
            size = c(60, 60, 60),
            title = c("", "", ""),
            x = c(-300, 0, 300),
            y = c(legend_y, legend_y, legend_y),
            stringsAsFactors = FALSE
          )
          nodes_df <- rbind(nodes_df, legend_nodes)

          # JavaScript configuration for interactive UI
          js_custom_panel <- paste0("
          function(el, x, data) {
            var wrapper = document.createElement('div');
            wrapper.style.position = 'absolute';
            wrapper.style.left = '20px';
            wrapper.style.bottom = '20px';
            wrapper.style.zIndex = '99999';
            wrapper.style.fontFamily = 'sans-serif';

            var toggleBtn = document.createElement('button');
            toggleBtn.innerHTML = 'Hide Controls';
            toggleBtn.style.padding = '8px 16px';
            toggleBtn.style.backgroundColor = '#f1f5f9';
            toggleBtn.style.color = '#0f172a';
            toggleBtn.style.border = '1px solid #cbd5e1';
            toggleBtn.style.borderRadius = '6px';
            toggleBtn.style.cursor = 'pointer';
            toggleBtn.style.fontWeight = 'bold';
            toggleBtn.style.boxShadow = '0 4px 6px rgba(0,0,0,0.1)';
            toggleBtn.style.display = 'block';
            toggleBtn.style.marginBottom = '10px';
            wrapper.appendChild(toggleBtn);

            var panel = document.createElement('div');
            panel.style.backgroundColor = 'rgba(255, 255, 255, 0.95)';
            panel.style.padding = '15px';
            panel.style.border = '1px solid #cbd5e1';
            panel.style.borderRadius = '8px';
            panel.style.boxShadow = '0 4px 6px rgba(0,0,0,0.1)';

            var isPanelOpen = true;
            toggleBtn.onclick = function() {
              isPanelOpen = !isPanelOpen;
              if (isPanelOpen) {
                panel.style.display = 'block';
                toggleBtn.innerHTML = 'Hide Controls';
              } else {
                panel.style.display = 'none';
                toggleBtn.innerHTML = 'Customize Network';
              }
            };

            var tip = document.createElement('div');
            tip.innerHTML = '<b>Tip:</b> Scroll to zoom. Hold Ctrl + drag to select multiple nodes. Drag nodes to perfect layout.<br><hr style=\"margin:10px 0; border:0; border-top:1px solid #e2e8f0;\">';
            tip.style.fontSize = '14px';
            tip.style.color = '#475569';
            tip.style.marginBottom = '10px';
            panel.appendChild(tip);

            var title = document.createElement('div');
            title.innerHTML = '<b>Customize Network</b>';
            title.style.fontSize = '14px';
            title.style.color = '#0f172a';
            title.style.marginBottom = '8px';
            panel.appendChild(title);

            var groups = [
              {name: 'Microbe', color: '#9aa374', shape: 'triangle'},
              {name: 'Pathway', color: '#c1abad', shape: 'dot'},
              {name: 'Metabolite', color: '#4e7286', shape: 'square'}
            ];

            // Target the visualization engine instance
            var widget = this;
            var visEngine = widget.network;

            groups.forEach(function(g) {
              var wrap = document.createElement('div');
              wrap.style.display = 'flex';
              wrap.style.alignItems = 'center';
              wrap.style.marginBottom = '6px';
              wrap.style.fontSize = '13px';

              var lbl = document.createElement('div');
              lbl.innerText = g.name + ':';
              lbl.style.width = '75px';
              lbl.style.fontWeight = 'bold';
              wrap.appendChild(lbl);

              var colorIn = document.createElement('input');
              colorIn.type = 'color';
              colorIn.value = g.color;
              colorIn.style.width = '25px';
              colorIn.style.height = '25px';
              colorIn.style.padding = '0';
              colorIn.style.border = 'none';
              colorIn.style.cursor = 'pointer';
              wrap.appendChild(colorIn);

              var shapeSel = document.createElement('select');
              var shapes = ['dot', 'triangle', 'square', 'diamond', 'star', 'hexagon'];
              shapes.forEach(function(s) {
                var opt = document.createElement('option');
                opt.value = s;
                opt.text = s;
                if(s === g.shape) opt.selected = true;
                shapeSel.appendChild(opt);
              });
              shapeSel.style.marginLeft = '10px';
              shapeSel.style.padding = '2px';
              shapeSel.style.borderRadius = '4px';
              shapeSel.style.border = '1px solid #cbd5e1';
              wrap.appendChild(shapeSel);

              panel.appendChild(wrap);

              // Update the group styling directly in the engine
              var updateGraph = function() {
                var newColor = colorIn.value;
                var newShape = shapeSel.value;

                var options = { groups: {} };
                options.groups[g.name] = {
                  color: { background: newColor, border: '#475569', highlight: newColor },
                  shape: newShape
                };

                if (visEngine && typeof visEngine.setOptions === 'function') {
                  visEngine.setOptions(options);
                }
              };

              colorIn.addEventListener('change', updateGraph);
              shapeSel.addEventListener('change', updateGraph);
            });

            var saveBtn = document.createElement('button');
            saveBtn.innerHTML = 'Save Network';
            saveBtn.style.marginTop = '15px';
            saveBtn.style.padding = '8px 16px';
            saveBtn.style.backgroundColor = '#f1f5f9';
            saveBtn.style.color = '#0f172a';
            saveBtn.style.border = '1px solid #cbd5e1';
            saveBtn.style.borderRadius = '6px';
            saveBtn.style.cursor = 'pointer';
            saveBtn.style.fontWeight = 'bold';
            saveBtn.style.width = '100%';

            saveBtn.onclick = function() {
              var originalCanvas = el.getElementsByTagName('canvas')[0];
              if (!originalCanvas) return;

              var tempCanvas = document.createElement('canvas');
              tempCanvas.width = originalCanvas.width;
              tempCanvas.height = originalCanvas.height;
              var ctx = tempCanvas.getContext('2d');
              ctx.fillStyle = '#ffffff';
              ctx.fillRect(0, 0, tempCanvas.width, tempCanvas.height);
              ctx.drawImage(originalCanvas, 0, 0);

              var safeName = data ? data : 'network';
              var link = document.createElement('a');
              link.download = 'fig_mln_' + safeName + '.png';
              link.href = tempCanvas.toDataURL('image/png');
              link.click();
            };
            panel.appendChild(saveBtn);

            wrapper.appendChild(panel);
            el.appendChild(wrapper);
          ", get_ctrl_drag_js(), "
          }
          ")

          # Render network
          vis_plot <- visNetwork::visNetwork(nodes_df, edges, width = "100%", height = "95vh") |>
            visNetwork::visNodes(font = list(color = "#0f172a", size = 35, face = "sans-serif", background = "rgba(255,255,255,0.85)"), borderWidth = 1.5, shadow = TRUE) |>
            visNetwork::visEdges(smooth = FALSE, color = list(color = "rgba(180, 180, 180, 0.4)", highlight = "#e11d48"), width = 1) |>
            visNetwork::visGroups(groupname = "Microbe", color = list(background = "#9AA374", border = "#7A825C", highlight = "#B4BE89"), shape = "triangle") |>
            visNetwork::visGroups(groupname = "Pathway", color = list(background = "#C1ABAD", border = "#9A898A", highlight = "#D8C5C7"), shape = "dot") |>
            visNetwork::visGroups(groupname = "Metabolite", color = list(background = "#4E7286", border = "#3A5565", highlight = "#6392AB"), shape = "square") |>
            visNetwork::visInteraction(navigationButtons = FALSE, dragNodes = TRUE, multiselect = TRUE, hover = TRUE) |>
            visNetwork::visPhysics(enabled = FALSE) |>
            htmlwidgets::onRender(js_custom_panel, data = base_name)

          vis_plot$x$background <- "#ffffff"
          save_widget_safe(vis_plot, file = out_html, title = "NUIMM")

        },
        error = function(e) {
          message("Visualization rendering failed: ", e$message)
        }
      )
    }
  }

  return(ifelse(file.exists(out_html), out_html, out_csv))
}
