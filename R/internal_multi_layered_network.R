utils::globalVariables(c("type", "weight", "name", "edge_score", "layer"))

# Helper to clean taxonomic name
clean_taxonomy <- function(tax_string) {
  if (is.na(tax_string) || tax_string == "") {
    return(tax_string)
  }
  
  parts <- unlist(strsplit(tax_string, ';\\s*|;'))
  s_part <- parts[grepl('^s__', parts)][1]
  g_part <- parts[grepl('^g__', parts)][1]
  f_part <- parts[grepl('^f__', parts)][1]
  o_part <- parts[grepl('^o__', parts)][1]
  c_part <- parts[grepl('^c__', parts)][1]
  p_part <- parts[grepl('^p__', parts)][1]
  d_part <- parts[grepl('^d__', parts)][1]
  
  if (!is.na(s_part) && nchar(s_part) > 3 && !is.na(g_part)) {
    return(paste(g_part, s_part, sep = ' '))
  } else if (!is.na(g_part) && nchar(g_part) > 3) {
    return(g_part)
  } else if (!is.na(f_part) && nchar(f_part) > 3) {
    return(f_part)
  } else if (!is.na(o_part) && nchar(o_part) > 3) {
    return(o_part)
  } else if (!is.na(c_part) && nchar(c_part) > 3) {
    return(c_part)
  } else if (!is.na(p_part) && nchar(p_part) > 3) {
    return(p_part)
  } else if (!is.na(d_part) && nchar(d_part) > 3) {
    return(d_part)
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
#' @noRd
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
    if (nrow(mpn_valid) > 0) edges <- rbind(edges, data.frame(from = mpn_valid$TaxonID, to = mpn_valid$FunctionID, value = mpn_valid$relative_contribution, type = "Microbe-Pathway", direction = NA_character_))
  }
  if (!is.null(ppn) && nrow(ppn) > 0) {
    ppn_valid <- ppn[ppn$FunctionID_1 %in% valid_paths & ppn$FunctionID_2 %in% valid_paths, , drop = FALSE]
    if (nrow(ppn_valid) > 0) {
      ppn_dir <- if ("direction" %in% colnames(ppn_valid)) ppn_valid$direction else ifelse(ppn_valid$jaccard_index > 0, "positive", "negative")
      edges <- rbind(edges, data.frame(from = ppn_valid$FunctionID_1, to = ppn_valid$FunctionID_2, value = ppn_valid$jaccard_index, type = "Pathway-Pathway", direction = ppn_dir))
    }
  }
  if (!is.null(pmn)) {
    pmn_valid <- pmn[pmn$FunctionID %in% valid_paths, ]
    if (nrow(pmn_valid) > 0) {
      # Use direction column if available, otherwise infer from sign
      pmn_dir <- if ("direction" %in% colnames(pmn_valid)) pmn_valid$direction else ifelse(pmn_valid$correlation > 0, "positive", "negative")
      edges <- rbind(edges, data.frame(from = pmn_valid$FunctionID, to = pmn_valid$MetaboliteID, value = pmn_valid$correlation, type = "Pathway-Metabolite", direction = pmn_dir))
    }
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

    # Per-edge color: PMN and PPN edges colored by direction, MPN edges neutral grey
    edge_default_color <- "rgba(160, 160, 160, 0.5)"
    edge_pos_color <- "#B5654A"  # terracotta (positive)
    edge_neg_color <- "#5A7D8B"  # slate blue (negative)
    edges$color <- ifelse(
      edges$type %in% c("Pathway-Metabolite", "Pathway-Pathway"),
      ifelse(!is.na(edges$direction) & edges$direction == "positive", edge_pos_color, edge_neg_color),
      edge_default_color
    )
    edges$title <- paste0("<div style='padding:10px; font-family:sans-serif;'><b>Value:</b> ", round(edges$value, 4), "</div>")
    # Use absolute value for edge width so negative edges render with proper thickness
    edges$width <- abs(edges$value) * 5 + 0.5

    if (visualize) {
      tryCatch(
        {
          g <- igraph::graph_from_data_frame(edges, directed = FALSE)

          node_groups <- as.character(determine_node_groups(igraph::V(g)$name, edges, "from", "to"))

          # Helper to darken a hex color for borders and lighten for highlights
          darken_color <- function(hex, factor = 0.75) {
            rgb_vals <- grDevices::col2rgb(hex)
            darkened <- pmax(0, round(rgb_vals * factor))
            grDevices::rgb(darkened[1], darkened[2], darkened[3], maxColorValue = 255)
          }
          lighten_color <- function(hex, factor = 0.25) {
            rgb_vals <- grDevices::col2rgb(hex)
            lightened <- pmin(255, round(rgb_vals + (255 - rgb_vals) * factor))
            grDevices::rgb(lightened[1], lightened[2], lightened[3], maxColorValue = 255)
          }

          # Build GSEA NES lookup for pathway border coloring
          nes_lookup <- setNames(gsea$NES, gsea$ID)
          # Translate pathway IDs to names if mapping exists
          if (length(pwy_names) > 0) {
            translated_names <- sapply(names(nes_lookup), function(pid) {
              if (pid %in% names(pwy_names) && !is.na(pwy_names[[pid]]) && pwy_names[[pid]] != "") pwy_names[[pid]] else pid
            })
            names(nes_lookup) <- translated_names
          }

          # Assign per-node border color based on GSEA NES direction
          border_pos_color <- "#B5654A"  # terracotta (enriched / positive NES)
          border_neg_color <- "#5A7D8B"  # slate blue (depleted / negative NES)

          # vis.js rule: setting ANY per-node color stops group color inheritance.
          # So we must set background + border + highlight explicitly for all nodes.
          node_bg_colors <- unname(node_colors[node_groups])
          node_hl_colors <- sapply(node_bg_colors, lighten_color)
          node_border_colors <- sapply(seq_along(igraph::V(g)$name), function(i) {
            nid <- igraph::V(g)$name[i]
            if (node_groups[i] == "Pathway" && nid %in% names(nes_lookup)) {
              if (nes_lookup[[nid]] > 0) border_pos_color else border_neg_color
            } else {
              darken_color(unname(node_colors[node_groups[i]]))
            }
          })

          nodes_df <- data.frame(
            id = igraph::V(g)$name,
            label = sapply(seq_along(igraph::V(g)$name), function(i) {
              if (node_groups[i] == "Microbe") clean_taxonomy(igraph::V(g)$name[i]) else igraph::V(g)$name[i]
            }),
            group = node_groups,
            size = c("Microbe" = 20, "Pathway" = 30, "Metabolite" = 40)[node_groups],
            title = paste0("<div style='padding:10px; font-family:sans-serif;'><b>ID:</b> ", igraph::V(g)$name, "</div>"),
            color.background = node_bg_colors,
            color.border = node_border_colors,
            color.highlight = node_hl_colors,
            borderWidth = ifelse(node_groups == "Pathway", 3, 1.5),
            stringsAsFactors = FALSE
          )

          nodes_df <- compute_circular_layout(nodes_df)

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
              {name: 'Microbe', color: '" , tolower(unname(node_colors['Microbe'])),  "', shape: '" , unname(node_shapes['Microbe']),  "'},
              {name: 'Pathway', color: '" , tolower(unname(node_colors['Pathway'])),  "', shape: '" , unname(node_shapes['Pathway']),  "'},
              {name: 'Metabolite', color: '" , tolower(unname(node_colors['Metabolite'])),  "', shape: '" , unname(node_shapes['Metabolite']),  "'}
            ];

            // Target the visualization engine instance
            var widget = this;
            var visEngine = widget.network;
            var legendIcons = {};

            // Helper: apply CSS shape to a legend icon div
            function applyShapeCSS(iconEl, shape, color) {
              iconEl.style.cssText = 'width:14px;height:14px;margin-right:8px;flex-shrink:0;';
              iconEl.style.backgroundColor = color;
              if (shape === 'dot' || shape === 'circle') {
                iconEl.style.borderRadius = '50%';
              } else if (shape === 'diamond') {
                iconEl.style.width = '12px'; iconEl.style.height = '12px';
                iconEl.style.transform = 'rotate(45deg)'; iconEl.style.borderRadius = '2px';
              } else if (shape === 'star') {
                iconEl.style.clipPath = 'polygon(50% 0%,61% 35%,98% 35%,68% 57%,79% 91%,50% 70%,21% 91%,32% 57%,2% 35%,39% 35%)';
              } else if (shape === 'hexagon') {
                iconEl.style.clipPath = 'polygon(25% 0%,75% 0%,100% 50%,75% 100%,25% 100%,0% 50%)';
              } else if (shape === 'triangle') {
                iconEl.style.backgroundColor = 'transparent';
                iconEl.style.width = '0'; iconEl.style.height = '0';
                iconEl.style.borderLeft = '7px solid transparent';
                iconEl.style.borderRight = '7px solid transparent';
                iconEl.style.borderBottom = '14px solid ' + color;
              } else {
                iconEl.style.borderRadius = '2px';
              }
            }

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

              // Update per-node colors via DataSet + sync legend
              var updateGraph = function() {
                var newColor = colorIn.value;
                var newShape = shapeSel.value;

                // Update individual nodes in the DataSet
                if (visEngine && visEngine.body && visEngine.body.data && visEngine.body.data.nodes) {
                  var nodesDS = visEngine.body.data.nodes;
                  var allNodes = nodesDS.get();
                  var updates = [];
                  allNodes.forEach(function(n) {
                    if (n.group === g.name) {
                      updates.push({
                        id: n.id, shape: newShape,
                        color: {
                          background: newColor,
                          border: (n.color && n.color.border) ? n.color.border : newColor,
                          highlight: newColor
                        }
                      });
                    }
                  });
                  if (updates.length > 0) nodesDS.update(updates);
                }

                // Sync legend icon
                if (legendIcons[g.name]) {
                  applyShapeCSS(legendIcons[g.name], newShape, newColor);
                }
              };

              colorIn.addEventListener('change', updateGraph);
              shapeSel.addEventListener('change', updateGraph);
            });

            // --- Scale Layout for Selected Nodes ---
            var scaleWrap = document.createElement('div');
            scaleWrap.style.marginTop = '15px';
            scaleWrap.style.marginBottom = '5px';
            
            var scaleLbl = document.createElement('div');
            scaleLbl.innerHTML = '<b>Adjust Spacing (Selected Nodes):</b> <span id=\"nuimm-scale-val\">1.0x</span>';
            scaleLbl.style.fontSize = '13px';
            scaleLbl.style.marginBottom = '4px';
            scaleWrap.appendChild(scaleLbl);

            var scaleSlider = document.createElement('input');
            scaleSlider.type = 'range';
            scaleSlider.min = '0.5';
            scaleSlider.max = '3.0';
            scaleSlider.step = '0.1';
            scaleSlider.value = '1.0';
            scaleSlider.style.width = '100%';
            scaleWrap.appendChild(scaleSlider);
            panel.appendChild(scaleWrap);

            var originalPositions = {};
            var lastSelectedIds = [];

            function cacheSelectedPositions() {
              lastSelectedIds = visEngine.getSelectedNodes();
              if (lastSelectedIds.length > 0) {
                originalPositions = visEngine.getPositions(lastSelectedIds);
              } else {
                originalPositions = {};
              }
              scaleSlider.value = '1.0';
              document.getElementById('nuimm-scale-val').innerText = '1.0x';
            }

            // Sync cache when selection changes or nodes are manually dragged
            visEngine.on('select', cacheSelectedPositions);
            visEngine.on('deselect', cacheSelectedPositions);
            visEngine.on('selectNode', cacheSelectedPositions);
            visEngine.on('deselectNode', cacheSelectedPositions);
            visEngine.on('dragEnd', function(params) {
              if (params.nodes && params.nodes.length > 0) {
                cacheSelectedPositions();
              }
            });

            function updateScaleVisuals(scale) {
              if (lastSelectedIds.length === 0) return;

              // Compute center of selected nodes
              var cx = 0, cy = 0, count = 0;
              for (var id in originalPositions) {
                cx += originalPositions[id].x;
                cy += originalPositions[id].y;
                count++;
              }
              if (count === 0) return;
              cx /= count;
              cy /= count;

              var updates = [];
              for (var id in originalPositions) {
                var pos = originalPositions[id];
                var dx = pos.x - cx;
                var dy = pos.y - cy;
                var newX = cx + dx * scale;
                var newY = cy + dy * scale;
                
                // Fast interactive mutation
                if (visEngine.body.nodes[id]) {
                  visEngine.body.nodes[id].x = newX;
                  visEngine.body.nodes[id].y = newY;
                  visEngine.body.nodes[id].options.x = newX;
                  visEngine.body.nodes[id].options.y = newY;
                }
                updates.push({ id: id, x: newX, y: newY });
              }
              
              if (typeof visEngine.redraw === 'function') {
                visEngine.redraw();
              }
              return updates;
            }

            scaleSlider.addEventListener('input', function() {
              var scale = parseFloat(this.value);
              document.getElementById('nuimm-scale-val').innerText = scale.toFixed(1) + 'x';
              // Only do visual update on drag (prevents dataset freeze)
              updateScaleVisuals(scale);
            });

            scaleSlider.addEventListener('change', function() {
              var scale = parseFloat(this.value);
              var updates = updateScaleVisuals(scale);
              // Save to dataset ONLY when user releases the slider
              if (updates && updates.length > 0 && visEngine.body && visEngine.body.data && visEngine.body.data.nodes) {
                visEngine.body.data.nodes.update(updates);
              }
            });

            var saveTableBtn = document.createElement('button');
            saveTableBtn.innerHTML = 'Save Network Table (CSV)';
            saveTableBtn.style.marginTop = '15px';
            saveTableBtn.style.padding = '8px 16px';
            saveTableBtn.style.backgroundColor = '#f8fafc';
            saveTableBtn.style.color = '#0f172a';
            saveTableBtn.style.border = '1px solid #cbd5e1';
            saveTableBtn.style.borderRadius = '6px';
            saveTableBtn.style.cursor = 'pointer';
            saveTableBtn.style.fontWeight = 'bold';
            saveTableBtn.style.width = '100%';
            saveTableBtn.style.marginBottom = '4px';
            saveTableBtn.onclick = function() {
              var edgesDS = visEngine.body.data.edges;
              var allEdges = edgesDS.get();
              var csvContent = 'data:text/csv;charset=utf-8,From,To,Value,Type,Direction\\n';
              allEdges.forEach(function(e) {
                 var fromLabel = e.from ? e.from.toString() : '';
                 var toLabel = e.to ? e.to.toString() : '';
                 if (fromLabel.indexOf(',') !== -1) fromLabel = '\"' + fromLabel + '\"';
                 if (toLabel.indexOf(',') !== -1) toLabel = '\"' + toLabel + '\"';
                 var edgeVal = e.value !== undefined ? e.value : '';
                 var edgeType = e.type !== undefined ? e.type : '';
                 var edgeDir = e.direction !== undefined ? e.direction : '';
                 csvContent += fromLabel + ',' + toLabel + ',' + edgeVal + ',' + edgeType + ',' + edgeDir + '\\n';
              });
              var encodedUri = encodeURI(csvContent);
              var link = document.createElement('a');
              link.setAttribute('href', encodedUri);
              link.setAttribute('download', 'network_table.csv');
              document.body.appendChild(link);
              link.click();
              document.body.removeChild(link);
            };
            panel.appendChild(saveTableBtn);

            var saveBtn = document.createElement('button');
            saveBtn.innerHTML = 'Save Network';
            saveBtn.style.marginTop = '4px';
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

              // Draw legend at bottom-right using Canvas 2D API
              var pad = 18;
              var lineH = 24;
              var iconSize = 14;
              var textX = 40;
              var sectionGap = 10;

              // Get current group colors/shapes from the customization inputs
              var currentGroups = [];
              groups.forEach(function(g, idx) {
                var inputs = panel.querySelectorAll('input[type=color]');
                var selects = panel.querySelectorAll('select');
                currentGroups.push({
                  name: g.name,
                  color: inputs[idx] ? inputs[idx].value : g.color,
                  shape: selects[idx] ? selects[idx].value : g.shape
                });
              });

              // Measure legend dimensions
              ctx.font = 'bold 14px sans-serif';
              var titles = ['Node Types', 'Pathway Border', 'Edge Color'];
              var nodeItems = currentGroups.map(function(g) { return g.name; });
              var borderItems = [{label:'Upregulated',color:'#B5654A'},{label:'Downregulated',color:'#5A7D8B'}];
              var edgeItems = [{label:'Positive correlation',color:'#B5654A'},{label:'Negative correlation',color:'#5A7D8B'}];

              var totalLines = 3 + nodeItems.length + borderItems.length + edgeItems.length + 2;
              var legW = 220;
              var legH = pad * 2 + totalLines * lineH + sectionGap * 2;

              var lx = tempCanvas.width - legW - 40;
              var ly = tempCanvas.height - legH - 40;

              // Legend background (no border in saved PNG)
              ctx.fillStyle = '#ffffff';
              ctx.fillRect(lx, ly, legW, legH);

              var curY = ly + pad;

              // Helper: draw shape icon
              function drawShapeIcon(cx, cy, shape, color, size) {
                ctx.fillStyle = color;
                ctx.strokeStyle = color;
                var hs = size / 2;
                if (shape === 'dot' || shape === 'circle') {
                  ctx.beginPath();
                  ctx.arc(cx + hs, cy + hs, hs, 0, 2 * Math.PI);
                  ctx.fill();
                } else if (shape === 'square') {
                  ctx.fillRect(cx, cy, size, size);
                } else if (shape === 'diamond') {
                  ctx.beginPath();
                  ctx.moveTo(cx + hs, cy);
                  ctx.lineTo(cx + size, cy + hs);
                  ctx.lineTo(cx + hs, cy + size);
                  ctx.lineTo(cx, cy + hs);
                  ctx.closePath();
                  ctx.fill();
                } else if (shape === 'triangle') {
                  ctx.beginPath();
                  ctx.moveTo(cx + hs, cy);
                  ctx.lineTo(cx + size, cy + size);
                  ctx.lineTo(cx, cy + size);
                  ctx.closePath();
                  ctx.fill();
                } else if (shape === 'star') {
                  ctx.beginPath();
                  for (var i = 0; i < 5; i++) {
                    var outerA = (i * 72 - 90) * Math.PI / 180;
                    var innerA = ((i * 72) + 36 - 90) * Math.PI / 180;
                    ctx.lineTo(cx + hs + hs * Math.cos(outerA), cy + hs + hs * Math.sin(outerA));
                    ctx.lineTo(cx + hs + hs * 0.4 * Math.cos(innerA), cy + hs + hs * 0.4 * Math.sin(innerA));
                  }
                  ctx.closePath();
                  ctx.fill();
                } else if (shape === 'hexagon') {
                  ctx.beginPath();
                  for (var i = 0; i < 6; i++) {
                    var a = (i * 60 - 30) * Math.PI / 180;
                    ctx.lineTo(cx + hs + hs * Math.cos(a), cy + hs + hs * Math.sin(a));
                  }
                  ctx.closePath();
                  ctx.fill();
                } else {
                  ctx.fillRect(cx, cy, size, size);
                }
              }

              // Section 1: Node Types
              ctx.fillStyle = '#0f172a';
              ctx.font = 'bold 13px sans-serif';
              ctx.fillText('Node Types', lx + pad, curY + 12);
              curY += lineH;

              ctx.font = '13px sans-serif';
              currentGroups.forEach(function(g) {
                drawShapeIcon(lx + pad, curY, g.shape, g.color, iconSize);
                ctx.fillStyle = '#0f172a';
                ctx.fillText(g.name, lx + textX, curY + 12);
                curY += lineH;
              });

              // Divider
              curY += sectionGap;

              // Section 2: Pathway Border
              ctx.fillStyle = '#0f172a';
              ctx.font = 'bold 13px sans-serif';
              ctx.fillText('Pathway Border', lx + pad, curY + 12);
              curY += lineH;

              ctx.font = '13px sans-serif';
              borderItems.forEach(function(bt) {
                ctx.strokeStyle = bt.color;
                ctx.lineWidth = 3;
                ctx.beginPath();
                ctx.arc(lx + pad + 7, curY + 7, 6, 0, 2 * Math.PI);
                ctx.stroke();
                ctx.fillStyle = '#0f172a';
                ctx.fillText(bt.label, lx + textX, curY + 12);
                curY += lineH;
              });

              // Divider
              curY += sectionGap;

              // Section 3: Edge Color
              ctx.fillStyle = '#0f172a';
              ctx.font = 'bold 13px sans-serif';
              ctx.fillText('Edge Color', lx + pad, curY + 12);
              curY += lineH;

              ctx.font = '13px sans-serif';
              edgeItems.forEach(function(et) {
                ctx.strokeStyle = et.color;
                ctx.fillStyle = et.color;
                ctx.lineWidth = 3;
                ctx.beginPath();
                ctx.moveTo(lx + pad, curY + 7);
                ctx.lineTo(lx + pad + 20, curY + 7);
                ctx.stroke();
                ctx.fillStyle = '#0f172a';
                ctx.fillText(et.label, lx + textX, curY + 12);
                curY += lineH;
              });

              var safeName = data ? data : 'network';
              var link = document.createElement('a');
              link.download = 'fig_mln_' + safeName + '.png';
              link.href = tempCanvas.toDataURL('image/png');
              link.click();
            };
            panel.appendChild(saveBtn);

            wrapper.appendChild(panel);
            el.appendChild(wrapper);

            // --- Legend Panel (bottom-right) ---
            var legend = document.createElement('div');
            legend.setAttribute('data-nuimm-legend', 'true');
            legend.style.position = 'absolute';
            legend.style.right = '20px';
            legend.style.bottom = '20px';
            legend.style.zIndex = '99998';
            legend.style.fontFamily = 'sans-serif';
            legend.style.backgroundColor = 'rgba(255, 255, 255, 0.95)';
            legend.style.padding = '14px 18px';
            legend.style.border = '1px solid #cbd5e1';
            legend.style.borderRadius = '8px';
            legend.style.boxShadow = '0 4px 6px rgba(0,0,0,0.1)';
            legend.style.fontSize = '13px';
            legend.style.color = '#0f172a';
            legend.style.lineHeight = '1.6';

            // Section 1: Node Types
            var s1Title = document.createElement('div');
            s1Title.innerHTML = '<b>Node Types</b>';
            s1Title.style.marginBottom = '6px';
            legend.appendChild(s1Title);

            groups.forEach(function(g) {
              var row = document.createElement('div');
              row.style.display = 'flex';
              row.style.alignItems = 'center';
              row.style.marginBottom = '4px';
              var icon = document.createElement('div');
              applyShapeCSS(icon, g.shape, g.color);
              legendIcons[g.name] = icon;
              var lbl = document.createElement('span');
              lbl.innerText = g.name;
              row.appendChild(icon);
              row.appendChild(lbl);
              legend.appendChild(row);
            });

            // Divider
            var div1 = document.createElement('div');
            div1.style.height = '16px';
            legend.appendChild(div1);

            // Section 2: Pathway Border (GSEA Direction)
            var s2Title = document.createElement('div');
            s2Title.innerHTML = '<b>Pathway Border</b>';
            s2Title.style.marginBottom = '6px';
            legend.appendChild(s2Title);

            var borderTypes = [
              {label: 'Upregulated', color: '#B5654A'},
              {label: 'Downregulated', color: '#5A7D8B'}
            ];
            borderTypes.forEach(function(bt) {
              var row = document.createElement('div');
              row.style.display = 'flex';
              row.style.alignItems = 'center';
              row.style.marginBottom = '4px';
              var icon = document.createElement('div');
              icon.style.width = '14px';
              icon.style.height = '14px';
              icon.style.borderRadius = '50%';
              icon.style.backgroundColor = 'transparent';
              icon.style.border = '3px solid ' + bt.color;
              icon.style.marginRight = '8px';
              icon.style.flexShrink = '0';
              icon.style.boxSizing = 'border-box';
              var lbl = document.createElement('span');
              lbl.innerText = bt.label;
              row.appendChild(icon);
              row.appendChild(lbl);
              legend.appendChild(row);
            });

            // Divider
            var div2 = document.createElement('div');
            div2.style.height = '16px';
            legend.appendChild(div2);

            // Section 3: Edge Color (PMN Correlation)
            var s3Title = document.createElement('div');
            s3Title.innerHTML = '<b>Edge Color</b>';
            s3Title.style.marginBottom = '6px';
            legend.appendChild(s3Title);

            var edgeTypes = [
              {label: 'Positive correlation', color: '#B5654A'},
              {label: 'Negative correlation', color: '#5A7D8B'}
            ];
            edgeTypes.forEach(function(et) {
              var row = document.createElement('div');
              row.style.display = 'flex';
              row.style.alignItems = 'center';
              row.style.marginBottom = '4px';
              var line = document.createElement('div');
              line.style.width = '20px';
              line.style.height = '3px';
              line.style.backgroundColor = et.color;
              line.style.marginRight = '8px';
              line.style.flexShrink = '0';
              line.style.borderRadius = '2px';
              var lbl = document.createElement('span');
              lbl.innerText = et.label;
              row.appendChild(line);
              row.appendChild(lbl);
              legend.appendChild(row);
            });

            el.appendChild(legend);
          ", get_ctrl_drag_js(), "
          }
          ")



          # Render network
          vis_plot <- visNetwork::visNetwork(nodes_df, edges, width = "100%", height = "95vh") |>
            visNetwork::visNodes(font = list(color = "#0f172a", size = 35, face = "sans-serif", background = "rgba(255,255,255,0.85)"), borderWidth = 1.5, shadow = TRUE) |>
            visNetwork::visEdges(smooth = FALSE, color = list(color = "rgba(160, 160, 160, 0.5)", highlight = "#e11d48", inherit = FALSE)) |>
            visNetwork::visGroups(groupname = "Microbe",   color = list(background = unname(node_colors['Microbe']),    border = darken_color(node_colors['Microbe']),    highlight = lighten_color(node_colors['Microbe'])),    shape = unname(node_shapes['Microbe'])) |>
            visNetwork::visGroups(groupname = "Pathway",    color = list(background = unname(node_colors['Pathway']),     border = darken_color(node_colors['Pathway']),     highlight = lighten_color(node_colors['Pathway'])),     shape = unname(node_shapes['Pathway'])) |>
            visNetwork::visGroups(groupname = "Metabolite", color = list(background = unname(node_colors['Metabolite']),  border = darken_color(node_colors['Metabolite']),  highlight = lighten_color(node_colors['Metabolite'])),  shape = unname(node_shapes['Metabolite'])) |>
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
    return(if (file.exists(out_html)) out_html else out_csv)
  } else {
    message("  No edges could be assembled for this comparison; skipping output.")
    return(NULL)
  }
}
