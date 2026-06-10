utils::globalVariables(c("from", "to", "edge_score", "x", "y", "xend", "yend", "type", "name", "is_path"))

#' Pathfinding using Dijkstra's algorithm
#'
#' @details
#' This function calculates the shortest path between two nodes using Dijkstra's
#' algorithm, contextualized within a comprehensive overlay visualization.
#'
#' @param multi_layered_network_file Path to input CSV/TSV.
#' @param source_node Name of the starting node (or NULL for interactive selection).
#' @param target_node Name of the ending node (or NULL for interactive selection).
#' @param output_directory Path to save results.
#' @param visualize Logical. If TRUE, generates interactive HTML. Defaults to TRUE.
#' @return Invisible NULL. Results are written to `output_directory`.
#' @export
find_path <- function(
  multi_layered_network_file,
  source_node = NULL,
  target_node = NULL,
  output_directory,
  visualize = TRUE
) {
  message("Initiating shortest path analysis (Dijkstra algorithm).")

  input_file_base_name <- tools::file_path_sans_ext(basename(multi_layered_network_file))
  cleaned_input_file_name <- gsub("[^A-Za-z0-9_]", "", input_file_base_name)

  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)

  message("[1/2] Loading network data.")
  network_data <- read_input_file(multi_layered_network_file, stringsAsFactors = FALSE)

  if (all(c("from", "to") %in% colnames(network_data))) {
    source_col <- "from"
    target_col <- "to"
  } else if (all(c("Feature1", "Feature2") %in% colnames(network_data))) {
    source_col <- "Feature1"
    target_col <- "Feature2"
  } else {
    stop("Network file must contain columns: 'from' and 'to', OR 'Feature1' and 'Feature2'")
  }

  g <- igraph::graph_from_data_frame(d = network_data, directed = FALSE)

  if ("edge_score" %in% colnames(network_data)) {
    network_data$edge_score <- as.numeric(network_data$edge_score)
    network_data$edge_score_abs <- abs(network_data$edge_score)
    igraph::E(g)$weight <- sapply(network_data$edge_score_abs, function(w) {
      if (is.na(w)) Inf else if (w < 1) 1 / w else 1 / (w + 0.1)
    })
  }

  message(sprintf("  Network graph instantiated: |V|=%d, |E|=%d.", igraph::vcount(g), igraph::ecount(g)))

  # Build nodes data frame for visualization
  nodes_df <- data.frame(
    id = igraph::V(g)$name,
    label = sapply(igraph::V(g)$name, function(x) {
      if (grepl("d__|p__|c__|o__|f__|g__|s__|Bacteria", x)) clean_taxonomy(x) else x
    }),
    title = paste0("<div style='padding:10px; font-family:sans-serif;'><b>", igraph::V(g)$name, "</b></div>"),
    stringsAsFactors = FALSE
  )

  nodes_df$group <- as.character(determine_node_groups(nodes_df$id, network_data, source_col, target_col))

  nodes_df$size <- c("Microbe" = 20, "Pathway" = 30, "Metabolite" = 40)[nodes_df$group]

  edges_df <- data.frame(
    from = network_data[[source_col]],
    to = network_data[[target_col]],
    stringsAsFactors = FALSE
  )
  edge_score_col <- if ("edge_score" %in% colnames(network_data)) "edge_score" else if ("value" %in% colnames(network_data)) "value" else NULL
  if (!is.null(edge_score_col)) {
    edges_df$value <- abs(as.numeric(network_data[[edge_score_col]]))
    edges_df$title <- paste0("<div style='padding:10px; font-family:sans-serif;'><b>Value:</b> ", round(as.numeric(network_data[[edge_score_col]]), 4), "</div>")
  } else {
    edges_df$title <- "<div style='padding:10px; font-family:sans-serif;'><b>Value:</b> 1</div>"
  }

  # Circle layout
  nodes_df$x <- 0; nodes_df$y <- 0
  idx_mic <- which(nodes_df$group == "Microbe")
  idx_path <- which(nodes_df$group == "Pathway")
  idx_met <- which(nodes_df$group == "Metabolite")

  r_mic <- 200 + length(idx_mic) * 15
  r_path <- 150 + length(idx_path) * 20
  r_met <- 100 + length(idx_met) * 25
  x_mic <- -(r_mic + r_path + 500); x_path <- 0; x_met <- r_path + r_met + 500

  if (length(idx_mic) > 0) { ang <- seq(0, 2*pi, length.out=length(idx_mic)+1)[1:length(idx_mic)]; nodes_df$x[idx_mic] <- x_mic + r_mic*cos(ang); nodes_df$y[idx_mic] <- r_mic*sin(ang) }
  if (length(idx_path) > 0) { ang <- seq(0, 2*pi, length.out=length(idx_path)+1)[1:length(idx_path)]; nodes_df$x[idx_path] <- x_path + r_path*cos(ang); nodes_df$y[idx_path] <- r_path*sin(ang) }
  if (length(idx_met) > 0) { ang <- seq(0, 2*pi, length.out=length(idx_met)+1)[1:length(idx_met)]; nodes_df$x[idx_met] <- x_met + r_met*cos(ang); nodes_df$y[idx_met] <- r_met*sin(ang) }

  max_y <- max(nodes_df$y, na.rm = TRUE)
  legend_y <- max_y + 400
  legend_nodes <- data.frame(
    id = c("LEG_MIC", "LEG_PATH", "LEG_MET"),
    label = c("Microbe", "Pathway", "Metabolite"),
    title = c("", "", ""),
    group = c("Microbe", "Pathway", "Metabolite"),
    size = c(60, 60, 60),
    x = c(-300, 0, 300),
    y = c(legend_y, legend_y, legend_y),
    stringsAsFactors = FALSE
  )
  nodes_df <- rbind(nodes_df, legend_nodes)

  if (visualize) {
    if (!requireNamespace("visNetwork", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) {
      stop("To use visualize=TRUE, you must install 'visNetwork' and 'htmlwidgets'.")
    }

    js_custom_panel <- paste0("
    function(el, x, data) {
      document.body.style.overflow = 'hidden';
      var wrapper = document.createElement('div');
      wrapper.style.cssText = 'position:absolute;left:20px;bottom:20px;z-index:99999;font-family:sans-serif;';

      var toggleBtn = document.createElement('button');
      toggleBtn.innerHTML = 'Hide Controls';
      toggleBtn.style.cssText = 'padding:8px 16px;background:#f1f5f9;color:#0f172a;border:1px solid #cbd5e1;border-radius:6px;cursor:pointer;font-weight:bold;display:block;margin-bottom:10px;box-shadow:0 4px 6px rgba(0,0,0,0.1);';
      wrapper.appendChild(toggleBtn);

      var panel = document.createElement('div');
      panel.style.cssText = 'background:rgba(255,255,255,0.95);padding:15px;border:1px solid #cbd5e1;border-radius:8px;box-shadow:0 4px 6px rgba(0,0,0,0.1);max-height:85vh;overflow-y:auto;';

      var isPanelOpen = true;
      toggleBtn.onclick = function() {
        isPanelOpen = !isPanelOpen;
        panel.style.display = isPanelOpen ? 'block' : 'none';
        toggleBtn.innerHTML = isPanelOpen ? 'Hide Controls' : 'Customize Network';
      };

      function addHTML(parent, htmlString) {
          var temp = document.createElement('div');
          temp.innerHTML = htmlString;
          while (temp.firstChild) {
              parent.appendChild(temp.firstChild);
          }
      }

      addHTML(panel, '<div style=\"font-size:14px;color:#475569;margin-bottom:10px\"><b>Tip:</b> Scroll to zoom. Hold Ctrl + drag to select multiple nodes. Drag nodes to perfect layout.<hr style=\"margin:10px 0;border:0;border-top:1px solid #e2e8f0\"></div>');
      addHTML(panel, '<div style=\"font-size:14px;color:#0f172a;margin-bottom:10px\"><b>Pathfinding Analysis</b></div>');

      var visEngine = this.network;
      var nodesDS = visEngine.body.data.nodes;
      var edgesDS = visEngine.body.data.edges;
      var allNodes = nodesDS.get();
      var allEdges = edgesDS.get();

      var initialCoords = {};
      allNodes.forEach(function(n) { initialCoords[n.id] = {x:n.x, y:n.y}; });

      var realNodes = allNodes.filter(function(n){ return n.id.indexOf('LEG_') !== 0; });
      realNodes.sort(function(a,b){ return a.label.localeCompare(b.label); });

      // Source Node dropdown
      addHTML(panel, '<div style=\"font-size:13px;margin-bottom:4px\"><b>Source Node (Microbe):</b></div>');
      var srcSel = document.createElement('select');
      srcSel.style.cssText = 'width:100%;padding:4px;margin-bottom:8px;border:1px solid #cbd5e1;border-radius:4px;font-size:12px;';
      var optSrc0 = document.createElement('option'); optSrc0.value = ''; optSrc0.text = '-- Select Source --'; srcSel.appendChild(optSrc0);
      realNodes.forEach(function(n) {
        if (n.group === 'Microbe') {
          var opt = document.createElement('option'); opt.value = n.id; opt.text = n.label;
          srcSel.appendChild(opt);
        }
      });
      panel.appendChild(srcSel);

      // Target Node dropdown
      addHTML(panel, '<div style=\"font-size:13px;margin-bottom:4px\"><b>Target Node (Metabolite):</b></div>');
      var tgtSel = document.createElement('select');
      tgtSel.style.cssText = 'width:100%;padding:4px;margin-bottom:10px;border:1px solid #cbd5e1;border-radius:4px;font-size:12px;';
      var optTgt0 = document.createElement('option'); optTgt0.value = ''; optTgt0.text = '-- Select Target --'; tgtSel.appendChild(optTgt0);
      realNodes.forEach(function(n) {
        if (n.group === 'Metabolite') {
          var opt = document.createElement('option'); opt.value = n.id; opt.text = n.label;
          tgtSel.appendChild(opt);
        }
      });
      panel.appendChild(tgtSel);

      // Set default values if provided
      var defaultSrc = '", if (!is.null(source_node)) source_node else "", "';
      var defaultTgt = '", if (!is.null(target_node)) target_node else "", "';
      if (defaultSrc) srcSel.value = defaultSrc;
      if (defaultTgt) tgtSel.value = defaultTgt;

      // Run and Reset buttons
      var btnWrap = document.createElement('div');
      btnWrap.style.cssText = 'display:flex;gap:8px;margin-bottom:10px;';
      var runBtn = document.createElement('button');
      runBtn.innerHTML = 'Run Analysis';
      runBtn.style.cssText = 'padding:6px 12px;cursor:pointer;font-weight:bold;';
      var resetBtn = document.createElement('button');
      resetBtn.innerHTML = 'Reset';
      resetBtn.style.cssText = 'padding:6px 12px;cursor:pointer;font-weight:bold;';
      btnWrap.appendChild(runBtn);
      btnWrap.appendChild(resetBtn);
      panel.appendChild(btnWrap);

      // Customize Network section
      addHTML(panel, '<hr style=\"margin:10px 0;border:0;border-top:1px solid #e2e8f0\">');
      addHTML(panel, '<div style=\"font-size:14px;color:#0f172a;margin-bottom:8px\"><b>Customize Network</b></div>');

      var groups = [{name:'Microbe',color:'#9aa374',shape:'triangle'},{name:'Pathway',color:'#c1abad',shape:'dot'},{name:'Metabolite',color:'#4e7286',shape:'square'}];
      groups.forEach(function(g) {
        var w = document.createElement('div');
        w.style.cssText = 'display:flex;align-items:center;margin-bottom:6px;font-size:13px;';
        var l = document.createElement('div'); l.innerText = g.name+':'; l.style.cssText = 'width:80px;font-weight:bold;'; w.appendChild(l);

        var colorIn = document.createElement('input');
        colorIn.type = 'color'; colorIn.value = g.color;
        colorIn.style.cssText = 'width:25px;height:25px;padding:0;border:none;cursor:pointer;margin-right:8px;';
        w.appendChild(colorIn);

        var ss = document.createElement('select');
        ['dot','triangle','square','diamond','star','hexagon'].forEach(function(s){ var o=document.createElement('option'); o.value=s; o.text=s; if(s===g.shape) o.selected=true; ss.appendChild(o); });
        ss.style.cssText = 'padding:2px;border-radius:4px;border:1px solid #cbd5e1;width:90px;'; w.appendChild(ss);
        panel.appendChild(w);

        var updateGraph = function() {
          var opts = {groups:{}}; opts.groups[g.name] = {color:{background:colorIn.value,border:'#475569',highlight:colorIn.value}, shape:ss.value};
          if(visEngine && visEngine.setOptions) visEngine.setOptions(opts);
        };
        colorIn.addEventListener('change', updateGraph);
        ss.addEventListener('change', updateGraph);
      });

      // Save Network button
      var saveNetBtn = document.createElement('button');
      saveNetBtn.innerHTML = 'Save Network';
      saveNetBtn.style.cssText = 'margin-top:10px;padding:8px 16px;background:#f1f5f9;color:#0f172a;border:1px solid #cbd5e1;border-radius:6px;cursor:pointer;font-weight:bold;width:100%;margin-bottom:4px;';
      saveNetBtn.onclick = function() {
        var c = el.getElementsByTagName('canvas')[0]; if(!c) return;
        var tc = document.createElement('canvas'); tc.width=c.width; tc.height=c.height;
        var ctx=tc.getContext('2d'); ctx.fillStyle='#fff'; ctx.fillRect(0,0,tc.width,tc.height); ctx.drawImage(c,0,0);
        var a=document.createElement('a'); a.download='path_network.png'; a.href=tc.toDataURL('image/png'); a.click();
      };
      panel.appendChild(saveNetBtn);

      wrapper.appendChild(panel);
      el.appendChild(wrapper);

      // Dijkstra Engine
      function runDijkstra(src, tgt) {
        var dist = {}, prev = {}, prevEdge = {}, visited = {};
        realNodes.forEach(function(n) { dist[n.id] = Infinity; prev[n.id] = null; prevEdge[n.id] = null; });
        dist[src] = 0;

        var adj = {};
        realNodes.forEach(function(n) { adj[n.id] = []; });
        allEdges.forEach(function(e) {
          var w = 1;
          if (e.value) { var v = Math.abs(e.value); w = v < 1 ? 1/v : 1/(v+0.1); }
          adj[e.from].push({to: e.to, edgeId: e.id, weight: w});
          adj[e.to].push({to: e.from, edgeId: e.id, weight: w});
        });

        for (var i = 0; i < realNodes.length; i++) {
          var u = null, minD = Infinity;
          for (var id in dist) { if (!visited[id] && dist[id] < minD) { minD = dist[id]; u = id; } }
          if (u === null || u === tgt) break;
          visited[u] = true;
          if (adj[u]) adj[u].forEach(function(nb) {
            var alt = dist[u] + nb.weight;
            if (alt < dist[nb.to]) { dist[nb.to] = alt; prev[nb.to] = u; prevEdge[nb.to] = nb.edgeId; }
          });
        }

        var pathNodes = [], pathEdgeIds = [], pathDirs = {};
        var cur = tgt;
        while (cur !== null) {
          pathNodes.push(cur);
          if (prevEdge[cur]) {
            pathEdgeIds.push(prevEdge[cur]);
            pathDirs[prevEdge[cur]] = {from: prev[cur], to: cur};
          }
          cur = prev[cur];
        }

        if (pathNodes.indexOf(src) === -1) return null;
        return {nodes: pathNodes, edges: pathEdgeIds, dirs: pathDirs};
      }

      // Event Handlers
      resetBtn.onclick = function() {
        // Re-add legend nodes if removed
        if (!nodesDS.get('LEG_MIC')) {
          var ly = Math.max.apply(null, allNodes.map(function(n){return n.y||0;})) + 400;
          nodesDS.add([
            {id:'LEG_MIC',label:'Microbe',title:'',group:'Microbe',size:60,x:-300,y:ly},
            {id:'LEG_PATH',label:'Pathway',title:'',group:'Pathway',size:60,x:0,y:ly},
            {id:'LEG_MET',label:'Metabolite',title:'',group:'Metabolite',size:60,x:300,y:ly}
          ]);
        }
        var nUpdates = [], eUpdates = [];
        allNodes.forEach(function(n) {
          if (n.id.indexOf('LEG_') === 0) return;
          var ic = initialCoords[n.id];
          nUpdates.push({id:n.id, hidden:false, color:null, x:ic?ic.x:0, y:ic?ic.y:0});
        });
        allEdges.forEach(function(e) {
          eUpdates.push({id:e.id, hidden:false, color:null, width:null, arrows:''});
        });
        nodesDS.update(nUpdates);
        edgesDS.update(eUpdates);
        if(visEngine && visEngine.fit) visEngine.fit({animation:{duration:1000,easingFunction:'easeInOutQuad'}});
      };

      runBtn.onclick = function() {
        var src = srcSel.value, tgt = tgtSel.value;
        if (!src || !tgt) { alert('Please select both source and target nodes.'); return; }
        if (src === tgt) { alert('Source and target must be different.'); return; }

        // Remove legend nodes
        if (nodesDS.get('LEG_MIC')) nodesDS.remove(['LEG_MIC', 'LEG_PATH', 'LEG_MET']);

        runBtn.innerHTML = 'Running...'; runBtn.disabled = true;
        setTimeout(function() {
          var result = runDijkstra(src, tgt);

          if (!result) {
            alert('No path found between the selected nodes.');
            runBtn.innerHTML = 'Run Analysis'; runBtn.disabled = false;
            return;
          }

          var pathNodeSet = {};
          result.nodes.forEach(function(id) { pathNodeSet[id] = true; });
          var pathEdgeSet = {};
          result.edges.forEach(function(id) { pathEdgeSet[id] = true; });

          // Arrange path nodes horizontally
          var orderedPath = result.nodes.slice().reverse();
          var spacing = 300;
          var totalW = (orderedPath.length - 1) * spacing;
          var startX = -totalW / 2;

          var pathCoords = {};
          orderedPath.forEach(function(id, i) {
            pathCoords[id] = {x: startX + i * spacing, y: 0};
          });

          // Show path nodes and hide others
          var nUpdates = [];
          allNodes.forEach(function(nd) {
            if (nd.id.indexOf('LEG_') === 0) return;
            if (pathNodeSet[nd.id]) {
              var pos = pathCoords[nd.id];
              nUpdates.push({id:nd.id, hidden:false, x:pos.x, y:pos.y, size:30});
            } else {
              nUpdates.push({id:nd.id, hidden:true});
            }
          });

          // Highlight path edges
          var eUpdates = [];
          allEdges.forEach(function(e) {
            if (pathEdgeSet[e.id]) {
              var dir = result.dirs[e.id];
              eUpdates.push({id:e.id, hidden:false, from:dir.from, to:dir.to, color:{color:'#94a3b8',highlight:'#64748b'}, width:3, value:null, arrows:'to'});
            } else {
              eUpdates.push({id:e.id, hidden:true});
            }
          });
          nodesDS.update(nUpdates);
          edgesDS.update(eUpdates);

          if(visEngine && visEngine.fit) visEngine.fit({animation:{duration:1000,easingFunction:'easeInOutQuad'}});
          runBtn.innerHTML = 'Run Analysis'; runBtn.disabled = false;
        }, 50);
      };

      // Auto-run if default nodes provided
      if (defaultSrc && defaultTgt) {
        setTimeout(function() { runBtn.click(); }, 500);
      }
    ", get_ctrl_drag_js(), "
    }
    ")

    vis_plot <- visNetwork::visNetwork(nodes_df, edges_df, width = "100%", height = "95vh") |>
      visNetwork::visNodes(font = list(color = "#0f172a", size = 35, face = "sans-serif", background = "rgba(255,255,255,0.85)"), borderWidth = 1.5, shadow = TRUE) |>
      visNetwork::visEdges(smooth = FALSE, color = list(color = "rgba(180, 180, 180, 0.4)", highlight = "#e11d48"), width = 1) |>
      visNetwork::visGroups(groupname = "Microbe", color = list(background = "#9AA374", border = "#7A825C", highlight = "#B4BE89"), shape = "triangle") |>
      visNetwork::visGroups(groupname = "Pathway", color = list(background = "#C1ABAD", border = "#9A898A", highlight = "#D8C5C7"), shape = "dot") |>
      visNetwork::visGroups(groupname = "Metabolite", color = list(background = "#4E7286", border = "#3A5565", highlight = "#6392AB"), shape = "square") |>
      visNetwork::visInteraction(navigationButtons = FALSE, dragNodes = TRUE, multiselect = TRUE, hover = TRUE) |>
      visNetwork::visPhysics(enabled = FALSE) |>
      htmlwidgets::onRender(js_custom_panel)

    vis_plot$x$background <- "#ffffff"
    out_html <- file.path(output_directory, paste0("interactive_path_", cleaned_input_file_name, ".html"))
    save_widget_safe(vis_plot, file = out_html, title = "NUIMM")
    message("  Visualization saved: ", basename(out_html))
  }

  message("[2/2] Pathfinding algorithm completed.")
  invisible(NULL)
}
