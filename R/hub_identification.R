utils::globalVariables(c("MCC_score", "mcc_score", "name"))

#' Hub Identification using Maximal Clique Centrality (MCC) algorithm
#'
#' @details
#' Identifies influential hub nodes in a multi-layered network using the
#' Maximal Clique Centrality (MCC) algorithm. Exports a ranked CSV of all
#' node scores and optionally generates an interactive HTML network colored
#' by MCC score with a customizable palette.
#'
#' @param multi_layered_network_file Path to the multi-layered network CSV/TSV
#'   (output from `con_mln`). Must contain 'from'/'to' or 'Feature1'/'Feature2' columns.
#' @param output_directory Path to save results (CSV + HTML).
#' @param top_n_hubs Integer. Number of top hub nodes to highlight in the HTML visualization.
#' @param visualize Logical. If TRUE, generates an interactive HTML network. Defaults to TRUE.
#' @return Invisible NULL. Results are written to `output_directory`.
#' @export
iden_hub <- function(
  multi_layered_network_file,
  output_directory,
  top_n_hubs = 5,
  visualize = TRUE
) {
  message("Initiating hub identification (Maximal Clique Centrality).")

  input_file_base_name <- tools::file_path_sans_ext(basename(multi_layered_network_file))
  cleaned_input_file_name <- gsub("[^A-Za-z0-9_]", "", input_file_base_name)

  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }

  message("[1/5] Loading network data.")
  if (!file.exists(multi_layered_network_file)) stop("File not found.")

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

  g <- igraph::graph_from_data_frame(d = network_data[, c(source_col, target_col)], directed = FALSE)
  message(sprintf("  Network graph instantiated: |V|=%d.", igraph::vcount(g)))

  message("[2/5] Computing maximal cliques.")
  cliques <- tryCatch(
    igraph::max_cliques(g),
    error = function(e) stop("Error finding cliques: ", e$message)
  )

  message("[3/5] Calculating MCC scores.")
  mcc_scores <- setNames(numeric(igraph::vcount(g)), igraph::V(g)$name)

  for (clique_nodes_indices in cliques) {
    clique_size <- length(clique_nodes_indices)
    clique_score <- factorial(clique_size - 1)
    for (node_index in clique_nodes_indices) {
      node_name <- igraph::V(g)$name[node_index]
      mcc_scores[node_name] <- mcc_scores[node_name] + clique_score
    }
  }

  node_degrees <- igraph::degree(g)
  node_clustering_coeffs <- igraph::transitivity(g, type = "local", vids = igraph::V(g))
  names(node_clustering_coeffs) <- igraph::V(g)$name

  for (node_name in igraph::V(g)$name) {
    current_degree <- node_degrees[node_name]
    is_no_edge <- FALSE
    if (current_degree <= 1) {
      is_no_edge <- TRUE
    } else {
      if (!is.na(node_clustering_coeffs[node_name]) && node_clustering_coeffs[node_name] == 0) {
        is_no_edge <- TRUE
      }
    }
    if (is_no_edge) {
      mcc_scores[node_name] <- current_degree
    }
  }

  message("[4/5] Ranking nodes and exporting results.")
  full_results_df <- dplyr::arrange(
    data.frame(Node = names(mcc_scores), MCC_score = mcc_scores, stringsAsFactors = FALSE),
    dplyr::desc(MCC_score)
  )

  output_filepath <- file.path(output_directory, paste0("hub_mcc_", cleaned_input_file_name, ".csv"))
  message("  Results saved: ", basename(output_filepath))
  write.csv(full_results_df, output_filepath, row.names = FALSE)

  if (visualize) {
    message("[5/5] Generating network visualization.")
    if (!requireNamespace("visNetwork", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) {
      stop("To use visualize=TRUE, you must install 'visNetwork' and 'htmlwidgets'.")
    }

    nodes_df <- data.frame(
      id = full_results_df$Node,
      label = sapply(full_results_df$Node, function(x) {
        if (grepl("d__|p__|c__|o__|f__|g__|s__|Bacteria", x)) clean_taxonomy(x) else x
      }),
      mcc_score = full_results_df$MCC_score,
      title = paste0("<div style='padding:10px; font-family:sans-serif;'><b>ID:</b> ", full_results_df$Node, "<br><b>MCC Score:</b> ", round(full_results_df$MCC_score, 2), "</div>"),
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

    # Initialize node coordinates for circular layout
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

    max_y <- max(nodes_df$y, na.rm = TRUE)
    legend_y <- max_y + 400

    legend_nodes <- data.frame(
      id = c("LEG_MIC", "LEG_PATH", "LEG_MET"),
      label = c("Microbe", "Pathway", "Metabolite"),
      mcc_score = c(0, 0, 0),
      title = c("", "", ""),
      group = c("Microbe", "Pathway", "Metabolite"),
      size = c(60, 60, 60),
      x = c(-300, 0, 300),
      y = c(legend_y, legend_y, legend_y),
      stringsAsFactors = FALSE
    )
    nodes_df <- rbind(nodes_df, legend_nodes)

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
      addHTML(panel, '<div style=\"font-size:14px;color:#0f172a;margin-bottom:10px\"><b>Hub Identification Analysis</b></div>');

      var visEngine = this.network;
      var nodesDS = visEngine.body.data.nodes;
      var edgesDS = visEngine.body.data.edges;
      var allNodes = nodesDS.get();
      var allEdges = edgesDS.get();

      var initialCoords = {};
      allNodes.forEach(function(n) { initialCoords[n.id] = {x:n.x, y:n.y}; });

      // Top Hubs input
      var topNWrap = document.createElement('div');
      topNWrap.style.cssText = 'margin-bottom:8px;font-size:13px;';
      topNWrap.innerHTML = '<b>Top Hubs:</b> ';
      var topNInput = document.createElement('input');
      topNInput.type = 'number'; topNInput.value = '", top_n_hubs, "';
      topNInput.style.cssText = 'width:50px;margin-right:10px;padding:2px;border:1px solid #cbd5e1;border-radius:4px;';
      topNWrap.appendChild(topNInput);
      panel.appendChild(topNWrap);

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

      // Palette selector
      var palWrap = document.createElement('div');
      palWrap.style.display = 'none';
      addHTML(palWrap, '<div style=\"font-size:13px;margin-bottom:4px\"><b>Hub Color Tone:</b></div>');
      var palSel = document.createElement('select');
      palSel.style.cssText = 'width:100%;padding:2px;border-radius:4px;border:1px solid #cbd5e1;margin-bottom:8px;';
      ['Plasma','Viridis','Blues','Reds'].forEach(function(p){ var o=document.createElement('option'); o.value=p; o.text=p; palSel.appendChild(o); });
      palWrap.appendChild(palSel);
      panel.appendChild(palWrap);

      // Legend
      var legendWrap = document.createElement('div');
      legendWrap.style.cssText = 'margin-bottom:8px;font-size:12px;display:none;';
      legendWrap.innerHTML = '<div style=\"text-align:center;margin-bottom:5px\"><b>MCC Score</b></div>';
      var gradBar = document.createElement('div');
      gradBar.style.cssText = 'height:15px;width:100%;border-radius:4px;border:1px solid #ccc;';
      legendWrap.appendChild(gradBar);
      var legLabels = document.createElement('div');
      legLabels.style.cssText = 'display:flex;justify-content:space-between;margin-top:2px;';
      var minLbl = document.createElement('span'); minLbl.innerHTML = 'Min';
      var maxLbl = document.createElement('span'); maxLbl.innerHTML = 'Max';
      legLabels.appendChild(minLbl); legLabels.appendChild(maxLbl);
      legendWrap.appendChild(legLabels);
      panel.appendChild(legendWrap);

      // Customize Network section
      addHTML(panel, '<hr style=\"margin:10px 0;border:0;border-top:1px solid #e2e8f0\">');
      addHTML(panel, '<div style=\"font-size:14px;color:#0f172a;margin-bottom:8px\"><b>Customize Network</b></div>');

      var groups = [{name:'Microbe',shape:'triangle'},{name:'Pathway',shape:'dot'},{name:'Metabolite',shape:'square'}];
      groups.forEach(function(g) {
        var w = document.createElement('div');
        w.style.cssText = 'display:flex;align-items:center;margin-bottom:6px;font-size:13px;';
        var l = document.createElement('div'); l.innerText = g.name+':'; l.style.cssText = 'width:80px;font-weight:bold;'; w.appendChild(l);
        var ss = document.createElement('select');
        ['dot','triangle','square','diamond','star','hexagon'].forEach(function(s){ var o=document.createElement('option'); o.value=s; o.text=s; if(s===g.shape) o.selected=true; ss.appendChild(o); });
        ss.style.cssText = 'padding:2px;border-radius:4px;border:1px solid #cbd5e1;width:120px;'; w.appendChild(ss);
        panel.appendChild(w);
        ss.addEventListener('change', function(){ var opts={groups:{}}; opts.groups[g.name]={shape:ss.value}; if(visEngine&&visEngine.setOptions) visEngine.setOptions(opts); });
      });

      // Save Network button
      var saveNetBtn = document.createElement('button');
      saveNetBtn.innerHTML = 'Save Network';
      saveNetBtn.style.cssText = 'margin-top:10px;padding:8px 16px;background:#f1f5f9;color:#0f172a;border:1px solid #cbd5e1;border-radius:6px;cursor:pointer;font-weight:bold;width:100%;margin-bottom:4px;';
      saveNetBtn.onclick = function() {
        var c = el.getElementsByTagName('canvas')[0]; if(!c) return;
        var tc = document.createElement('canvas'); tc.width=c.width; tc.height=c.height;
        var ctx=tc.getContext('2d'); ctx.fillStyle='#fff'; ctx.fillRect(0,0,tc.width,tc.height); ctx.drawImage(c,0,0);

        if (legendWrap.style.display !== 'none') {
            var pal = palSel.value;
            var scale = Math.max(1, tc.width / 1200);
            var w = 250 * scale, h = 15 * scale, px = tc.width - w - (40 * scale), py = tc.height - (60 * scale);
            ctx.fillStyle = '#0f172a'; ctx.font = 'bold ' + (16 * scale) + 'px sans-serif'; ctx.textAlign = 'center';
            ctx.fillText('MCC Score', px + w/2, py - (10 * scale));
            var grd = ctx.createLinearGradient(px, 0, px+w, 0);
            var cm = colorMap[pal];
            for(var i=0; i<cm.length; i++) grd.addColorStop(i/(cm.length-1), 'rgb('+cm[i].join(',')+')');
            ctx.fillStyle = grd; ctx.fillRect(px, py, w, h);
            var minM = parseFloat(minLbl.innerHTML), maxM = parseFloat(maxLbl.innerHTML);
            ctx.fillStyle = '#475569'; ctx.font = (12 * scale) + 'px sans-serif';
            ctx.textAlign = 'left'; ctx.fillText(minM.toFixed(1), px, py + h + (18 * scale));
            ctx.textAlign = 'right'; ctx.fillText(maxM.toFixed(1), px+w, py + h + (18 * scale));
        }

        var a=document.createElement('a'); a.download='hub_network.png'; a.href=tc.toDataURL('image/png'); a.click();
      };
      panel.appendChild(saveNetBtn);

      wrapper.appendChild(panel);
      el.appendChild(wrapper);

      // Color maps
      var colorMap = {
        'Plasma': [[13,8,135],[126,3,168],[204,71,120],[248,149,64],[240,249,33]],
        'Viridis': [[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]],
        'Blues': [[247,251,255],[198,219,239],[107,174,214],[33,113,181],[8,48,107]],
        'Reds': [[255,245,240],[254,224,210],[252,146,114],[222,45,38],[165,15,21]]
      };
      function getColor(val,mn,mx,pal){
        var pct=(mx===mn)?0.5:(val-mn)/(mx-mn); pct=Math.max(0,Math.min(1,pct));
        var s=colorMap[pal],ns=s.length,idx=pct*(ns-1),i=Math.floor(idx),f=idx-i;
        if(i>=ns-1) return 'rgb('+s[ns-1].join(',')+')';
        return 'rgb('+Math.round(s[i][0]+f*(s[i+1][0]-s[i][0]))+','+Math.round(s[i][1]+f*(s[i+1][1]-s[i][1]))+','+Math.round(s[i][2]+f*(s[i+1][2]-s[i][2]))+')'
      }
      function getGrad(pal){
        var s=colorMap[pal],p=[];
        for(var i=0;i<s.length;i++) p.push('rgb('+s[i].join(',')+') '+(i/(s.length-1)*100)+'%');
        return 'linear-gradient(to right,'+p.join(',')+')';
      }

      // Hub Analysis Engine
      function applyAnalysis() {
        var n = parseInt(topNInput.value);
        if(isNaN(n) || n <= 0) { alert('Please enter a valid number of top hubs.'); return; }

        // Remove legend nodes from canvas
        if (nodesDS.get('LEG_MIC')) nodesDS.remove(['LEG_MIC', 'LEG_PATH', 'LEG_MET']);

        var realNodes = allNodes.filter(function(nd){ return nd.id.indexOf('LEG_') !== 0; });
        realNodes.sort(function(a, b) { return b.mcc_score - a.mcc_score; });
        var topNodes = realNodes.slice(0, Math.min(n, realNodes.length));
        var minMcc = topNodes[topNodes.length-1] ? topNodes[topNodes.length-1].mcc_score : 0;
        var maxMcc = topNodes[0] ? topNodes[0].mcc_score : 1;

        var micN=[],pathN=[],metN=[];
        topNodes.forEach(function(s){
          if(s.group==='Microbe') micN.push(s);
          else if(s.group==='Pathway') pathN.push(s);
          else if(s.group==='Metabolite') metN.push(s);
        });

        var rm=200+(micN.length*15), rp=150+(pathN.length*20), rme=100+(metN.length*25);
        var xm=-(rm+rp+500), xp=0, xme=(rp+rme+500);
        function assignC(arr,cx,r){ for(var i=0;i<arr.length;i++){ var a=(i/arr.length)*2*Math.PI; arr[i].nx=cx+r*Math.cos(a); arr[i].ny=r*Math.sin(a); }}
        assignC(micN,xm,rm); assignC(pathN,xp,rp); assignC(metN,xme,rme);

        var topMap = {};
        topNodes.forEach(function(s){ topMap[s.id]=s; });

        var pal = palSel.value, updates = [];
        allNodes.forEach(function(nd) {
          if(nd.id.indexOf('LEG_')===0) return;
          if(topMap[nd.id]) {
            var s=topMap[nd.id], c=getColor(s.mcc_score,minMcc,maxMcc,pal);
            updates.push({id:nd.id, hidden:false, x:s.nx, y:s.ny, color:{background:c,border:'#475569',highlight:c}});
          } else {
            updates.push({id:nd.id, hidden:true});
          }
        });
        nodesDS.update(updates);

        palWrap.style.display = 'block';
        legendWrap.style.display = 'block';
        gradBar.style.background = getGrad(pal);
        minLbl.innerHTML = minMcc.toFixed(1);
        maxLbl.innerHTML = maxMcc.toFixed(1);

        if(visEngine && visEngine.fit) visEngine.fit({animation:{duration:1000,easingFunction:'easeInOutQuad'}});
      }

      runBtn.onclick = applyAnalysis;
      palSel.onchange = function() { applyAnalysis(); };

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
        var nUpdates = [];
        allNodes.forEach(function(n) {
          if (n.id.indexOf('LEG_') === 0) return;
          var ic = initialCoords[n.id];
          nUpdates.push({id:n.id, hidden:false, color:null, x:ic?ic.x:0, y:ic?ic.y:0});
        });
        nodesDS.update(nUpdates);
        palWrap.style.display = 'none';
        legendWrap.style.display = 'none';
        if(visEngine && visEngine.fit) visEngine.fit({animation:{duration:1000,easingFunction:'easeInOutQuad'}});
      };
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
    out_html <- file.path(output_directory, paste0("interactive_hub_", cleaned_input_file_name, ".html"))
    save_widget_safe(vis_plot, file = out_html, title = "NUIMM")
  }

  message("Hub identification completed.")
  invisible(NULL)
}
