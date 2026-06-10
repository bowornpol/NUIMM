utils::globalVariables(c(
  "edge_type", "Feature2", "Feature1", "from", "to",
  "Node_A", "Node_B", "edge_score", "Heat_score", "Time", "Correlation", "type"
))

#' Node Prioritization using Laplacian Heat Diffusion (LHD)
#'
#' @details
#' Implements LHD on a multi-layered network to prioritize nodes.
#' Generates an interactive HTML with client-side heat diffusion,
#' seed selection, edge filtering, and subnetwork visualization.
#'
#' @param multi_layered_network_file Path to the multi-layered network CSV/TSV.
#' @param output_directory Path to save results.
#' @param time_step_interval Numeric. Interval between time steps (default 0.01).
#' @param stabilization_threshold Numeric. Threshold for stabilization detection.
#' @param stabilization_window_size Integer. Window size for stabilization check.
#' @param visualize Logical. If TRUE, generates interactive HTML.
#' @return Invisible NULL.
#' @export
node_prior <- function(
  multi_layered_network_file,
  output_directory,
  time_step_interval = 0.01,
  stabilization_threshold = 0.0001,
  stabilization_window_size = 10,
  visualize = TRUE
) {
  message("Initiating node prioritization (Laplacian Heat Diffusion).")

  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)

  input_file_base_name <- tools::file_path_sans_ext(basename(multi_layered_network_file))
  cleaned_input_file_name <- gsub("[^A-Za-z0-9_]", "", input_file_base_name)

  message("[1/1] Loading network data.")
  if (!file.exists(multi_layered_network_file)) stop("File not found.")

  network_data <- read_input_file(multi_layered_network_file, stringsAsFactors = FALSE)

  if (all(c("from", "to") %in% colnames(network_data))) {
    source_col <- "from"
    target_col <- "to"
  } else if (all(c("Feature1", "Feature2") %in% colnames(network_data))) {
    source_col <- "Feature1"
    target_col <- "Feature2"
  } else {
    stop("Network file must have 'from'/'to' or 'Feature1'/'Feature2' columns.")
  }

  edge_score_col <- if ("edge_score" %in% colnames(network_data)) "edge_score" else if ("value" %in% colnames(network_data)) "value" else NULL

  edges_df <- data.frame(
    from = network_data[[source_col]],
    to = network_data[[target_col]],
    stringsAsFactors = FALSE
  )
  if (!is.null(edge_score_col)) {
    edges_df$value <- abs(as.numeric(network_data[[edge_score_col]]))
    edges_df$title <- paste0("<div style='padding:10px; font-family:sans-serif;'><b>Value:</b> ", round(as.numeric(network_data[[edge_score_col]]), 4), "</div>")
  } else {
    edges_df$title <- "<div style='padding:10px; font-family:sans-serif;'><b>Value:</b> 1</div>"
  }
  edges_df$edge_weight <- if (!is.null(edge_score_col)) abs(as.numeric(network_data[[edge_score_col]])) else rep(1, nrow(edges_df))

  unique_nodes <- unique(c(edges_df$from, edges_df$to))

  nodes_df <- data.frame(
    id = unique_nodes,
    label = sapply(unique_nodes, function(x) {
      if (grepl("d__|p__|c__|o__|f__|g__|s__|Bacteria", x)) clean_taxonomy(x) else x
    }),
    title = paste0("<div style='padding:10px; font-family:sans-serif;'><b>ID:</b> ", unique_nodes, "</div>"),
    stringsAsFactors = FALSE
  )

  nodes_df$group <- as.character(determine_node_groups(nodes_df$id, network_data, source_col, target_col))
  nodes_df$size <- c("Microbe" = 20, "Pathway" = 30, "Metabolite" = 40)[nodes_df$group]

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
    title = c("", "", ""),
    group = c("Microbe", "Pathway", "Metabolite"),
    size = c(60, 60, 60),
    x = c(-300, 0, 300), y = c(legend_y, legend_y, legend_y),
    stringsAsFactors = FALSE
  )
  nodes_df <- rbind(nodes_df, legend_nodes)

  if (visualize) {
    if (!requireNamespace("visNetwork", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) {
      stop("Install 'visNetwork' and 'htmlwidgets' for visualize=TRUE.")
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
      addHTML(panel, '<div style=\"font-size:14px;color:#0f172a;margin-bottom:10px\"><b>Node Prioritization Analysis</b></div>');

      var visEngine = this.network;
      var nodesDS = visEngine.body.data.nodes;
      var edgesDS = visEngine.body.data.edges;
      var allNodes = nodesDS.get();
      var allEdges = edgesDS.get();

      var initialCoords = {};
      allNodes.forEach(function(n) { initialCoords[n.id] = {x:n.x, y:n.y}; });

      var metNodes = allNodes.filter(function(n){ return n.group === 'Metabolite' && n.id.indexOf('LEG_') !== 0; });
      metNodes.sort(function(a,b){ return a.label.localeCompare(b.label); });

      // Seed Selection
      addHTML(panel, '<div style=\"font-size:13px;margin-bottom:4px\"><b>Select Seed Node (Metabolite):</b></div>');
      var seedContainer = document.createElement('div');
      seedContainer.style.cssText = 'width:100%;height:100px;overflow-y:auto;border:1px solid #cbd5e1;padding:4px;margin-bottom:8px;font-size:12px;background:#fff;';
      var seedCheckboxes = [];
      metNodes.forEach(function(n) {
        var lbl = document.createElement('label');
        lbl.style.cssText = 'display:block; cursor:pointer; margin-bottom:2px;';
        var cb = document.createElement('input');
        cb.type = 'checkbox';
        cb.value = n.id;
        cb.style.marginRight = '5px';
        seedCheckboxes.push(cb);
        lbl.appendChild(cb);
        lbl.appendChild(document.createTextNode(n.label));
        seedContainer.appendChild(lbl);
      });
      panel.appendChild(seedContainer);

      // Filter Toggle
      var filterWrap = document.createElement('div');
      filterWrap.style.cssText = 'font-size:13px;margin-bottom:8px;';
      filterWrap.innerHTML = '<b>Filter Other Metabolite Edges:</b> ';
      var filterSel = document.createElement('select');
      filterSel.style.cssText = 'padding:2px;border-radius:4px;border:1px solid #cbd5e1;';
      ['YES','NO'].forEach(function(v){ var o=document.createElement('option'); o.value=v; o.text=v; filterSel.appendChild(o); });
      filterWrap.appendChild(filterSel);
      panel.appendChild(filterWrap);

      // Run and Reset Buttons
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

      // Post-analysis controls
      var postWrap = document.createElement('div');
      postWrap.style.display = 'none';

      addHTML(postWrap, '<hr style=\"margin:10px 0;border:0;border-top:1px solid #e2e8f0\">');
      addHTML(postWrap, '<div style=\"font-size:13px;margin-bottom:4px\"><b>Show Top % Nodes:</b></div>');
      var pctSlider = document.createElement('input');
      pctSlider.type = 'range'; pctSlider.min = 1; pctSlider.max = 100; pctSlider.value = 100;
      pctSlider.style.cssText = 'width:100%;margin-bottom:2px;';
      postWrap.appendChild(pctSlider);
      var pctLabel = document.createElement('div');
      pctLabel.style.cssText = 'font-size:12px;text-align:center;margin-bottom:8px;';
      pctLabel.innerHTML = '100%';
      postWrap.appendChild(pctLabel);
      // Top Microbes input
      var topNWrap = document.createElement('div');
      topNWrap.style.cssText = 'margin-bottom:8px;font-size:13px;';
      topNWrap.innerHTML = '<b>Top Microbes:</b> ';
      var microbeInput = document.createElement('input');
      microbeInput.type = 'number'; microbeInput.value = '5';
      microbeInput.min = 1;
      microbeInput.style.cssText = 'width:50px;margin-right:10px;padding:2px;border:1px solid #cbd5e1;border-radius:4px;';
      topNWrap.appendChild(microbeInput);
      postWrap.appendChild(topNWrap);

      addHTML(postWrap, '<div style=\"font-size:13px;margin-bottom:4px\"><b>Heat Color Tone:</b></div>');
      var palSel = document.createElement('select');
      palSel.style.cssText = 'width:100%;padding:2px;border-radius:4px;border:1px solid #cbd5e1;margin-bottom:8px;';
      ['Plasma','Viridis','Blues','Reds'].forEach(function(p){ var o=document.createElement('option'); o.value=p; o.text=p; palSel.appendChild(o); });
      postWrap.appendChild(palSel);

      var legendWrap = document.createElement('div');
      legendWrap.style.cssText = 'margin-bottom:8px;font-size:12px;';
      legendWrap.innerHTML = '<div style=\"text-align:center;margin-bottom:5px\"><b>Heat Score</b></div>';
      var gradBar = document.createElement('div');
      gradBar.style.cssText = 'height:15px;width:100%;border-radius:4px;border:1px solid #ccc;';
      legendWrap.appendChild(gradBar);
      var legLabels = document.createElement('div');
      legLabels.style.cssText = 'display:flex;justify-content:space-between;margin-top:2px;';
      var minLbl = document.createElement('span'); minLbl.innerHTML = 'Min';
      var maxLbl = document.createElement('span'); maxLbl.innerHTML = 'Max';
      legLabels.appendChild(minLbl); legLabels.appendChild(maxLbl);
      legendWrap.appendChild(legLabels);
      postWrap.appendChild(legendWrap);

      // Correlation plot canvas
      addHTML(postWrap, '<div style=\"font-size:13px;margin-bottom:4px\"><b>Stabilization Curve</b></div>');
      var corrCanvas = document.createElement('canvas');
      corrCanvas.width = 1120; corrCanvas.height = 560;
      corrCanvas.style.cssText = 'width:100%;height:140px;border:1px solid #e2e8f0;border-radius:4px;margin-bottom:8px;background:#fff;';
      postWrap.appendChild(corrCanvas);
      var stabInfo = document.createElement('div');
      stabInfo.style.cssText = 'font-size:11px;color:#475569;margin-bottom:8px;';
      postWrap.appendChild(stabInfo);

      panel.appendChild(postWrap);

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
        var upd = function(){ var opts={groups:{}}; opts.groups[g.name]={shape:ss.value}; if(visEngine&&visEngine.setOptions) visEngine.setOptions(opts); };
        ss.addEventListener('change',upd);
      });

      // Save buttons
      var saveNetBtn = document.createElement('button');
      saveNetBtn.innerHTML = 'Save Network';
      saveNetBtn.style.cssText = 'margin-top:10px;padding:8px 16px;background:#f1f5f9;color:#0f172a;border:1px solid #cbd5e1;border-radius:6px;cursor:pointer;font-weight:bold;width:100%;margin-bottom:4px;';
      saveNetBtn.onclick = function() {
        var c = el.getElementsByTagName('canvas')[0]; if(!c) return;
        var tc = document.createElement('canvas'); tc.width=c.width; tc.height=c.height;
        var ctx=tc.getContext('2d'); ctx.fillStyle='#fff'; ctx.fillRect(0,0,tc.width,tc.height); ctx.drawImage(c,0,0);

        if (lastHeatResult && postWrap.style.display !== 'none') {
            var pal = palSel.value;
            var scale = Math.max(1, tc.width / 1200);
            var w = 250 * scale, h = 15 * scale, x = tc.width - w - (40 * scale), y = tc.height - (60 * scale);
            ctx.fillStyle = '#0f172a'; ctx.font = 'bold ' + (16 * scale) + 'px sans-serif'; ctx.textAlign = 'center';
            ctx.fillText('Heat Score', x + w/2, y - (10 * scale));
            var grd = ctx.createLinearGradient(x, 0, x+w, 0);
            var cm = colorMap[pal];
            for(var i=0; i<cm.length; i++) grd.addColorStop(i/(cm.length-1), 'rgb('+cm[i].join(',')+')');
            ctx.fillStyle = grd; ctx.fillRect(x, y, w, h);
            var minH = parseFloat(minLbl.innerHTML), maxH = parseFloat(maxLbl.innerHTML);
            ctx.fillStyle = '#475569'; ctx.font = (12 * scale) + 'px sans-serif';
            ctx.textAlign = 'left'; ctx.fillText(minH.toFixed(4), x, y + h + (18 * scale));
            ctx.textAlign = 'right'; ctx.fillText(maxH.toFixed(4), x+w, y + h + (18 * scale));
        }

        var a=document.createElement('a'); a.download='node_prior_network.png'; a.href=tc.toDataURL('image/png'); a.click();
      };
      panel.appendChild(saveNetBtn);

      var saveCsvBtn = document.createElement('button');
      saveCsvBtn.innerHTML = 'Save Heat Scores (CSV)';
      saveCsvBtn.style.cssText = 'padding:8px 16px;background:#f1f5f9;color:#0f172a;border:1px solid #cbd5e1;border-radius:6px;cursor:pointer;font-weight:bold;width:100%;margin-bottom:4px;';
      saveCsvBtn.style.display = 'none';
      panel.appendChild(saveCsvBtn);

      var savePlotBtn = document.createElement('button');
      savePlotBtn.innerHTML = 'Save Correlation Plot';
      savePlotBtn.style.cssText = 'padding:8px 16px;background:#f1f5f9;color:#0f172a;border:1px solid #cbd5e1;border-radius:6px;cursor:pointer;font-weight:bold;width:100%;';
      savePlotBtn.style.display = 'none';
      panel.appendChild(savePlotBtn);

      wrapper.appendChild(panel);
      el.appendChild(wrapper);

      // Color Maps
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
        return 'rgb('+Math.round(s[i][0]+f*(s[i+1][0]-s[i][0]))+','+Math.round(s[i][1]+f*(s[i+1][1]-s[i][1]))+','+Math.round(s[i][2]+f*(s[i+1][2]-s[i][2]))+')';
      }
      function getGrad(pal){
        var s=colorMap[pal],p=[];
        for(var i=0;i<s.length;i++) p.push('rgb('+s[i].join(',')+') '+(i/(s.length-1)*100)+'%');
        return 'linear-gradient(to right,'+p.join(',')+')';
      }

      // LHD Engine
      var DT = ", time_step_interval, ";
      var STAB_THRESH = ", stabilization_threshold, ";
      var STAB_WIN = ", stabilization_window_size, ";

      function rank(a){
        var s=a.map(function(v,i){return{v:v,i:i};}).sort(function(a,b){return a.v-b.v;});
        var r=new Array(a.length),i=0;
        while(i<s.length){ var j=i; while(j<s.length&&s[j].v===s[i].v) j++;
          var avg=(i+j+1)/2; for(var k=i;k<j;k++) r[s[k].i]=avg; i=j; }
        return r;
      }
      function pearson(a,b){
        var n=a.length,sa=0,sb=0,sab=0,sa2=0,sb2=0;
        for(var i=0;i<n;i++){sa+=a[i];sb+=b[i];sab+=a[i]*b[i];sa2+=a[i]*a[i];sb2+=b[i]*b[i];}
        var d=Math.sqrt((n*sa2-sa*sa)*(n*sb2-sb*sb));
        return d===0?1:(n*sab-sa*sb)/d;
      }
      function spearman(a,b){ return pearson(rank(a),rank(b)); }

      var lastHeatResult = null;

      function runLHD(seeds, filterMet) {
        var realNodes = allNodes.filter(function(n){ return n.id.indexOf('LEG_')!==0; });
        var nodeIds = realNodes.map(function(n){return n.id;});
        var nodeIdx = {}; nodeIds.forEach(function(id,i){nodeIdx[id]=i;});
        var n = nodeIds.length;

        var metIds = {};
        metNodes.forEach(function(m){ metIds[m.id]=true; });

        var adj = []; for(var i=0;i<n;i++){ adj[i]=[]; for(var j=0;j<n;j++) adj[i][j]=0; }

        allEdges.forEach(function(e) {
          if(filterMet) {
            var fromIsMet = metIds[e.from] && !seeds.some(function(s){return s===e.from;});
            var toIsMet = metIds[e.to] && !seeds.some(function(s){return s===e.to;});
            if(fromIsMet || toIsMet) return;
          }
          var fi = nodeIdx[e.from], ti = nodeIdx[e.to];
          if(fi!==undefined && ti!==undefined) {
            var w = e.edge_weight || e.value || 1;
            adj[fi][ti] = w; adj[ti][fi] = w;
          }
        });

        var deg = new Array(n).fill(0);
        for(var i=0;i<n;i++) for(var j=0;j<n;j++) deg[i]+=adj[i][j];

        var H = new Array(n).fill(0);
        var hps = 1.0/seeds.length;
        seeds.forEach(function(s){ var idx=nodeIdx[s]; if(idx!==undefined) H[idx]=hps; });

        var maxSteps = Math.ceil(1.0/DT);
        var corrs=[], times=[], stabT=maxSteps*DT;

        for(var step=1; step<=maxSteps; step++){
          var nH = new Array(n);
          for(var i=0;i<n;i++){
            var lh = deg[i]*H[i];
            for(var j=0;j<n;j++) lh -= adj[i][j]*H[j];
            nH[i] = H[i] - DT*lh;
          }
          corrs.push(spearman(nH, H));
          times.push(step*DT);
          H = nH;

          if(corrs.length >= STAB_WIN){
            var ok=true;
            for(var k=corrs.length-STAB_WIN; k<corrs.length-1; k++){
              if(Math.abs(corrs[k+1]-corrs[k])>=STAB_THRESH){ok=false;break;}
            }
            if(ok){ stabT=times[times.length-1]; break; }
          }
        }

        var result = [];
        for(var i=0;i<n;i++) {
          var nid = nodeIds[i];
          if(filterMet && metIds[nid] && !seeds.some(function(s){return s===nid;})) continue;
          result.push({id:nid, heat:H[i], group:realNodes[i].group, label:realNodes[i].label});
        }
        result.sort(function(a,b){return b.heat-a.heat;});

        return {scores:result, corrs:corrs, times:times, stabT:stabT, seedLabels:seeds.map(function(s){
          var nd=allNodes.find(function(nn){return nn.id===s;}); return nd?nd.label:s;
        })};
      }

      function drawCorrPlot(corrs, times, stabT, seedLabels) {
        var ctx = corrCanvas.getContext('2d');
        var W=corrCanvas.width, H=corrCanvas.height;
        ctx.clearRect(0,0,W,H);
        ctx.fillStyle='#fff'; ctx.fillRect(0,0,W,H);

        ctx.save();
        ctx.scale(4, 4);
        var logicalW = 280, logicalH = 140;

        var pad={t:25,r:30,b:35,l:40};
        var pw=logicalW-pad.l-pad.r, ph=logicalH-pad.t-pad.b;
        var minC=Math.min.apply(null,corrs), maxC=Math.max.apply(null,corrs);
        if(maxC===minC){minC-=0.01;maxC+=0.01;}
        var maxT=times[times.length-1];

        ctx.fillStyle='#0f172a'; ctx.font='bold 10px sans-serif'; ctx.textAlign='center';
        ctx.fillText('Seeds: '+seedLabels.join(', ').substring(0,40), logicalW/2, 12);

        ctx.strokeStyle='#94a3b8'; ctx.lineWidth=1;
        ctx.beginPath(); ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,pad.t+ph); ctx.lineTo(pad.l+pw,pad.t+ph); ctx.stroke();

        ctx.fillStyle='#475569'; ctx.font='9px sans-serif'; ctx.textAlign='center';
        ctx.fillText('Time Step',pad.l+pw/2,logicalH-3);
        ctx.save(); ctx.translate(10,pad.t+ph/2); ctx.rotate(-Math.PI/2); ctx.fillText('Spearman Corr',0,0); ctx.restore();

        ctx.font='8px sans-serif'; ctx.textAlign='right';
        for(var i=0;i<=4;i++){var v=minC+i*(maxC-minC)/4; var y=pad.t+ph-(v-minC)/(maxC-minC)*ph; ctx.fillText(v.toFixed(3),pad.l-3,y+3);}
        ctx.textAlign='center';
        for(var i=0;i<=4;i++){var v=i*maxT/4; var xx=pad.l+v/maxT*pw; ctx.fillText(v.toFixed(2),xx,pad.t+ph+12);}

        ctx.strokeStyle='#0f172a'; ctx.lineWidth=1.5; ctx.beginPath();
        for(var i=0;i<corrs.length;i++){
          var xx=pad.l+(times[i]/maxT)*pw, yy=pad.t+ph-(corrs[i]-minC)/(maxC-minC)*ph;
          if(i===0) ctx.moveTo(xx,yy); else ctx.lineTo(xx,yy);
        }
        ctx.stroke();

        var sx=pad.l+(stabT/maxT)*pw;
        ctx.strokeStyle='#D55E00'; ctx.lineWidth=1; ctx.setLineDash([4,3]);
        ctx.beginPath(); ctx.moveTo(sx,pad.t); ctx.lineTo(sx,pad.t+ph); ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle='#D55E00'; ctx.font='8px sans-serif'; ctx.textAlign='left';
        ctx.fillText('t='+stabT.toFixed(3), sx+2, pad.t+10);

        ctx.restore();
      }

      function applyVisualization(result, topPct, pal, topMicrobeCount) {
        var scores = result.scores;
        var count = Math.max(1, Math.ceil(scores.length * topPct / 100));
        var topScores = scores.slice(0, count);
        var topMap = {}; topScores.forEach(function(s){topMap[s.id]=s;});

        var minH = topScores[topScores.length-1] ? topScores[topScores.length-1].heat : 0;
        var maxH = topScores[0] ? topScores[0].heat : 1;

        // Highlight top microbes
        var microbeScores = scores.filter(function(s) { return s.group === 'Microbe'; });
        var topMicrobeCountVal = topMicrobeCount !== null && topMicrobeCount !== undefined && !isNaN(topMicrobeCount) ? topMicrobeCount : 5;
        var topMicrobeIds = {};
        microbeScores.slice(0, topMicrobeCountVal).forEach(function(s) {
          topMicrobeIds[s.id] = true;
        });

        var micN=[],pathN=[],metN2=[];
        topScores.forEach(function(s){
          if(s.group==='Microbe') micN.push(s);
          else if(s.group==='Pathway') pathN.push(s);
          else if(s.group==='Metabolite') metN2.push(s);
        });

        var rm=200+(micN.length*15), rp=150+(pathN.length*20), rme=100+(metN2.length*25);
        var xm=-(rm+rp+500), xp=0, xme2=(rp+rme+500);
        function assignC(arr,cx,r){ for(var i=0;i<arr.length;i++){ var a=(i/arr.length)*2*Math.PI; arr[i].nx=cx+r*Math.cos(a); arr[i].ny=r*Math.sin(a); }}
        assignC(micN,xm,rm); assignC(pathN,xp,rp); assignC(metN2,xme2,rme);

        var updates=[];
        allNodes.forEach(function(nd){
          if(nd.id.indexOf('LEG_')===0) return;
          if(topMap[nd.id]){
            var s=topMap[nd.id], c=getColor(s.heat,minH,maxH,pal);
            var isHighlighted = (s.group === 'Microbe' && topMicrobeIds[nd.id]);
            updates.push({
              id: nd.id, 
              hidden: false, 
              x: s.nx, 
              y: s.ny, 
              borderWidth: isHighlighted ? 4 : 1.5,
              color: {
                background: c, 
                border: isHighlighted ? '#22c55e' : '#475569', 
                highlight: {
                  background: c,
                  border: isHighlighted ? '#22c55e' : '#475569'
                }
              }
            });
          } else {
            updates.push({id:nd.id, hidden:true});
          }
        });
        nodesDS.update(updates);

        gradBar.style.background = getGrad(pal);
        minLbl.innerHTML = minH.toFixed(4);
        maxLbl.innerHTML = maxH.toFixed(4);

        if(visEngine && visEngine.fit) visEngine.fit({animation:{duration:1000,easingFunction:'easeInOutQuad'}});
      }

      // Event Handlers
      runBtn.onclick = function() {
        var selected = [];
        seedCheckboxes.forEach(function(cb) {
           if (cb.checked) selected.push(cb.value);
        });
        if(selected.length===0){ alert('Please select at least one seed node.'); return; }

        // Remove legend nodes
        if (nodesDS.get('LEG_MIC')) nodesDS.remove(['LEG_MIC', 'LEG_PATH', 'LEG_MET']);

        runBtn.innerHTML = 'Running...'; runBtn.disabled = true;
        setTimeout(function(){
          lastHeatResult = runLHD(selected, filterSel.value==='YES');
          drawCorrPlot(lastHeatResult.corrs, lastHeatResult.times, lastHeatResult.stabT, lastHeatResult.seedLabels);
          stabInfo.innerHTML = 'Stabilization at t = '+lastHeatResult.stabT.toFixed(4)+' | Seeds: '+lastHeatResult.seedLabels.join(', ');
          
          var val = parseInt(microbeInput.value);
          applyVisualization(lastHeatResult, parseInt(pctSlider.value), palSel.value, isNaN(val) ? null : val);
          
          postWrap.style.display = 'block';
          saveCsvBtn.style.display = 'block';
          savePlotBtn.style.display = 'block';
          runBtn.innerHTML = 'Run Analysis'; runBtn.disabled = false;
        }, 50);
      };

      pctSlider.oninput = function() {
        pctLabel.innerHTML = pctSlider.value + '%';
        if(lastHeatResult) {
          var val = parseInt(microbeInput.value);
          applyVisualization(lastHeatResult, parseInt(pctSlider.value), palSel.value, isNaN(val) ? null : val);
        }
      };
      
      microbeInput.oninput = function() {
        if(lastHeatResult) {
          var val = parseInt(microbeInput.value);
          applyVisualization(lastHeatResult, parseInt(pctSlider.value), palSel.value, isNaN(val) ? null : val);
        }
      };
      
      palSel.onchange = function() {
        if(lastHeatResult) {
          var val = parseInt(microbeInput.value);
          applyVisualization(lastHeatResult, parseInt(pctSlider.value), palSel.value, isNaN(val) ? null : val);
        }
      };

      saveCsvBtn.onclick = function() {
        if(!lastHeatResult) return;
        var csv = 'Node,Heat_score,Group\\n';
        lastHeatResult.scores.forEach(function(s){ csv += s.id+','+s.heat+','+s.group+'\\n'; });
        var blob = new Blob([csv], {type:'text/csv'});
        var a = document.createElement('a'); a.download='node_prior_heat_scores.csv'; a.href=URL.createObjectURL(blob); a.click();
      };

      savePlotBtn.onclick = function() {
        var a = document.createElement('a'); a.download='node_prior_stabilization_curve.png'; a.href=corrCanvas.toDataURL('image/png'); a.click();
      };

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
          nUpdates.push({id:n.id, hidden:false, color:null, borderWidth:1.5, x:ic?ic.x:0, y:ic?ic.y:0});
        });
        allEdges.forEach(function(e) {
          eUpdates.push({id:e.id, hidden:false, color:null, width:null});
        });
        nodesDS.update(nUpdates);
        edgesDS.update(eUpdates);
        postWrap.style.display = 'none';
        saveCsvBtn.style.display = 'none';
        savePlotBtn.style.display = 'none';
        microbeInput.value = '5';
        lastHeatResult = null;
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
    out_html <- file.path(output_directory, paste0("interactive_node_prior_", cleaned_input_file_name, ".html"))
    save_widget_safe(vis_plot, file = out_html, title = "NUIMM")
  }

  message("Node prioritization completed.")
  invisible(NULL)
}
