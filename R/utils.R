# R/utils.R

#' Internal helper function to read CSV or TSV files
#'
#' @param file_path A character string specifying the path to the input file.
#' @param file_type A character string indicating the file type ("csv" or "tsv").
#' @param ... Additional arguments to pass to the reader.
#' @return A data frame containing the data from the specified file.
#' @keywords internal
#' @noRd
read_input_file <- function(file_path, file_type = NULL, ...) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }

  # Capture extra arguments
  args <- list(...)

  # Extract and remove row.names
  row_col <- args$row.names
  args$row.names <- NULL

  if (requireNamespace("data.table", quietly = TRUE)) {
    # Fast read using data.table
    data <- do.call(data.table::fread, c(list(file = file_path, data.table = FALSE), args))
  } else {
    # Fallback to base R and guess separator
    ext <- tolower(tools::file_ext(file_path))
    sep_char <- if (ext == "tsv" || (!is.null(file_type) && file_type == "tsv")) "\t" else ","
    data <- read.delim(file_path, sep = sep_char, ...)
  }

  # Manually apply row names
  if (!is.null(row_col) && row_col == 1) {
    rownames(data) <- as.character(data[[1]])
    data <- data[, -1, drop = FALSE]
  }

  return(data)
}

#' Get JavaScript for Ctrl-drag Multi-Node Selection
#' @keywords internal
#' @noRd
get_ctrl_drag_js <- function() {
  "
  // Ctrl-Drag Multi-Node Selection
  (function() {
    var network = this.network || (typeof widget !== 'undefined' ? widget.network : null);
    if (!network) return;

    var isDragging = false;
    var startX, startY;
    var selectionBox = null;
    var ctrlPressed = false;

    window.addEventListener('keydown', function(e) {
      if (e.key === 'Control') ctrlPressed = true;
    });
    window.addEventListener('keyup', function(e) {
      if (e.key === 'Control') ctrlPressed = false;
    });

    el.addEventListener('mousedown', function(e) {
      if (e.ctrlKey || ctrlPressed) {
        isDragging = true;
        startX = e.pageX;
        startY = e.pageY;

        if (!selectionBox) {
          selectionBox = document.createElement('div');
          selectionBox.style.cssText = 'position:absolute;border:1px dashed #3b82f6;background:rgba(59,130,246,0.15);pointer-events:none;z-index:999999;';
          document.body.appendChild(selectionBox);
        }

        selectionBox.style.left = startX + 'px';
        selectionBox.style.top = startY + 'px';
        selectionBox.style.width = '0px';
        selectionBox.style.height = '0px';
        selectionBox.style.display = 'block';

        network.setOptions({ interaction: { dragView: false } });
        e.stopPropagation();
      }
    });

    window.addEventListener('mousemove', function(e) {
      if (isDragging && selectionBox) {
        var currentX = e.pageX;
        var currentY = e.pageY;

        var left = Math.min(startX, currentX);
        var top = Math.min(startY, currentY);
        var width = Math.abs(startX - currentX);
        var height = Math.abs(startY - currentY);

        selectionBox.style.left = left + 'px';
        selectionBox.style.top = top + 'px';
        selectionBox.style.width = width + 'px';
        selectionBox.style.height = height + 'px';
      }
    });

    window.addEventListener('mouseup', function(e) {
      if (isDragging) {
        isDragging = false;
        if (selectionBox) selectionBox.style.display = 'none';

        network.setOptions({ interaction: { dragView: true } });

        var rect = el.getBoundingClientRect();
        var scrollLeft = window.pageXOffset || document.documentElement.scrollLeft;
        var scrollTop = window.pageYOffset || document.documentElement.scrollTop;

        var domStartX = startX - rect.left - scrollLeft;
        var domStartY = startY - rect.top - scrollTop;
        var domEndX = e.pageX - rect.left - scrollLeft;
        var domEndY = e.pageY - rect.top - scrollTop;

        var canvasStart = network.DOMtoCanvas({ x: domStartX, y: domStartY });
        var canvasEnd = network.DOMtoCanvas({ x: domEndX, y: domEndY });

        var minX = Math.min(canvasStart.x, canvasEnd.x);
        var maxX = Math.max(canvasStart.x, canvasEnd.x);
        var minY = Math.min(canvasStart.y, canvasEnd.y);
        var maxY = Math.max(canvasStart.y, canvasEnd.y);

        var allNodes = network.body.nodeIndices;
        var nodesToSelect = [];
        allNodes.forEach(function(nodeId) {
          if (typeof nodeId === 'string' && nodeId.indexOf('LEG_') === 0) return;
          var pos = network.getPosition(nodeId);
          if (pos.x >= minX && pos.x <= maxX && pos.y >= minY && pos.y <= maxY) {
            nodesToSelect.push(nodeId);
          }
        });

        if (nodesToSelect.length > 0) {
          network.selectNodes(nodesToSelect);
        } else {
          network.unselectAll();
        }
      }
    });
  }).call(this);
  "
}

#' Determine node groups from network data and falls back to regex
#' @keywords internal
#' @noRd
determine_node_groups <- function(nodes, network_data, source_col, target_col) {
  groups <- rep("Metabolite", length(nodes))
  names(groups) <- nodes

  type_col <- if ("type" %in% colnames(network_data)) "type" else if ("edge_type" %in% colnames(network_data)) "edge_type" else NULL

  if (!is.null(type_col)) {
    microbe_pathway_edges <- network_data[network_data[[type_col]] == "Microbe-Pathway", ]
    pathway_metabolite_edges <- network_data[network_data[[type_col]] == "Pathway-Metabolite", ]

    microbe_nodes <- unique(microbe_pathway_edges[[source_col]])
    metabolite_nodes <- unique(pathway_metabolite_edges[[target_col]])

    for (node in nodes) {
      if (node %in% microbe_nodes) {
        groups[node] <- "Microbe"
      } else if (node %in% metabolite_nodes) {
        groups[node] <- "Metabolite"
      } else {
        groups[node] <- "Pathway"
      }
    }
    return(groups)
  }

  for (node in nodes) {
    if (grepl("d__|p__|c__|o__|f__|g__|s__|Bacteria", node)) {
      groups[node] <- "Microbe"
    } else if (grepl("ko[0-9]+|PATH|Pwy|pwy|Pathway|pathway|degradation|biosynthesis|metabolism|synthesis", node, ignore.case = TRUE)) {
      groups[node] <- "Pathway"
    } else {
      groups[node] <- "Metabolite"
    }
  }

  return(groups)
}

#' Save an htmlwidget safely, falling back to non-selfcontained if pandoc is missing
#' @param widget The htmlwidget to save.
#' @param file Output file path.
#' @param title HTML page title.
#' @keywords internal
#' @noRd
save_widget_safe <- function(widget, file, title = "NUIMM") {
  has_pandoc <- tryCatch({
    info <- rmarkdown::find_pandoc()
    !is.null(info$dir) && nzchar(info$dir)
  }, error = function(e) FALSE)

  if (has_pandoc) {
    htmlwidgets::saveWidget(widget, file = file, selfcontained = TRUE, title = title)
  } else {
    htmlwidgets::saveWidget(widget, file = file, selfcontained = FALSE, title = title)
    message("Pandoc unavailable. Exporting HTML with external dependencies.")
  }
}

