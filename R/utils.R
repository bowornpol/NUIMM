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
    if (!"fill" %in% names(args)) args$fill <- Inf
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

#' Strip group suffix from sample IDs (regex-safe)
#'
#' Uses an anchored, escaped pattern so group names containing regex
#' metacharacters (e.g. ".", "+") are matched literally.
#'
#' @param ids Character vector of sample IDs.
#' @param group Character string: the group label to strip.
#' @return Character vector with the trailing `_<group>` removed.
#' @keywords internal
#' @noRd
strip_group_suffix <- function(ids, group) {
  # Escape regex metacharacters so the group name is matched literally
  meta_chars <- c("\\", "[", "]", "(", ")", "{", "}", ".", "*", "+", "?", "^", "$", "|")
  esc <- group
  for (ch in meta_chars) {
    esc <- gsub(ch, paste0("\\", ch), esc, fixed = TRUE)
  }
  sub(paste0("_", esc, "$"), "", ids)
}

#' Derive pairwise comparisons from metadata or use user-supplied list
#'
#' Centralizes the repeated pattern of `sort(unique(class)); combn(...,2)`.
#'
#' @param meta Data frame with a `class` column.
#' @param comparisons_list Optional user-supplied list of comparisons.
#' @return A list of length-2 character vectors.
#' @keywords internal
#' @noRd
derive_comparisons <- function(meta, comparisons_list = NULL) {
  if (!is.null(comparisons_list)) {
    return(comparisons_list)
  }
  conditions <- sort(unique(meta$class))
  combn(conditions, 2, simplify = FALSE)
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

  warning("No 'type' or 'edge_type' column found in network data; falling back to regex-based node classification. ",
    "This heuristic may misclassify nodes whose names contain pathway-related keywords ",
    "(e.g. 'degradation', 'biosynthesis'). Consider adding a 'type' column.",
    call. = FALSE
  )

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
  has_pandoc <- tryCatch(
    {
      info <- rmarkdown::find_pandoc()
      !is.null(info$dir) && nzchar(info$dir)
    },
    error = function(e) FALSE
  )

  if (has_pandoc) {
    htmlwidgets::saveWidget(widget, file = file, selfcontained = TRUE, title = title)
  } else {
    htmlwidgets::saveWidget(widget, file = file, selfcontained = FALSE, title = title)
    message("Pandoc unavailable. Exporting HTML with external dependencies.")
  }
}

#' Compute three-cluster circular layout for nodes by group
#' @keywords internal
#' @noRd
compute_circular_layout <- function(nodes_df) {
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

  place_circle <- function(idx, cx, r) {
    if (length(idx) > 0) {
      ang <- seq(0, 2 * pi, length.out = length(idx) + 1)[seq_along(idx)]
      nodes_df$x[idx] <<- cx + r * cos(ang)
      nodes_df$y[idx] <<- r * sin(ang)
    }
  }

  place_circle(idx_mic, x_mic, r_mic)
  place_circle(idx_path, x_path, r_path)
  place_circle(idx_met, x_met, r_met)

  return(nodes_df)
}

#' Add legend nodes to a nodes data.frame, auto-adapting to its columns
#' @keywords internal
#' @noRd
add_legend_nodes <- function(nodes_df) {
  max_y <- max(nodes_df$y, na.rm = TRUE)
  legend_y <- max_y + 400
  legend_list <- list()
  for (col in colnames(nodes_df)) {
    legend_list[[col]] <- switch(col,
      id = c("LEG_MIC", "LEG_PATH", "LEG_MET"),
      label = c("Microbe", "Pathway", "Metabolite"),
      group = c("Microbe", "Pathway", "Metabolite"),
      size = c(60, 60, 60),
      x = c(-300, 0, 300),
      y = rep(legend_y, 3),
      title = c("", "", ""),
      {
        if (is.numeric(nodes_df[[col]])) rep(0, 3) else rep("", 3)
      }
    )
  }
  legend_df <- as.data.frame(legend_list, stringsAsFactors = FALSE)
  rbind(nodes_df, legend_df)
}

#' Validate comparisons_list structure and optionally check group existence
#' @param comparisons_list The list to validate.
#' @param metadata Optional data frame with a `class` column for group-existence checks.
#' @keywords internal
#' @noRd
validate_comparisons_structure <- function(comparisons_list, metadata = NULL) {
  if (is.null(comparisons_list)) {
    return(invisible(NULL))
  }
  if (!is.list(comparisons_list)) {
    stop("'comparisons_list' must be a list of character vectors, e.g., list(c('V1', 'V3')). Got: ", class(comparisons_list)[1])
  }
  for (i in seq_along(comparisons_list)) {
    comp <- comparisons_list[[i]]
    if (!is.character(comp) || length(comp) != 2) {
      stop(sprintf("comparisons_list[[%d]] must be a character vector of length 2.", i))
    }
  }
  # Check that specified group names actually exist in metadata
  if (!is.null(metadata) && "class" %in% colnames(metadata)) {
    available_groups <- unique(metadata$class)
    requested_groups <- unique(unlist(comparisons_list))
    missing <- setdiff(requested_groups, available_groups)
    if (length(missing) > 0) {
      stop(sprintf(
        "comparisons_list references groups not found in metadata$class: %s. Available groups: %s",
        paste(missing, collapse = ", "),
        paste(available_groups, collapse = ", ")
      ))
    }
  }
  invisible(NULL)
}
