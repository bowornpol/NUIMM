#' Pathfinding using Dijkstra's algorithm
#'
#' @details
#' This function calculates the shortest path between two nodes using Dijkstra's
#' algorithm.
#'
#' @param multi_layered_network_file Path to input CSV/TSV.
#' @param source_node Name of the starting node.
#' @param target_node Name of the ending node.
#' @param output_directory Path to save results.
#' @param file_type "csv" or "tsv".
#' @param visualize Logical. If TRUE, saves PDF and PNG plots. Defaults to TRUE.
#' @export
find_path <- function(
  multi_layered_network_file,
  source_node,
  target_node,
  output_directory,
  file_type = c("csv", "tsv"),
  visualize = TRUE
) {
  file_type <- match.arg(file_type)
  message("Starting shortest pathfinding using Dijkstra's algorithm.")

  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)

  # 1. Load network
  message("\n1. Loading network...")
  if (file_type == "csv") {
    network_data <- read.csv(multi_layered_network_file, stringsAsFactors = FALSE)
  } else {
    network_data <- read.delim(multi_layered_network_file, stringsAsFactors = FALSE)
  }

  # Ensure columns exist
  req <- c("Feature1", "Feature2", "edge_score")
  if (!all(req %in% colnames(network_data))) stop("Missing columns: ", paste(req, collapse=", "))

  # Prepare numeric weights
  network_data$edge_score <- as.numeric(network_data$edge_score)
  network_data$edge_score_abs <- abs(network_data$edge_score)

  # 2. Create graph
  g <- igraph::graph_from_data_frame(d = network_data, directed = FALSE)

  # Transform weights (Higher score = Shorter path)
  igraph::E(g)$weight <- sapply(igraph::E(g)$edge_score_abs, function(w) {
    if (is.na(w)) Inf else if (w < 1) 1/w else 1/(w + 0.1)
  })

  # 3. Check nodes
  if (!source_node %in% igraph::V(g)$name) stop("Source node not found.")
  if (!target_node %in% igraph::V(g)$name) stop("Target node not found.")

  # 4. Find Path
  message("\n2. Finding shortest path...")
  result <- igraph::shortest_paths(g, from = source_node, to = target_node, weights = igraph::E(g)$weight, output = "both")
  path_vertices <- result$vpath[[1]]
  path_edges <- result$epath[[1]]

  if (length(path_vertices) > 1) {
    message("  Path found with ", length(path_edges), " steps.")

    # Extract Ordered Nodes
    node_names_ordered <- names(path_vertices)

    # Create Coordinates Dataframe for Visualization (Straight Line)
    # x = 1, 2, 3... (Sequence)
    # y = 0 (Flat line)
    coords <- data.frame(
      name = node_names_ordered,
      x = seq_along(node_names_ordered),
      y = rep(0, length(node_names_ordered))
    )

    # Assign Node Categories
    coords$type <- sapply(coords$name, function(x) {
      if (grepl("d__|p__|c__|o__|f__|g__|s__|Bacteria", x)) return("Microbe")
      if (grepl("^ko[0-9]+", x)) return("Pathway")
      return("Metabolite")
    })

    # Set Factor levels to control Legend Order
    coords$type <- factor(coords$type, levels = c("Microbe", "Pathway", "Metabolite"))

    # Prepare Edge Coordinates
    edges_plot <- data.frame(
      x = coords$x[1:(nrow(coords)-1)],
      y = coords$y[1:(nrow(coords)-1)],
      xend = coords$x[2:nrow(coords)],
      yend = coords$y[2:nrow(coords)]
    )

    # Save CSV
    edge_df <- igraph::as_data_frame(g, what = "edges")[path_edges, ]
    path_df <- dplyr::select(edge_df, Source = from, Target = to, edge_score = edge_score)

    safe_from <- gsub("[^A-Za-z0-9]", "_", source_node)
    safe_to <- gsub("[^A-Za-z0-9]", "_", target_node)
    base_name <- paste0("path_", safe_from, "_to_", safe_to)
    write.csv(path_df, file.path(output_directory, paste0(base_name, ".csv")), row.names = FALSE)

    # 5. Shortest path visualization
    if (visualize) {
      message("\n3. Generating network visualization...")
      if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install 'ggplot2'.")

      p <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = edges_plot,
                              ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                              color = "gray70", size = 1.2, alpha = 0.8) +
        ggplot2::geom_point(data = coords,
                            ggplot2::aes(x = x, y = y, color = type),
                            size = 12) +
        ggplot2::scale_color_manual(
          name = "Node Type", # Legend Title
          values = c(
            "Microbe" = "#C1ABAD",
            "Pathway" = "#9AA374",
            "Metabolite" = "#4E7286"
          )
        ) +
        ggplot2::geom_text(data = coords,
                           ggplot2::aes(x = x, y = y, label = name),
                           vjust = -2.0, # Push text up
                           hjust = 0.5,
                           fontface = "bold",
                           size = 4) +
        ggplot2::theme_void() +
        ggplot2::theme(
          legend.position = "bottom",
          legend.title = ggplot2::element_text(face = "bold"),
          plot.margin = ggplot2::unit(c(2, 1, 1, 1), "cm") # Top margin for labels
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::expand_limits(y = c(-1, 1), x = c(0.5, max(coords$x) + 0.5))

      # Save
      pdf_path <- file.path(output_directory, paste0(base_name, ".pdf"))
      png_path <- file.path(output_directory, paste0(base_name, ".png"))
      ggplot2::ggsave(pdf_path, plot = p, width = 12, height = 6)
      ggplot2::ggsave(png_path, plot = p, width = 12, height = 6, dpi = 600, bg = "white")
      message("  Saved visualization to: ", png_path)
    }

    return(invisible(NULL))
  } else {
    message("  No path found.")
    return(invisible(NULL))
  }
}
