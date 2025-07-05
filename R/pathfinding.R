#' Pathfinding using Dijkstra's algorithm
#'
#' This function identifies the shortest path between a specified source node
#' and a target node within a multi-layered network using Dijkstra's algorithm.
#'
#' @details
#' The edge weights are transformed from `edge_score` as follows:
#' - `edge_score` is first converted to its absolute value to ensure non-negative weights.
#' - If the absolute `edge_score` is `NA`, weight is `Inf`.
#' - If the absolute `edge_score` is less than 1, weight is `1 / absolute_edge_score`.
#' - If the absolute `edge_score` is exactly 1, weight is `1 / (absolute_edge_score + 0.1)` to avoid `1/1=1` which might not be ideal for shortest path.
#' - If the absolute `edge_score` is greater than 1, weight is `absolute_edge_score`.
#' This transformation aims to represent stronger connections (higher `edge_score`)
#' as shorter paths (lower weights), ensuring all weights are non-negative for Dijkstra's algorithm.
#'
#' @param multi_layered_network_file A character string specifying the path to the
#'    integrated multi-layered network data file (e.g., output
#'    from `construct_multi_layered_network`). Expected columns: 'Feature1',
#'    'Feature2', 'edge_score' (numeric), and 'edge_type'.
#' @param source_node A character string specifying the name of the starting node
#'    for pathfinding. This node must exist in the network.
#' @param target_node A character string specifying the name of the ending node
#'    for pathfinding. This node must exist in the network.
#' @param output_directory A character string specifying the path to the directory
#'    where the shortest path results will be saved as a CSV file. The directory
#'    will be created if it does not exist.
#' @param file_type A character string indicating the type of input file.
#'    Must be "csv" (for comma-separated) or "tsv" (for tab-separated).
#' @return The function's primary output is a CSV file saved
#'    to the specified `output_directory`, detailing the edges of the shortest path
#'    found between the source and target nodes. If no path is found, `NULL` is returned
#'    invisibly and a message is printed.
#' @references
#' Dijkstra EW. A note on two problems in connexion with graphs.  Edsger Wybe Dijkstra: his life, work, and legacy2022. p. 287-90.
#' @export
pathfinding <- function(
  multi_layered_network_file,
  source_node,
  target_node,
  output_directory,
  file_type = c("csv", "tsv")
) {
  file_type <- match.arg(file_type)
  message("Starting shortest pathfinding using Dijkstra's algorithm.")

  # Create output directory if it doesn't exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    message("Created output directory: ", output_directory)
  } else {
    message("Output directory already exists: ", output_directory)
  }

  # 1. Load network
  message("\n1. Loading multi-layered network from: ", multi_layered_network_file)

  required_cols <- c("Feature1", "Feature2", "edge_score", "edge_type")
  if (!file.exists(multi_layered_network_file)) {
    stop("Network file not found: '", multi_layered_network_file, "'")
  }

  network_data <- read_input_file(multi_layered_network_file, file_type = file_type, stringsAsFactors = FALSE)
  if (!all(required_cols %in% colnames(network_data))) {
    stop("Network file must contain columns: ", paste(required_cols, collapse = ", "))
  }
  message("  Successfully loaded and validated network file.")

  # Ensure edge_score is numeric
  network_data$edge_score <- as.numeric(network_data$edge_score)

  # --- FIX: Convert edge_score to absolute value before weight transformation ---
  network_data$edge_score_abs <- abs(network_data$edge_score)
  message("  Converted 'edge_score' to its absolute value for weight calculation.")

  # 2. Create graph directly with edge attributes
  message("  Creating graph and assigning weights.")
  g <- igraph::graph_from_data_frame(d = network_data, directed = FALSE)

  # Apply weight transformation using the *absolute* edge_score
  igraph::E(g)$weight <- sapply(igraph::E(g)$edge_score_abs, function(w) { # Changed to edge_score_abs
    if (is.na(w)) {
      Inf
    } else if (w < 1) {
      1 / w
    } else if (w == 1) {
      1 / (w + 0.1) # Add a small epsilon to avoid 1/1=1 if 1 is a special case
    } else {
      w
    }
  })

  # 3. Validate nodes
  message("\n2. Validating source and target nodes against the network.")
  if (!source_node %in% igraph::V(g)$name) stop("Source node not found in network: ", source_node)
  if (!target_node %in% igraph::V(g)$name) stop("Target node not found in network: ", target_node)
  message("  Both source and target nodes found.")

  # 4. Find shortest path
  message("\n3. Finding shortest path from '", source_node, "' to '", target_node, "'...")
  result <- igraph::shortest_paths(g, from = source_node, to = target_node, weights = igraph::E(g)$weight, output = "both")

  path_vertices <- result$vpath[[1]]
  path_edges <- result$epath[[1]]

  if (length(path_vertices) > 1 && length(path_edges) > 0) {
    message("  Path found with ", length(path_edges), " steps.")

    edge_df <- igraph::as_data_frame(g, what = "edges")[path_edges, ]
    path_df <- dplyr::select(edge_df,
                             Source = from,
                             Target = to,
                             original_edge_score = edge_score,       # This refers to the original, potentially negative, edge_score
                             transformed_weight = weight,
                             edge_type = edge_type)

    # Sanitize node names for filename
    safe_from <- gsub("[^A-Za-z0-9_.-]", "_", source_node)
    safe_to <- gsub("[^A-Za-z0-9_.-]", "_", target_node)
    output_file <- file.path(output_directory, paste0("path_", safe_from, "_to_", safe_to, ".csv"))

    write.csv(path_df, output_file, row.names = FALSE)
    message("  Saved path to: ", output_file)

    message("Pathfinding complete.")
    return(invisible(NULL))
  } else {
    message("  No finite path found between source and target.")
    message("Pathfinding complete (no path found).")
    return(invisible(NULL))
  }
}
