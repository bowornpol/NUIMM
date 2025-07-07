#' Hub identification using Maximal Clique Centrality (MCC) algorithm
#'
#' This function identifies hub nodes in a multi-layered network using the
#' Maximal Clique Centrality (MCC) algorithm.
#'
#' @details
#' The function performs the following steps:
#' 1. Loads the integrated multi-layered network from the specified CSV/TSV file.
#' 2. Constructs an undirected `igraph` object from the network data.
#' 3. Finds all maximal cliques within the graph.
#' 4. Calculates the MCC score for each node:
#'    - It sums `factorial(clique_size - 1)` for each maximal clique that contains the node.
#'      (Note: This means cliques of any size will contribute to the sum).
#'    - If there is no edge between the neighbors of a node, its MCC is then set to its degree.
#'      This rule applies if the node has degree 0 (MCC becomes 0), degree 1 (MCC becomes 1),
#'      or degree >= 2 with a local clustering coefficient of 0.
#' 5. Ranks nodes by their MCC score in descending order.
#' 6. Optionally filters to return only the top N hub nodes.
#' 7. Saves the ranked hub identification results to a CSV file.
#'
#' @param multi_layered_network_file A character string specifying the path to the
#'    integrated multi-layered network data file (output
#'    from `con.mln` or `con.mln.all`). Expected columns: 'Feature1' and
#'    'Feature2'. 'edge_score' and 'edge_type' are loaded but not directly used
#'    in the MCC calculation itself (only edge presence matters for cliques).
#' @param output_directory A character string specifying the path to the directory
#'    where the hub identification results will be saved as a CSV file.
#'    The directory will be created if it does not exist.
#' @param file_type A character string indicating the type of input file.
#'    Must be "csv" (for comma-separated) or "tsv" (for tab-separated).
#' @param top_n_hubs An optional integer specifying the number of top hub nodes
#'    to return. If `NULL` (default), all nodes will be returned, ranked by MCC score.
#' @return The function's primary output is a CSV file saved
#'    to the specified `output_directory`, containing the ranked list of nodes
#'    and their MCC scores.
#' @references
#' Chin CH, Chen SH, Wu HH, Ho CW, Ko MT, Lin CY. cytoHubba: identifying hub objects and sub-networks from complex interactome. BMC systems biology. 2014 Dec;8:1-7.
#' @export
iden_hub <- function(
  multi_layered_network_file,
  output_directory,
  file_type = c("csv", "tsv"),
  top_n_hubs = NULL
) {
  file_type <- match.arg(file_type)
  message("Starting hub identification using Maximal Clique Centrality (MCC) algorithm.")

  # Create output directory if it doesn't exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    message("Created output directory: ", output_directory)
  } else {
    message("Output directory already exists: ", output_directory)
  }

  # Extract base name from input file for output specificity
  input_file_base_name <- tools::file_path_sans_ext(basename(multi_layered_network_file))
  # Clean the base name for use in filenames (e.g., replace non-alphanumeric with underscore)
  cleaned_input_file_name <- gsub("[^A-Za-z0-9_]", "", input_file_base_name)

  # 1. Load the multi-layered network file
  message("\n1. Loading multi-layered network from: ", multi_layered_network_file)
  required_cols <- c("Feature1", "Feature2") # Only these are strictly needed for graph structure
  if (!file.exists(multi_layered_network_file)) {
    stop("Network file not found: '", multi_layered_network_file, "'")
  }

  network_data <- read_input_file(multi_layered_network_file, file_type = file_type, stringsAsFactors = FALSE)
  if (!all(required_cols %in% colnames(network_data))) {
    stop("Network file must contain columns: ", paste(required_cols, collapse = ", "))
  }
  message("  Successfully loaded and validated network file.")

  # 2. Create an undirected graph
  message("  Creating undirected graph for MCC calculation.")
  g <- igraph::graph_from_data_frame(d = network_data[, required_cols], directed = FALSE)
  message("  Graph created with ", igraph::vcount(g), " nodes and ", igraph::ecount(g), " edges.")

  # 3. Find all maximal cliques
  message("\n2. Finding all maximal cliques. This may take a while for large or dense networks...")
  cliques <- tryCatch(
    igraph::max_cliques(g),
    error = function(e) {
      stop(paste("Error finding maximal cliques: ", e$message, ". This can be computationally intensive for large graphs.", sep = ""))
    }
  )
  message("  Found ", length(cliques), " maximal cliques.")

  # 4. Calculate MCC scores for each node based on the provided formula
  message("\n3. Calculating MCC scores for each node.")
  # Initialize MCC scores to 0 for all nodes
  mcc_scores <- setNames(numeric(igraph::vcount(g)), igraph::V(g)$name)

  # Sum (clique_size - 1)! for all maximal cliques containing the node
  # Removed the 'clique_size >= 3' filter as per exact definition
  for (clique_nodes_indices in cliques) {
    clique_size <- length(clique_nodes_indices)
    clique_score <- factorial(clique_size - 1) # factorial(0) = 1, factorial(1) = 1
    for (node_index in clique_nodes_indices) {
      node_name <- igraph::V(g)$name[node_index]
      mcc_scores[node_name] <- mcc_scores[node_name] + clique_score
    }
  }

  # Pre-calculate degrees and local clustering coefficients
  node_degrees <- igraph::degree(g)
  # transitivity(type="local") returns NaN for nodes with degree < 2.
  # We will handle these cases based on the definition "no edge between neighbors".
  node_clustering_coeffs <- igraph::transitivity(g, type = "local", vids = igraph::V(g))
  names(node_clustering_coeffs) <- igraph::V(g)$name

  # Apply the special case: If there is no edge between the neighbors of the node v, MCC(v) = degree(v)
  for (node_name in igraph::V(g)$name) {
    current_degree <- node_degrees[node_name]

    # Determine if "no edge between neighbors" condition is met
    is_no_edge_between_neighbors <- FALSE
    if (current_degree == 0) {
      # Isolated node: no neighbors, so no edges between them. MCC should be 0 (its degree).
      is_no_edge_between_neighbors <- TRUE
    } else if (current_degree == 1) {
      # Node with one neighbor: trivially, no edges between multiple neighbors. MCC should be 1 (its degree).
      is_no_edge_between_neighbors <- TRUE
    } else { # current_degree >= 2
      # For nodes with 2 or more neighbors, check if local clustering coefficient is 0
      if (!is.na(node_clustering_coeffs[node_name]) && node_clustering_coeffs[node_name] == 0) {
        is_no_edge_between_neighbors <- TRUE
      }
    }

    # If the special condition is met, override the MCC score with the node's degree
    if (is_no_edge_between_neighbors) {
      # Only log if the score is actually changing from the sum calculation.
      # For degree 0 and 1, the summation might already result in degree, so no change.
      if (mcc_scores[node_name] != current_degree) {
        message("  Node '", node_name, "' has no edge between neighbors (or degree < 2); MCC overridden to its degree (", current_degree, ").")
      }
      mcc_scores[node_name] <- current_degree
    }
  }

  # 5. Create a data frame of results and rank
  message("\n4. Ranking nodes by MCC score.")
  hub_results_df <- dplyr::arrange(
    data.frame(
      Node = names(mcc_scores),
      MCC_score = mcc_scores,
      stringsAsFactors = FALSE
    ),
    dplyr::desc(MCC_score) # Sort in descending order of MCC_score
  )

  # 6. Apply top_n_hubs filter if specified
  if (!is.null(top_n_hubs) && is.numeric(top_n_hubs) && top_n_hubs > 0) {
    if (top_n_hubs > nrow(hub_results_df)) {
      warning("Requested top_n_hubs (", top_n_hubs, ") is greater than total nodes (", nrow(hub_results_df), "). Returning all nodes.")
    } else {
      hub_results_df <- head(hub_results_df, n = top_n_hubs)
      message("  Filtered to top ", top_n_hubs, " hub nodes.")
    }
  }

  # 7. Save the results
  output_filename <- paste0("hub_mcc_", cleaned_input_file_name,
                            if(!is.null(top_n_hubs)) paste0("_top", top_n_hubs) else "", ".csv")
  output_filepath <- file.path(output_directory, output_filename)
  message("\n5. Saving hub identification results to: ", output_filepath)
  write.csv(hub_results_df, output_filepath, row.names = FALSE)
  message("Hub identification complete.")

  return(invisible(NULL))
}
