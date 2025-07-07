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
#'    - It sums `factorial(clique_size - 1)` for each maximal clique of size
#'      3 or more that contains the node.
#'    - After this summation, if a node's calculated MCC score is still zero
#'      (meaning it wasn't part of any clique of size >= 3) AND it has a
#'      degree greater than zero, its MCC is then set to its degree.
#' 5. Ranks nodes by their MCC score in descending order.
#' 6. Optionally filters to return only the top N hub nodes.
#' 7. Saves the ranked hub identification results to a CSV file.
#'
#' @param multi_layered_network_file A character string specifying the path to the
#'    integrated multi-layered network data file (output
#'    from `con.mln` or `con.mln.all`). Expected columns: 'Feature1' and
#'    'Feature2'. 'edge_score' and 'edge_type' are loaded but not directly used
#'    in the MCC calculation itself (only edge presence matters for cliques).
#' @param output_directory A character string specifying the path to the directory
#'    where the hub identification results will be saved as a CSV file.
#'    The directory will be created if it does not exist.
#' @param file_type A character string indicating the type of input file.
#'    Must be "csv" (for comma-separated) or "tsv" (for tab-separated).
#' @param top_n_hubs An optional integer specifying the number of top hub nodes
#'    to return. If `NULL` (default), all nodes will be returned, ranked by MCC score.
#' @return The function's primary output is a CSV file saved
#'    to the specified `output_directory`, containing the ranked list of nodes
#'    and their MCC scores.
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
  message("  Successfully loaded and validated network file.")

  # 2. Create an undirected graph
  message("  Creating undirected graph for MCC calculation.")
  # MCC is a structural centrality, typically applied to unweighted, undirected graphs.
  # We only need Feature1 and Feature2 for graph structure.
  g <- igraph::graph_from_data_frame(d = network_data[, required_cols], directed = FALSE)
  message("  Graph created with ", igraph::vcount(g), " nodes and ", igraph::ecount(g), " edges.")

  # 3. Find all maximal cliques
  message("\n2. Finding all maximal cliques. This may take a while for large or dense networks...")
  cliques <- tryCatch(
    igraph::max_cliques(g),
    error = function(e) {
      stop(paste("Error finding maximal cliques: ", e$message, ". This can be computationally intensive for large graphs.", sep = ""))
    }
  )
  message("  Found ", length(cliques), " maximal cliques.")

  # 4. Calculate MCC scores for each node based on the provided formula
  message("\n3. Calculating MCC scores for each node.")
  # Initialize MCC scores to 0 for all nodes
  mcc_scores <- setNames(numeric(igraph::vcount(g)), igraph::V(g)$name)

  # Sum (clique_size - 1)! for all maximal cliques of size >= 3
  for (clique_nodes_indices in cliques) {
    clique_size <- length(clique_nodes_indices)
    if (clique_size >= 3) { # Only cliques of size 3 or more contribute to the sum for the primary formula
      clique_score <- factorial(clique_size - 1)
      for (node_index in clique_nodes_indices) {
        node_name <- igraph::V(g)$name[node_index]
        mcc_scores[node_name] <- mcc_scores[node_name] + clique_score
      }
    }
  }

  # Apply the special case: if MCC is 0 (from cliques >=3) but degree > 0, MCC = degree
  node_degrees <- igraph::degree(g)
  for (node_name in igraph::V(g)$name) {
    if (mcc_scores[node_name] == 0 && node_degrees[node_name] > 0) {
      mcc_scores[node_name] <- node_degrees[node_name]
      message("  Node '", node_name, "' not in cliques of size >=3; MCC set to its degree (", node_degrees[node_name], ").")
    }
  }
  # Note: Isolated nodes (degree 0) will correctly retain MCC of 0.

  # 5. Create a data frame of results and rank
  message("\n4. Ranking nodes by MCC score.")
  hub_results_df <- dplyr::arrange(
    data.frame(
      Node = names(mcc_scores),
      MCC_score = mcc_scores, # Changed MCC_Score to MCC_score (PascalCase to CamelCase)
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
      message("  Filtered to top ", top_n_hubs, " hub nodes.")
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
