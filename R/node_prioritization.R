#' Node Prioritization using Laplacian Heat Diffusion (LHD) algorithm
#'
#' This function implements Laplacian Heat Diffusion (LHD) algorithm
#' on a multi-layered network to prioritize nodes based on their
#' association from specific metabolite seed nodes.
#'
#' @details
#' The process involves:
#' 1. Loading the integrated multi-layered network data.
#' 2. Identifying all unique metabolite nodes, which serve as seed points for diffusion.
#' 3. For each metabolite seed:
#'    a. Optionally filters the network to exclude all other metabolite nodes and their connected edges
#'       if `filter_other_metabolite_edges` is `TRUE`.
#'    b. Converts edge scores to absolute values and constructs an `igraph` object.
#'    c. Computes the graph Laplacian matrix.
#'    d. Initializes a heat vector with heat concentrated on the current seed node.
#'    e. Simulates heat diffusion over time steps using matrix exponentials.
#'    f. Determines the stabilization time by monitoring the Spearman correlation
#'       between successive heat vectors.
#'    g. Calculates final heat scores for all nodes at the stabilization time.
#'    h. Saves the heat scores and the correlation-over-time data to CSV files.
#'    i. Generates and saves a ggplot2 visualization of the Spearman correlation
#'       over time, indicating the stabilization point.
#'
#' @param multi_layered_network_file A character string specifying the path to the
#'    integrated multi-layered network data file (e.g., output
#'    from `construct_multi_layered_network`). Expected columns: 'Feature1',
#'    'Feature2', 'edge_score' (numeric), and 'edge_type'. # Changed Edge_Score and Edge_Type
#' @param output_directory A character string specifying the path to the directory
#'    where the output CSV files (heat scores, correlation data) and plots (correlation
#'    plots) will be saved. The directory will be created if it does not exist.
#' @param file_type A character string indicating the type of input file.
#'    Must be "csv" (for comma-separated) or "tsv" (for tab-separated).
#' @param time_step_interval A numeric value representing the interval between
#'    time steps for the heat diffusion simulation (e.g., 0.01).
#' @param stabilization_threshold A numeric value defining the threshold for
#'    determining stabilization. If the absolute difference in Spearman correlation
#'    between successive heat vectors falls below this threshold for a specified
#'    window, the process is considered stabilized.
#' @param stabilization_window_size An integer specifying the number of consecutive
#'    time steps over which the correlation difference must remain below the
#'    `stabilization_threshold` for stabilization to be declared.
#' @param filter_other_metabolite_edges A logical value. If `TRUE`, when a
#'    metabolite is used as a seed, all other edges connected to *other* metabolite
#'    nodes (i.e., not the current seed) are excluded from the network for that
#'    specific diffusion run. This isolates the diffusion to pathways and other
#'    non-metabolite nodes. If `FALSE`, the full network is used for each seed.
#' @return The function's primary output consists of multiple
#'    CSV files (heat scores and correlation data) and JPG plots (correlation plots),
#'    saved to the specified `output_directory`, one set for each metabolite seed node.
#' @references
#' Carlin DE, Demchak B, Pratt D, Sage E, Ideker T. Network propagation in the cytoscape cyberinfrastructure. PLoS computational biology. 2017;13(10):e1005598.
#' @export
node_prior <- function(
  multi_layered_network_file,
  output_directory,
  file_type = c("csv", "tsv"),
  time_step_interval = 0.01,
  stabilization_threshold = 0.0001, # This is 1e-04
  stabilization_window_size = 10,
  filter_other_metabolite_edges
) {
  file_type <- match.arg(file_type)
  message("Starting node prioritization using LHD algorithm.")

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

  # 1. Load the multi-layered network file (FULL network data)
  message("\n1. Loading multi-layered network from: ", multi_layered_network_file)
  if (!file.exists(multi_layered_network_file)) {
    stop("Input network file not found: '", multi_layered_network_file, "'. Cannot proceed.")
  }

  combined_data_full <- tryCatch( # Renamed to combined_data_full to denote original data
    read_input_file(multi_layered_network_file, file_type = file_type, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error reading multi-layered network file '", multi_layered_network_file, "': ", e$message, sep = ""))
    }
  )

  # Validate required columns in combined_data_full (now expecting snake_case)
  required_cols_network <- c("Feature1", "Feature2", "edge_score", "edge_type") # Changed Edge_Score and Edge_Type
  if (!all(required_cols_network %in% colnames(combined_data_full))) {
    stop(paste("Input network file '", multi_layered_network_file, "' must contain columns: ", paste(required_cols_network, collapse = ", "), ". Please ensure it's the output from Module 2, Step 4.", sep = ""))
  }

  # Ensure edge_score is numeric
  if (!is.numeric(combined_data_full$edge_score)) { # Changed Edge_Score
    stop("Column 'edge_score' in the input network file must be numeric.") # Changed Edge_Score
  }

  # 2. Identify all unique metabolite nodes from the FULL network (for looping and for filtering)
  all_metabolite_nodes <- unique(
    dplyr::pull(
      dplyr::filter(combined_data_full, edge_type == "Pathway-Metabolite"), # Changed Edge_Type
      Feature2
    )
  )
  all_metabolite_nodes <- as.character(all_metabolite_nodes)

  if (length(all_metabolite_nodes) == 0) {
    stop("No metabolite nodes found in the network (based on 'Pathway-Metabolite' edge_type). Cannot proceed diffusion.") # Changed Edge_Type
  }
  message("Identified ", length(all_metabolite_nodes), " unique metabolite nodes for seeding and potential filtering.")

  # Define the function for heat vector at time t
  H_vector_func <- function(t_val, laplacian_matrix, initial_heat_vector) {
    # Ensure dimensions match before multiplication
    if (nrow(laplacian_matrix) != length(initial_heat_vector)) {
      stop("Dimension mismatch: Laplacian matrix rows (", nrow(laplacian_matrix), ") and initial heat vector length (", length(initial_heat_vector), ") do not match.")
    }
    expm::expm(-laplacian_matrix * t_val) %*% initial_heat_vector
  }

  # Function to find stabilization time and return correlations
  find_stabilization_data <- function(current_L, current_H0, time_interval, threshold, window_size) {
    time_steps_local <- seq(0, 1, by = time_interval)

    heat_vectors_over_time <- sapply(time_steps_local, function(t) {
      H_vector_func(t, current_L, current_H0)
    })

    if (ncol(heat_vectors_over_time) < 2) {
      # If only one time step or no progression, return empty correlations and last time step
      return(list(stabilization_time = time_steps_local[length(time_steps_local)], correlations_df = data.frame(Time = numeric(), Correlation = numeric())))
    }

    correlations <- numeric(ncol(heat_vectors_over_time) - 1)

    for (i in seq_along(correlations)) {
      H_current <- heat_vectors_over_time[, i + 1]
      H_prev <- heat_vectors_over_time[, i]
      correlations[i] <- ifelse(sd(H_current) == 0 || sd(H_prev) == 0, 1, stats::cor(H_current, H_prev, method = "spearman", use = "complete.obs"))
    }

    correlation_df <- data.frame(Time = time_steps_local[-1], Correlation = correlations) # Time corresponds to the later time point of the pair

    stabilization_t_found <- time_steps_local[length(correlations) + 1] # Default to last time if no stabilization

    # Ensure window_size is not larger than available correlations
    if (window_size >= length(correlations)) {
      stabilization_t_found <- time_steps_local[length(correlations) + 1]
    } else {
      for (i in seq(length(correlations) - window_size + 1)) {
        diffs <- abs(diff(correlations[i:(i + window_size - 1)]))
        if (all(diffs < threshold)) {
          stabilization_t_found <- time_steps_local[i+1] # Time at the start of the stabilization *window* of the second vector in the pair.
          break
        }
      }
    }

    return(list(stabilization_time = stabilization_t_found, correlations_df = correlation_df))
  }

  # 3. Loop through each metabolite node to perform heat diffusion
  message("\nPerforming heat diffusion for each metabolite seed node...")
  for (seed_metabolite_id in all_metabolite_nodes) { # Loop through all identified metabolites as seeds
    message("  Processing seed metabolite: '", seed_metabolite_id, "'")

    # --- Determine network for current seed based on filter_other_metabolite_edges ---
    if (filter_other_metabolite_edges) {
      message("    Applying filtering: Excluding all nodes that are metabolites EXCEPT the current seed ('", seed_metabolite_id, "') and their connected edges.")

      # Identify all other metabolite nodes to exclude (all metabolites MINUS the current seed)
      other_metabolite_nodes_to_exclude <- setdiff(all_metabolite_nodes, seed_metabolite_id)

      # Filter edges: keep only edges where *neither* Feature1 nor Feature2 is one of the 'other metabolite' nodes
      current_combined_data <- dplyr::filter(
        combined_data_full,
        !(Feature1 %in% other_metabolite_nodes_to_exclude) &
          !(Feature2 %in% other_metabolite_nodes_to_exclude)
      )

      # Check if the seed node itself is still present after filtering
      if (!seed_metabolite_id %in% unique(c(current_combined_data$Feature1, current_combined_data$Feature2))) {
        warning("  Seed metabolite '", seed_metabolite_id, "' is not present in the graph after filtering. This implies it was only connected to other metabolites, which were removed. Skipping diffusion for this seed.")
        next
      }

    } else {
      # Use the full network data (original behavior)
      current_combined_data <- combined_data_full
      message("    Using the full network for diffusion.")
    }

    # --- Convert edge_score to absolute value ---
    current_combined_data$edge_score <- abs(current_combined_data$edge_score) # Changed Edge_Score
    message("    edge_scores converted to absolute values.") # Changed Edge_Scores

    # --- Create graph and Laplacian matrix for the CURRENT network configuration ---
    # This block is now inside the loop, as the network structure changes per seed if filtered

    # Check for valid edges/nodes to create a graph
    if (nrow(current_combined_data) == 0) {
      warning("  No edges found for graph construction for seed '", seed_metabolite_id, "' after filtering. Skipping diffusion.")
      next
    }

    # Get all unique nodes from the current filtered data to ensure the graph includes them
    all_nodes_in_current_data <- unique(c(current_combined_data$Feature1, current_combined_data$Feature2))

    g_current <- igraph::graph_from_data_frame(d = current_combined_data[, c("Feature1", "Feature2")], directed = FALSE,
                                               vertices = all_nodes_in_current_data) # Explicitly define vertices

    # Assign weights (handling potential mismatches as before)
    graph_edges_for_weighting <- dplyr::ungroup(
      dplyr::mutate(
        dplyr::rowwise(igraph::as_data_frame(g_current, what = "edges")),
        Node_A = min(from, to), Node_B = max(from, to)
      )
    )

    current_combined_data_sorted <- dplyr::ungroup(
      dplyr::distinct(
        dplyr::select(
          dplyr::mutate(
            dplyr::rowwise(current_combined_data),
            Node_A = min(Feature1, Feature2), Node_B = max(Feature1, Feature2)
          ),
          Node_A, Node_B, edge_score # Changed Edge_Score
        )
      )
    )

    # Ensure correct matching for weights
    # Create a unique identifier for edges in both dataframes for robust matching
    edge_id_graph <- paste0(graph_edges_for_weighting$Node_A, "_", graph_edges_for_weighting$Node_B)
    edge_id_data <- paste0(current_combined_data_sorted$Node_A, "_", current_combined_data_sorted$Node_B)

    matched_weights <- current_combined_data_sorted$edge_score[match(edge_id_graph, edge_id_data)] # Changed Edge_Score

    igraph::E(g_current)$weight <- matched_weights

    if (any(is.na(igraph::E(g_current)$weight))) {
      warning("  Some edges in the graph for seed '", seed_metabolite_id, "' could not be matched to an 'edge_score' in the input data. Assigning 0 weight to unmatched edges.") # Changed Edge_Score
      igraph::E(g_current)$weight[is.na(igraph::E(g_current)$weight)] <- 0
    }

    # If the graph has become too small (e.g., only 1 node, no edges, or disconnected)
    if (igraph::vcount(g_current) < 2 || igraph::ecount(g_current) == 0 || !igraph::is.connected(g_current)) {
      warning("  Graph for seed '", seed_metabolite_id, "' is too sparse, disconnected, or lacks sufficient nodes/edges after filtering (nodes: ", igraph::vcount(g_current), ", edges: ", igraph::ecount(g_current), ", connected: ", igraph::is.connected(g_current), "). Skipping diffusion for this seed as Laplacian cannot be computed meaningfully.")
      next
    }

    L_current <- igraph::laplacian_matrix(g_current, weights = igraph::E(g_current)$weight)
    L_current <- as.matrix(L_current)

    # Initialize the heat vector H_0 for the current seed within the CURRENT graph's nodes
    H_0_current_graph <- numeric(igraph::vcount(g_current))
    names(H_0_current_graph) <- igraph::V(g_current)$name
    seed_index_current_graph <- which(igraph::V(g_current)$name == seed_metabolite_id)

    if (length(seed_index_current_graph) == 0) {
      warning("  Seed metabolite '", seed_metabolite_id, "' not found in the *filtered* graph's node set. This should not happen if previous checks are correct. Skipping diffusion.")
      next
    }
    H_0_current_graph[seed_index_current_graph] <- 1.0

    # Find stabilization time and get correlations data
    stabilization_data_result <- find_stabilization_data(L_current, H_0_current_graph, time_step_interval, stabilization_threshold, stabilization_window_size)
    stabilization_t <- stabilization_data_result$stabilization_time
    correlation_df <- stabilization_data_result$correlations_df

    message("    Stabilization time for '", seed_metabolite_id, "': t = ", round(stabilization_t, 4))

    # Calculate final heat scores at stabilization time
    final_heat_scores <- H_vector_func(stabilization_t, L_current, H_0_current_graph)

    # Create output data frame for this metabolite
    output_df <- dplyr::arrange(
      data.frame(
        Node = igraph::V(g_current)$name, # Nodes from the current filtered graph
        Heat_score = round(final_heat_scores, 10), # Changed Heat_Score to Heat_score
        stringsAsFactors = FALSE
      ),
      dplyr::desc(Heat_score) # Sort by heat score (descending)
    )

    # Generate and save output files with cleaned input file name
    cleaned_seed_id <- gsub("[^A-Za-z0-9_]", "", seed_metabolite_id) # Clean seed ID for filename

    # Save heat scores to CSV
    output_file_name_heat <- paste0("heat_scores_", cleaned_seed_id, "_", cleaned_input_file_name, ".csv")
    output_path_heat <- file.path(output_directory, output_file_name_heat)
    write.csv(output_df, file = output_path_heat, row.names = FALSE)
    message("    Saved heat scores for '", seed_metabolite_id, "' to: ", output_path_heat)

    # Save correlation data to CSV
    output_file_name_correlation_csv <- paste0("spearman_correlations_", cleaned_seed_id, "_", cleaned_input_file_name, ".csv")
    output_path_correlation_csv <- file.path(output_directory, output_file_name_correlation_csv)
    write.csv(correlation_df, file = output_path_correlation_csv, row.names = FALSE)
    message("    Saved Spearman correlation data for '", seed_metabolite_id, "' to: ", output_path_correlation_csv)

    # Generate and save correlation plot
    correlation_plot <- ggplot2::ggplot(correlation_df, ggplot2::aes(x = Time, y = Correlation)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = stabilization_t, linetype = "dashed", color = "#a62140") +
      ggplot2::geom_text(
        ggplot2::aes(
          x = stabilization_t,
          y = max(Correlation, na.rm = TRUE) * 0.5, # Position text dynamically
          label = paste("t =", round(stabilization_t, 4))
        ),
        color = "#a62140", hjust = -0.1, vjust = 0.5, size = 5, fontface = "bold"
      ) +
      ggplot2::xlab("Time step") +
      ggplot2::ylab("Spearman correlation with previous time step") +
      ggplot2::scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 1)
      )

    output_file_name_plot <- paste0("correlation_plot_", cleaned_seed_id, "_", cleaned_input_file_name, ".jpg")
    output_path_plot <- file.path(output_directory, output_file_name_plot)
    ggplot2::ggsave(output_path_plot, plot = correlation_plot, width = 8, height = 5, dpi = 600)
    message("    Saved correlation plot for '", seed_metabolite_id, "' to: ", output_path_plot)
  }

  message("Node prioritization complete.")
  return(invisible(NULL))
}
