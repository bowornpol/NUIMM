#' Hub identification using Maximal Clique Centrality (MCC) algorithm
#'
#' This function identifies hub nodes in a multi-layered network using the
#' Maximal Clique Centrality (MCC) algorithm and optionally visualizes the results.
#'
#' @details
#' [Standard details...]
#'
#' @param multi_layered_network_file Path to input CSV/TSV.
#' @param output_directory Path to save results.
#' @param file_type "csv" or "tsv".
#' @param top_n_hubs Integer. Number of top hubs to filter and visualize.
#' @param visualize Logical. If TRUE, saves a PDF plot of the top hubs and their connections. Defaults to FALSE.
#' @export
iden_hub <- function(
  multi_layered_network_file,
  output_directory,
  file_type = c("csv", "tsv"),
  top_n_hubs = NULL,
  visualize = TRUE
) {
  file_type <- match.arg(file_type)
  message("Starting hub identification using Maximal Clique Centrality (MCC) algorithm.")

  # Create output directory if it doesn't exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    message("Created output directory: ", output_directory)
  }

  # Extract base name for filenames
  input_file_base_name <- tools::file_path_sans_ext(basename(multi_layered_network_file))
  cleaned_input_file_name <- gsub("[^A-Za-z0-9_]", "", input_file_base_name)

  # 1. Load the network
  message("\n1. Loading network from: ", multi_layered_network_file)
  if (!file.exists(multi_layered_network_file)) stop("File not found.")

  # Standardized reading based on file_type (Replacing custom read_input_file for portability)
  if (file_type == "csv") {
    network_data <- read.csv(multi_layered_network_file, stringsAsFactors = FALSE)
  } else {
    network_data <- read.delim(multi_layered_network_file, stringsAsFactors = FALSE)
  }

  required_cols <- c("Feature1", "Feature2")
  if (!all(required_cols %in% colnames(network_data))) {
    stop("Network file must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # 2. Create graph
  g <- igraph::graph_from_data_frame(d = network_data[, required_cols], directed = FALSE)
  message("  Graph created with ", igraph::vcount(g), " nodes.")

  # 3. Find maximal cliques
  message("\n2. Finding maximal cliques...")
  cliques <- tryCatch(
    igraph::max_cliques(g),
    error = function(e) stop("Error finding cliques: ", e$message)
  )

  # 4. Calculate MCC scores
  message("\n3. Calculating MCC scores...")
  mcc_scores <- setNames(numeric(igraph::vcount(g)), igraph::V(g)$name)

  for (clique_nodes_indices in cliques) {
    clique_size <- length(clique_nodes_indices)
    clique_score <- factorial(clique_size - 1)
    for (node_index in clique_nodes_indices) {
      node_name <- igraph::V(g)$name[node_index]
      mcc_scores[node_name] <- mcc_scores[node_name] + clique_score
    }
  }

  # Handle special case: No edges between neighbors
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

  # 5. Rank Results
  message("\n4. Ranking nodes.")
  hub_results_df <- dplyr::arrange(
    data.frame(Node = names(mcc_scores), MCC_score = mcc_scores, stringsAsFactors = FALSE),
    dplyr::desc(MCC_score)
  )

  # 6. Filter Top N
  if (!is.null(top_n_hubs) && top_n_hubs > 0) {
    if (top_n_hubs > nrow(hub_results_df)) {
      warning("Requested top_n_hubs > total nodes. Returning all.")
    } else {
      hub_results_df <- head(hub_results_df, n = top_n_hubs)
      message("  Filtered to top ", top_n_hubs, " hub nodes.")
    }
  }

  # --- NEW: VISUALIZATION BLOCK ---
  if (visualize) {
    message("\n  Generating visualization...")

    # 1. Identify nodes to plot
    nodes_to_plot <- hub_results_df$Node

    # 2. Create a subgraph containing ONLY these nodes and connections between them
    sub_g <- igraph::induced_subgraph(g, vids = nodes_to_plot)

    # 3. Scale node sizes based on MCC Score (Min size 5, Max size 20)
    scores <- hub_results_df$MCC_score
    if (max(scores) == min(scores)) {
      sizes <- rep(15, length(scores))
    } else {
      # Normalize between 0 and 1, then scale
      sizes <- 5 + 15 * ((scores - min(scores)) / (max(scores) - min(scores)))
    }
    # Match sizes to the graph node order
    graph_order_names <- igraph::V(sub_g)$name
    matched_sizes <- sizes[match(graph_order_names, hub_results_df$Node)]

    # 4. Save to PDF
    plot_filename <- paste0("hub_plot_", cleaned_input_file_name,
                            if(!is.null(top_n_hubs)) paste0("_top", top_n_hubs) else "", ".pdf")
    plot_path <- file.path(output_directory, plot_filename)

    pdf(plot_path, width = 10, height = 10)
    plot(sub_g,
         vertex.label = igraph::V(sub_g)$name,
         vertex.label.cex = 0.8,       # Font size
         vertex.label.color = "black",
         vertex.size = matched_sizes,  # Dynamic size
         vertex.color = "tomato",      # Node color
         vertex.frame.color = "darkred",
         edge.color = "gray70",
         edge.width = 1.5,
         layout = igraph::layout_nicely(sub_g),
         main = paste("Top", length(nodes_to_plot), "Hubs (MCC Algorithm)"))
    dev.off()
    message("  Saved visualization to: ", plot_path)
  }
  # --------------------------------

  # 7. Save CSV
  output_filename <- paste0("hub_mcc_", cleaned_input_file_name,
                            if(!is.null(top_n_hubs)) paste0("_top", top_n_hubs) else "", ".csv")
  output_filepath <- file.path(output_directory, output_filename)
  message("\n5. Saving CSV results to: ", output_filepath)
  write.csv(hub_results_df, output_filepath, row.names = FALSE)

  message("Hub identification complete.")
  return(invisible(NULL))
}
