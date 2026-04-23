#' Hub identification using Maximal Clique Centrality (MCC) algorithm
#'
#' @details
#' This function identifies influential "hub" nodes in a multi-layered network
#' using the Maximal Clique Centrality (MCC) algorithm. The results are exported
#' as a CSV file, with optional generation of network plot to visualize the top-ranked hubs.
#'
#' @param multi_layered_network_file Path to input CSV/TSV.
#' @param output_directory Path to save results.
#' @param file_type "csv" or "tsv".
#' @param top_n_hubs Integer. Number of top hubs to filter and visualize.
#' @param visualize Logical. If TRUE, saves PDF and PNG plots of the top hubs. Defaults to FALSE.
#' @param layout_method Character. Graph layout algorithm (e.g., "fr" for Fruchterman-Reingold).
#' @param color_palette Character. Viridis colorblind-safe palette mapping ("plasma", "viridis", etc).
#' @param base_node_size Numeric Vector. Range limits for scaled nodes c(min, max).
#' @param plot_width Numeric. Width of output plot in inches.
#' @param plot_height Numeric. Height of output plot in inches.
#' @param plot_dpi Numeric. DPI for PNG rendering.
#' @export
utils::globalVariables(c("MCC_score", "mcc_score", "name"))

iden_hub <- function(
  multi_layered_network_file,
  output_directory,
  file_type = c("csv", "tsv"),
  top_n_hubs = 20,
  visualize = TRUE,
  layout_method = "fr", 
  color_palette = "plasma",
  base_node_size = c(4, 12),
  plot_width = 10, 
  plot_height = 8, 
  plot_dpi = 600
) {
  file_type <- match.arg(file_type)
  message("Starting hub identification using Maximal Clique Centrality (MCC) algorithm.")

  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    message("Created output directory: ", output_directory)
  }

  input_file_base_name <- tools::file_path_sans_ext(basename(multi_layered_network_file))
  cleaned_input_file_name <- gsub("[^A-Za-z0-9_]", "", input_file_base_name)

  message("\n1. Loading network from: ", multi_layered_network_file)
  if (!file.exists(multi_layered_network_file)) stop("File not found.")

  network_data <- read_input_file(multi_layered_network_file, file_type = file_type, stringsAsFactors = FALSE)

  required_cols <- c("Feature1", "Feature2")
  if (!all(required_cols %in% colnames(network_data))) {
    stop("Network file must contain columns: ", paste(required_cols, collapse = ", "))
  }

  g <- igraph::graph_from_data_frame(d = network_data[, required_cols], directed = FALSE)
  message("  Graph created with ", igraph::vcount(g), " nodes.")

  message("\n2. Finding maximal cliques...")
  cliques <- tryCatch(
    igraph::max_cliques(g),
    error = function(e) stop("Error finding cliques: ", e$message)
  )

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

  message("\n4. Ranking nodes...")
  hub_results_df <- dplyr::arrange(
    data.frame(Node = names(mcc_scores), MCC_score = mcc_scores, stringsAsFactors = FALSE),
    dplyr::desc(MCC_score)
  )

  if (!is.null(top_n_hubs) && top_n_hubs > 0) {
    if (top_n_hubs > nrow(hub_results_df)) {
      warning("Requested top_n_hubs > total nodes. Returning all.")
    } else {
      hub_results_df <- head(hub_results_df, n = top_n_hubs)
      message("  Filtered to top ", top_n_hubs, " hub nodes.")
    }
  }

  if (visualize) {
    message("\n5. Generating network visualization...")
    if (!requireNamespace("ggraph", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
      stop("To use visualize=TRUE, you must install 'ggraph' and 'ggplot2'.")
    }

    set.seed(42)
    sub_g <- igraph::induced_subgraph(g, vids = hub_results_df$Node)
    matched_scores <- hub_results_df$MCC_score[match(igraph::V(sub_g)$name, hub_results_df$Node)]
    igraph::V(sub_g)$mcc_score <- matched_scores

    p <- ggraph::ggraph(sub_g, layout = layout_method) +
      ggraph::geom_edge_link(alpha = 0.4, color = "gray70") +
      ggraph::geom_node_point(ggplot2::aes(color = mcc_score, size = mcc_score)) +
      ggplot2::scale_size_continuous(range = base_node_size, guide = "none") +
      ggplot2::scale_color_viridis_c(option = color_palette, name = "MCC Score") +
      ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 3.5, fontface = "bold", family = "sans") +
      ggplot2::theme_void(base_family = "sans") +
      ggplot2::theme(
        legend.position = "right",
        legend.title = ggplot2::element_text(face = "bold", size = 12),
        legend.text = ggplot2::element_text(size = 10),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )

    base_plot_name <- paste0("hub_plot_", cleaned_input_file_name, if (!is.null(top_n_hubs)) paste0("_top", top_n_hubs) else "")

    ggplot2::ggsave(file.path(output_directory, paste0(base_plot_name, ".pdf")), plot = p, width = plot_width, height = plot_height, device = cairo_pdf)
    message("  Saved PDF visualization to: ", file.path(output_directory, paste0(base_plot_name, ".pdf")))

    ggplot2::ggsave(file.path(output_directory, paste0(base_plot_name, ".png")), plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
    message("  Saved PNG visualization to: ", file.path(output_directory, paste0(base_plot_name, ".png")))
  }

  output_filename <- paste0("hub_mcc_", cleaned_input_file_name, if (!is.null(top_n_hubs)) paste0("_top", top_n_hubs) else "", ".csv")
  output_filepath <- file.path(output_directory, output_filename)
  message("\n6. Saving CSV results to: ", output_filepath)
  write.csv(hub_results_df, output_filepath, row.names = FALSE)

  message("Hub identification complete.")
  invisible(NULL)
}