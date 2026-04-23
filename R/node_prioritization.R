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
#'    integrated multi-layered network data file (output
#'    from `con.mln` or `con.mln.all`). Expected columns: 'Feature1',
#'    'Feature2', 'edge_score' (numeric), and 'edge_type'.
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
#'    specific diffusion run. If `FALSE`, the full network is used for each seed.
#' @param visualize Logical. If TRUE, generates stabilization curve and lollipop rank plots. Defaults to TRUE.
#' @param top_n_plot Integer. Top N nodes to display in the lollipop rank chart.
#' @param node_colors Named character vector defining plotting colors for node types.
#' @param plot_width Numeric. Figure width in inches.
#' @param plot_height Numeric. Figure height in inches.
#' @param plot_dpi Numeric. Output image resolution (DPI).
#' @return The function's primary output consists of multiple
#'    CSV files (heat scores and correlation data) and JPG plots (correlation plots),
#'    saved to the specified `output_directory`, one set for each metabolite seed node.
#' @references
#' Carlin DE, Demchak B, Pratt D, Sage E, Ideker T. Network propagation in the cytoscape cyberinfrastructure. PLoS computational biology. 2017;13(10):e1005598.
#' @export
utils::globalVariables(c("edge_type", "Feature2", "Feature1", "from", "to", "Node_A", "Node_B", "edge_score", "Heat_score", "Time", "Correlation", "type"))

node_prior <- function(
  multi_layered_network_file,
  output_directory,
  file_type = c("csv", "tsv"),
  time_step_interval = 0.01,
  stabilization_threshold = 0.0001,
  stabilization_window_size = 10,
  filter_other_metabolite_edges,
  visualize = TRUE,
  top_n_plot = 20, 
  node_colors = c("Microbe" = "#D55E00", "Pathway" = "#0072B2", "Metabolite" = "#009E73"),
  plot_width = 10, 
  plot_height = 8, 
  plot_dpi = 600
) {
  file_type <- match.arg(file_type)
  message("Starting node prioritization using LHD algorithm.")

  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    message("Created output directory: ", output_directory)
  }

  input_file_base_name <- tools::file_path_sans_ext(basename(multi_layered_network_file))
  cleaned_input_file_name <- gsub("[^A-Za-z0-9_]", "", input_file_base_name)

  message("\n1. Loading multi-layered network from: ", multi_layered_network_file)
  if (!file.exists(multi_layered_network_file)) stop("Input network file not found.")

  combined_data_full <- read_input_file(multi_layered_network_file, file_type = file_type, stringsAsFactors = FALSE)

  required_cols_network <- c("Feature1", "Feature2", "edge_score", "edge_type")
  if (!all(required_cols_network %in% colnames(combined_data_full))) {
    stop("Input network file missing required columns.")
  }

  all_metabolite_nodes <- as.character(unique(dplyr::pull(dplyr::filter(combined_data_full, edge_type == "Pathway-Metabolite"), Feature2)))

  if (length(all_metabolite_nodes) == 0) stop("No metabolite nodes found in the network.")
  
  H_vector_func_fast <- function(t_val, U, lambda, U_T_H0) { as.vector(U %*% (exp(-lambda * t_val) * U_T_H0)) }

  message("\nPerforming heat diffusion for each metabolite seed node...")
  for (seed_metabolite_id in all_metabolite_nodes) {
    message("  Processing seed metabolite: '", seed_metabolite_id, "'")

    if (filter_other_metabolite_edges) {
      exclude <- setdiff(all_metabolite_nodes, seed_metabolite_id)
      current_combined_data <- dplyr::filter(combined_data_full, !(Feature1 %in% exclude) & !(Feature2 %in% exclude))
      if (!seed_metabolite_id %in% unique(c(current_combined_data$Feature1, current_combined_data$Feature2))) next
    } else {
      current_combined_data <- combined_data_full
    }

    current_combined_data$edge_score <- abs(current_combined_data$edge_score)
    all_nodes <- unique(c(current_combined_data$Feature1, current_combined_data$Feature2))
    
    g_current <- igraph::graph_from_data_frame(d = current_combined_data[, c("Feature1", "Feature2")], directed = FALSE, vertices = all_nodes)
    
    df_edges <- dplyr::mutate(dplyr::rowwise(current_combined_data), Node_A = min(Feature1, Feature2), Node_B = max(Feature1, Feature2))
    g_edges <- dplyr::mutate(dplyr::rowwise(igraph::as_data_frame(g_current, what = "edges")), Node_A = min(from, to), Node_B = max(from, to))
    igraph::E(g_current)$weight <- df_edges$edge_score[match(paste0(g_edges$Node_A, "_", g_edges$Node_B), paste0(df_edges$Node_A, "_", df_edges$Node_B))]
    igraph::E(g_current)$weight[is.na(igraph::E(g_current)$weight)] <- 0

    L_current <- as.matrix(igraph::laplacian_matrix(g_current, weights = igraph::E(g_current)$weight))
    eig <- eigen(L_current, symmetric = TRUE)
    
    H_0 <- numeric(igraph::vcount(g_current))
    names(H_0) <- igraph::V(g_current)$name
    H_0[which(names(H_0) == seed_metabolite_id)] <- 1.0
    U_T_H0 <- t(eig$vectors) %*% H_0

    t_seq <- seq(0, 1, by = time_step_interval)
    h_mat <- sapply(t_seq, function(t) H_vector_func_fast(t, eig$vectors, eig$values, U_T_H0))
    corrs <- sapply(1:(ncol(h_mat)-1), function(i) stats::cor(h_mat[, i+1], h_mat[, i], method = "spearman"))
    
    stab_t <- t_seq[length(t_seq)]
    for (i in seq_len(length(corrs) - stabilization_window_size + 1)) {
      if (all(abs(diff(corrs[i:(i + stabilization_window_size - 1)])) < stabilization_threshold)) {
        stab_t <- t_seq[i + 1]; break
      }
    }

    message("    Stabilization time for '", seed_metabolite_id, "': t = ", round(stab_t, 4))
    final_heat <- H_vector_func_fast(stab_t, eig$vectors, eig$values, U_T_H0)
    output_df <- dplyr::arrange(data.frame(Node = names(H_0), Heat_score = final_heat, stringsAsFactors = FALSE), dplyr::desc(Heat_score))
    
    cln_seed <- gsub("[^A-Za-z0-9_]", "", seed_metabolite_id)
    write.csv(output_df, file.path(output_directory, paste0("heat_scores_", cln_seed, ".csv")), row.names = FALSE)

    if (visualize) {
      library(ggplot2)
      
      corr_df <- data.frame(Time = t_seq[-1], Correlation = corrs)
      p1 <- ggplot2::ggplot(corr_df, ggplot2::aes(x = Time, y = Correlation)) +
        ggplot2::geom_line(linewidth = 1, color = "black") +
        ggplot2::geom_vline(xintercept = stab_t, linetype = "dashed", color = "#D55E00", linewidth = 1) +
        ggplot2::annotate("text", x = stab_t, y = min(corrs), label = paste("Stab t =", stab_t), color = "#D55E00", hjust = -0.1, family = "sans") +
        ggplot2::labs(x = "Time Step", y = "Spearman Correlation") +
        ggplot2::theme_classic(base_family = "sans") +
        ggplot2::theme(axis.title = ggplot2::element_text(face = "bold"), plot.margin = ggplot2::margin(10, 10, 10, 10))
      
      ggplot2::ggsave(file.path(output_directory, paste0("lhd_curve_", cln_seed, ".pdf")), plot = p1, width = 6, height = 4, device = cairo_pdf)

      top_df <- head(output_df, top_n_plot)
      top_df$type <- sapply(top_df$Node, function(x) {
        if (grepl("d__|p__|c__|o__|f__|g__|s__|Bacteria", x)) "Microbe" else if (grepl("^ko[0-9]+", x)) "Pathway" else "Metabolite"
      })
      top_df$type <- factor(top_df$type, levels = c("Microbe", "Pathway", "Metabolite"))

      p2 <- ggplot2::ggplot(top_df, ggplot2::aes(x = reorder(Node, Heat_score), y = Heat_score, color = type)) +
        ggplot2::geom_segment(ggplot2::aes(xend = Node, yend = 0), linewidth = 1.2) +
        ggplot2::geom_point(size = 4) +
        ggplot2::scale_color_manual(name = "Node Type", values = node_colors) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = "Network Node", y = "Diffusion Heat Score") +
        ggplot2::theme_classic(base_family = "sans") +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(size = 10, face = "bold", color = "black"),
          axis.title = ggplot2::element_text(face = "bold", size = 12),
          legend.position = "right",
          plot.margin = ggplot2::margin(10, 10, 10, 10)
        )

      ggplot2::ggsave(file.path(output_directory, paste0("prioritization_rank_", cln_seed, ".pdf")), plot = p2, width = plot_width, height = plot_height, device = cairo_pdf)
      ggplot2::ggsave(file.path(output_directory, paste0("prioritization_rank_", cln_seed, ".png")), plot = p2, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
    }
  }
  invisible(NULL)
}