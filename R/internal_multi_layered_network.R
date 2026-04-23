#' Internal MLN Integration & Visualization Helper
#'
#' @details
#' Assembles the final network by merging edge lists into an `igraph` object and 
#' plots the hierarchical structure (Microbes, Pathways, Metabolites) using `ggraph`.
#'
#' @param gsea_file Character path to GSEA file.
#' @param mpn_file Character path to MPN file.
#' @param ppn_file Character path to PPN file.
#' @param pmn_file Character path to PMN file.
#' @param output_dir Character path to output directory.
#' @param visualize Logical internal flag to trigger plot rendering.
#' @param layout_method Internal layout algorithm parameter.
#' @param node_colors Internal color mappings.
#' @param node_shapes Internal shape mappings.
#' @param base_node_size Internal node size parameter.
#' @param plot_width Internal plot width.
#' @param plot_height Internal plot height.
#' @param plot_dpi Internal plot dpi.
#' @return Final network file path.
#' @keywords internal
utils::globalVariables(c("type", "weight", "name", "edge_score"))

con_mln_int <- function(
  gsea_file, mpn_file, ppn_file, pmn_file, output_dir,
  visualize, layout_method, node_colors, node_shapes, base_node_size, plot_width, plot_height, plot_dpi
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  mpn <- read_input_file(mpn_file, file_type = "csv", stringsAsFactors = FALSE)
  ppn <- if (!is.na(ppn_file) && file.exists(ppn_file)) read_input_file(ppn_file, file_type = "csv") else NULL
  pmn <- if (!is.null(pmn_file) && file.exists(pmn_file)) read_input_file(pmn_file, file_type = "csv", stringsAsFactors = FALSE) else NULL
  gsea <- read_input_file(gsea_file, file_type = "csv")

  valid_paths <- gsea$ID
  edges <- data.frame()

  mpn <- mpn[mpn$FunctionID %in% valid_paths, ]
  if (nrow(mpn) > 0) edges <- rbind(edges, data.frame(Feature1 = mpn$TaxonID, Feature2 = mpn$FunctionID, edge_score = mpn$relative_contribution, edge_type = "Microbe-Pathway"))

  if (!is.null(ppn)) {
    ppn <- ppn[ppn$FunctionID_1 %in% valid_paths & ppn$FunctionID_2 %in% valid_paths, ]
    if (nrow(ppn) > 0) edges <- rbind(edges, data.frame(Feature1 = ppn$FunctionID_1, Feature2 = ppn$FunctionID_2, edge_score = ppn$jaccard_index, edge_type = "Pathway-Pathway"))
  }

  if (!is.null(pmn)) {
    pmn <- pmn[pmn$FunctionID %in% valid_paths, ]
    if (nrow(pmn) > 0) edges <- rbind(edges, data.frame(Feature1 = pmn$FunctionID, Feature2 = pmn$MetaboliteID, edge_score = abs(pmn$correlation), edge_type = "Pathway-Metabolite"))
  }

  out_path <- file.path(output_dir, paste0("final_mln_", basename(gsea_file)))
  write.csv(edges, out_path, row.names = FALSE)

  if (visualize && nrow(edges) > 0) {
    tryCatch({
      library(igraph)
      library(ggraph)
      library(ggplot2)
      set.seed(42) # Deterministic layout for reproducibility

      g <- igraph::graph_from_data_frame(edges, directed = FALSE)
      igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% edges$Feature1[edges$edge_type == "Microbe-Pathway"], "Microbe",
                           ifelse(igraph::V(g)$name %in% edges$Feature2[edges$edge_type == "Pathway-Metabolite"], "Metabolite", "Pathway"))
      
      # Assign hierarchical layers for structural grouping
      igraph::V(g)$layer <- ifelse(igraph::V(g)$type == "Microbe", 1, ifelse(igraph::V(g)$type == "Pathway", 2, 3))
      igraph::V(g)$type <- factor(igraph::V(g)$type, levels = c("Microbe", "Pathway", "Metabolite"))

      lay <- if (layout_method == "sugiyama") igraph::layout_with_sugiyama(g, layers = igraph::V(g)$layer)$layout else layout_method

      p <- ggraph::ggraph(g, layout = lay) +
        ggraph::geom_edge_link(ggplot2::aes(width = edge_score, alpha = edge_score), color = "gray50", show.legend = c(alpha = FALSE)) +
        ggraph::scale_edge_width(range = c(0.2, 1.5), name = "Edge Score") +
        ggraph::geom_node_point(ggplot2::aes(fill = type, shape = type), size = base_node_size, color = "black", stroke = 0.5) +
        ggplot2::scale_fill_manual(name = "Node Type", values = node_colors) +
        ggplot2::scale_shape_manual(name = "Node Type", values = node_shapes) +
        ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 3, family = "sans") +
        ggplot2::theme_void(base_family = "sans") +
        ggplot2::theme(
          legend.position = "right",
          legend.title = ggplot2::element_text(face = "bold", size = 12),
          legend.text = ggplot2::element_text(size = 10),
          plot.margin = ggplot2::margin(10, 10, 10, 10)
        )

      ggplot2::ggsave(filename = gsub(".csv", ".pdf", out_path), plot = p, width = plot_width, height = plot_height, device = cairo_pdf)
      ggplot2::ggsave(filename = gsub(".csv", ".png", out_path), plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
    }, error = function(e) warning("Visualization skipped: ", e$message))
  }
  out_path
}