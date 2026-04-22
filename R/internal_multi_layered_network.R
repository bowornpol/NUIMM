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
#' @return Final network file path.
#' @keywords internal
utils::globalVariables(c("type", "weight", "name"))

con_mln_int <- function(
  gsea_file,
  mpn_file,
  ppn_file,
  pmn_file,
  output_dir
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  mpn <- read_input_file(mpn_file, file_type = "csv", stringsAsFactors = FALSE)
  ppn <- if (!is.na(ppn_file) && file.exists(ppn_file)) read_input_file(ppn_file, file_type = "csv") else NULL
  pmn <- if (file.exists(pmn_file)) read_input_file(pmn_file, file_type = "csv", stringsAsFactors = FALSE) else NULL
  gsea <- read_input_file(gsea_file, file_type = "csv")

  valid_paths <- gsea$ID
  edges <- data.frame()

  mpn <- mpn[mpn$FunctionID %in% valid_paths, ]
  if (nrow(mpn) > 0) {
    edges <- rbind(edges, data.frame(Feature1 = mpn$TaxonID, Feature2 = mpn$FunctionID, edge_score = mpn$relative_contribution, edge_type = "Microbe-Pathway"))
  }

  if (!is.null(ppn)) {
    ppn <- ppn[ppn$FunctionID_1 %in% valid_paths & ppn$FunctionID_2 %in% valid_paths, ]
    if (nrow(ppn) > 0) {
      edges <- rbind(edges, data.frame(Feature1 = ppn$FunctionID_1, Feature2 = ppn$FunctionID_2, edge_score = ppn$jaccard_index, edge_type = "Pathway-Pathway"))
    }
  }

  if (!is.null(pmn)) {
    pmn <- pmn[pmn$FunctionID %in% valid_paths, ]
    if (nrow(pmn) > 0) {
      edges <- rbind(edges, data.frame(Feature1 = pmn$FunctionID, Feature2 = pmn$MetaboliteID, edge_score = pmn$correlation, edge_type = "Pathway-Metabolite"))
    }
  }

  out_path <- file.path(output_dir, paste0("final_mln_", basename(gsea_file)))
  write.csv(edges, out_path, row.names = FALSE)

  if (nrow(edges) > 0) {
    tryCatch(
      {
        library(igraph)
        library(ggraph)
        library(ggplot2)
        g <- igraph::graph_from_data_frame(edges, directed = FALSE)
        igraph::V(g)$type <- ifelse(igraph::V(g)$name %in% edges$Feature1[edges$edge_type == "Microbe-Pathway"], "Microbe",
          ifelse(igraph::V(g)$name %in% edges$Feature2[edges$edge_type == "Pathway-Metabolite"], "Metabolite", "Pathway")
        )

        # Assign hierarchical layers for structural grouping
        igraph::V(g)$layer <- ifelse(igraph::V(g)$type == "Microbe", 1,
          ifelse(igraph::V(g)$type == "Pathway", 2, 3)
        )

        # Use Sugiyama layout to stack layers hierarchically (Microbes top, Pathways middle, Metabolites bottom)
        lay <- igraph::layout_with_sugiyama(g, layers = igraph::V(g)$layer)$layout

        p <- ggraph::ggraph(g, layout = lay) +
          ggraph::geom_edge_link(ggplot2::aes(width = edge_score, alpha = edge_score), color = "gray60", show.legend = c(alpha = FALSE)) +
          ggraph::scale_edge_width(range = c(0.5, 2.0)) +
          ggraph::geom_node_point(ggplot2::aes(color = type), size = 6) +
          ggplot2::scale_color_manual(
            name = "Node Type",
            values = c("Microbe" = "#E64B35", "Pathway" = "#4DBBD5", "Metabolite" = "#00A087")
          ) +
          ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 3.5, fontface = "bold", family = "sans") +
          ggplot2::theme_void(base_family = "sans") +
          ggplot2::theme(
            legend.position = "right",
            legend.title = ggplot2::element_text(face = "bold", size = 12, color = "black"),
            legend.text = ggplot2::element_text(size = 10, color = "black"),
            plot.title = ggplot2::element_blank()
          ) +
          ggplot2::labs(edge_color = "Edge Type")

        ggplot2::ggsave(filename = gsub(".csv", ".png", out_path), plot = p, width = 14, height = 10, dpi = 600, bg = "white")
        ggplot2::ggsave(filename = gsub(".csv", ".pdf", out_path), plot = p, width = 14, height = 10)
      },
      error = function(e) warning("Visualization skipped: ", e$message)
    )
  }
  out_path
}
