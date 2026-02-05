#' Internal MLN Integration & Visualization Helper
#'
#' Merges MPN, PPN, and PMN layers and visualizes.
#'
#' @param gsea_file Character path to GSEA file.
#' @param mpn_file Character path to MPN file.
#' @param ppn_file Character path to PPN file.
#' @param pmn_file Character path to PMN file.
#' @param output_dir Character path to output directory.
#' @return Final network file path.
#' @keywords internal
con_mln_int <- function(
    gsea_file,
    mpn_file,
    ppn_file,
    pmn_file,
    output_dir
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  mpn <- read.csv(mpn_file, stringsAsFactors = FALSE)
  ppn <- if(!is.na(ppn_file) && file.exists(ppn_file)) read.csv(ppn_file) else NULL
  pmn <- if(file.exists(pmn_file)) read.csv(pmn_file, stringsAsFactors = FALSE) else NULL
  gsea <- read.csv(gsea_file)

  valid_paths <- gsea$ID
  edges <- data.frame()

  mpn <- mpn[mpn$FunctionID %in% valid_paths, ]
  if (nrow(mpn) > 0) {
    edges <- rbind(edges, data.frame(from = mpn$TaxonID, to = mpn$FunctionID, weight = mpn$relative_contribution, type = "Microbe-Pathway"))
  }

  if (!is.null(ppn)) {
    ppn <- ppn[ppn$FunctionID_1 %in% valid_paths & ppn$FunctionID_2 %in% valid_paths, ]
    if (nrow(ppn) > 0) {
      edges <- rbind(edges, data.frame(from = ppn$FunctionID_1, to = ppn$FunctionID_2, weight = ppn$jaccard_index, type = "Pathway-Pathway"))
    }
  }

  if (!is.null(pmn)) {
    pmn <- pmn[pmn$FunctionID %in% valid_paths, ]
    if (nrow(pmn) > 0) {
      edges <- rbind(edges, data.frame(from = pmn$FunctionID, to = pmn$MetaboliteID, weight = pmn$correlation, type = "Pathway-Metabolite"))
    }
  }

  out_path <- file.path(output_dir, paste0("final_mln_", basename(gsea_file)))
  write.csv(edges, out_path, row.names = FALSE)

  if (nrow(edges) > 0) {
    tryCatch({
      library(igraph)
      library(ggraph)
      library(ggplot2)
      g <- graph_from_data_frame(edges, directed = FALSE)
      V(g)$type <- ifelse(V(g)$name %in% edges$from[edges$type == "Microbe-Pathway"], "Microbe", ifelse(V(g)$name %in% edges$to[edges$type == "Pathway-Metabolite"], "Metabolite", "Pathway"))

      p <- ggraph(g, layout = "fr") +
        geom_edge_link(aes(color = type, width = weight), alpha = 0.6) +
        geom_node_point(aes(color = type), size = 4) +
        geom_node_text(aes(label = name), repel = TRUE, size = 3) +
        theme_void() +
        labs(title = "Multi-Layered Network", color = "Node Type", edge_color = "Edge Type")

      ggsave(filename = gsub(".csv", ".png", out_path), plot = p, width = 12, height = 10)
    }, error = function(e) message("Visualization skipped: ", e$message))
  }
  return(out_path)
}
