#' Multi-layered Network Construction (End-to-End)
#'
#' @details
#' Constructs an integrated multi-layered network by sequentially building
#' three sub-network layers (Microbe-Pathway, Pathway-Pathway, Pathway-Metabolite),
#' then merging them into a single unified network with an interactive HTML visualization.
#'
#' @param gene_abun_file Path to gene/KO abundance table (samples as columns).
#' @param path_abun_file Path to pathway abundance table (samples as columns).
#' @param path_con_file Path to pathway contribution file linking taxa to functions.
#' @param met_con_file Path to metabolite concentration table (samples as rows).
#' @param metadata_file Path to sample metadata CSV with 'SampleID' and 'class' columns.
#' @param taxonomy_file Path to taxonomy mapping file. Required for PICRUSt format.
#' @param map_file Path to pathway-to-gene mapping file for GSEA.
#' @param output_dir Path to output directory.
#' @param format Input data format: "universal", "humann", or "picrust".
#' @param ppn_da_method Differential abundance method: "deseq2", "edger", "maaslin2", or "simple".
#' @param ppn_map_database Pathway database: "kegg", "metacyc", or "custom".
#' @param ppn_rank_by Gene ranking metric for GSEA: "signed_log_pvalue", "log2foldchange", or "pvalue".
#' @param ppn_p_adjust_method P-value adjustment method for GSEA.
#' @param ppn_pvalue_cutoff P-value cutoff for GSEA significance.
#' @param ppn_jaccard_cutoff Minimum Jaccard index to retain pathway-pathway edges.
#' @param ppn_jaccard_method Jaccard calculation source: "gsea_core" or "map_file".
#' @param comparisons_list Optional list of pairwise group comparisons.
#' @param mpn_filtering Microbe-pathway filtering: "unfiltered", "mean", "median", or "topN\%".
#' @param pmn_corr_method Correlation method for pathway-metabolite: "spearman", "pearson", or "kendall".
#' @param pmn_filter_by Significance filter: "none", "p_value", or "q_value".
#' @param pmn_corr_cutoff Minimum absolute correlation for pathway-metabolite edges.
#' @param pmn_pvalue_cutoff P-value cutoff for pathway-metabolite edges.
#' @param pmn_q_value_cutoff Q-value cutoff for pathway-metabolite edges.
#' @param pmn_p_adjust_method P-value adjustment method for pathway-metabolite correlations.
#' @param visualize Logical. If TRUE, generates interactive HTML visualization.
#' @param layout_method Network layout algorithm.
#' @param node_colors Named character vector of colors for each node group.
#' @param node_shapes Named character vector of shapes for each node group.
#' @param base_node_size Base size for network nodes.
#' @param plot_width Figure width in inches.
#' @param plot_height Figure height in inches.
#' @param plot_dpi Output image resolution (DPI).
#' @return Character vector of output file paths.
#' @export
con_mln <- function(
  gene_abun_file, path_abun_file, path_con_file, met_con_file, metadata_file,
  taxonomy_file = NULL, map_file, output_dir,
  format = c("universal", "humann", "picrust"),
  ppn_da_method = c("deseq2", "edger", "maaslin2", "simple"),
  ppn_map_database = c("kegg", "metacyc", "custom"),
  ppn_rank_by = c("signed_log_pvalue", "log2foldchange", "pvalue"),
  ppn_p_adjust_method = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none"),
  ppn_pvalue_cutoff = 0.05, ppn_jaccard_cutoff = 0.2,
  ppn_jaccard_method = c("gsea_core", "map_file"), comparisons_list = NULL,
  mpn_filtering = c("unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%"),
  pmn_corr_method = c("spearman", "pearson", "kendall"),
  pmn_filter_by = c("none", "p_value", "q_value"),
  pmn_corr_cutoff = 0.3, pmn_pvalue_cutoff = 0.05, pmn_q_value_cutoff = 0.05,
  pmn_p_adjust_method = "fdr", visualize = TRUE, layout_method = "sugiyama",
  node_colors = c("Microbe" = "#9AA374", "Pathway" = "#C1ABAD", "Metabolite" = "#4E7286"),
  node_shapes = c("Microbe" = "hexagon", "Pathway" = "dot", "Metabolite" = "diamond"),
  base_node_size = 6, plot_width = 12, plot_height = 10, plot_dpi = 600
) {
  format <- match.arg(format)
  ppn_da_method <- match.arg(ppn_da_method)
  ppn_map_database <- match.arg(ppn_map_database)
  ppn_jaccard_method <- match.arg(ppn_jaccard_method)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  message("--- Preparing Data Processing Pipeline (Format: ", format, ") ---")
  processed_contrib_file <- file.path(output_dir, "processed_contribution.csv")

  # Standardized Preprocessing
  if (format == "picrust") {
    if (is.null(taxonomy_file)) stop("Taxonomy file is required for PICRUSt format.")
    contrib <- read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE)
    tax <- read_input_file(taxonomy_file, file_type = "csv", stringsAsFactors = FALSE)
    merged <- merge(contrib, tax, by = "FeatureID", all = FALSE)
    keep_cols <- c("SampleID", "PathwayID", "FeatureID", "TaxonID", "taxon_function_abun")
    if (!all(keep_cols %in% colnames(merged))) stop("Missing required columns after taxonomy merge.")
    final_df <- merged[, keep_cols]
    colnames(final_df)[colnames(final_df) == "PathwayID"] <- "FunctionID"
    write.csv(final_df, processed_contrib_file, row.names = FALSE)
  } else if (format == "universal") {
    df <- read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE)
    colnames(df)[colnames(df) == "contribution"] <- "taxon_function_abun"
    colnames(df)[colnames(df) == "PathwayID"] <- "FunctionID"
    colnames(df)[colnames(df) == "TaxonID"] <- "FeatureID"
    df$TaxonID <- df$FeatureID
    write.csv(df, processed_contrib_file, row.names = FALSE)
  } else if (format == "humann") {
    df <- read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE, check.names = FALSE)
    long_df <- tidyr::pivot_longer(df, cols = -1, names_to = "SampleID", values_to = "taxon_function_abun")
    split_data <- stringr::str_split_fixed(long_df[[1]], "\\|", 2)
    long_df$FunctionID <- split_data[, 1]
    long_df$FeatureID <- split_data[, 2]
    long_df$TaxonID <- long_df$FeatureID
    write.csv(long_df[, c("SampleID", "FunctionID", "FeatureID", "TaxonID", "taxon_function_abun")], processed_contrib_file, row.names = FALSE)
  }

  # Execute Layers safely
  ppn_results <- con_ppn_int(
    gene_abun_file = gene_abun_file, metadata_file = metadata_file, map_file = map_file,
    output_dir = file.path(output_dir, "ppn_output"), ppn_da_method = ppn_da_method,
    ppn_map_database = ppn_map_database, ppn_rank_by = ppn_rank_by,
    ppn_p_adjust_method = ppn_p_adjust_method, ppn_pvalue_cutoff = ppn_pvalue_cutoff,
    ppn_jaccard_cutoff = ppn_jaccard_cutoff, ppn_jaccard_method = ppn_jaccard_method,
    comparisons_list = comparisons_list
  )

  final_outputs <- c()
  for (i in seq_along(ppn_results$gsea_paths)) {
    curr_gsea <- ppn_results$gsea_paths[i]
    curr_jaccard <- ppn_results$jaccard_paths[i]
    message("\nInitiating multi-layered integration for: ", basename(curr_gsea))

    curr_mpn <- con_mpn_int(processed_contrib_file, metadata_file, NULL, file.path(output_dir, "mpn_output"), mpn_filtering)[1]
    curr_pmn <- con_pmn_int(path_abun_file, met_con_file, curr_gsea, metadata_file, file.path(output_dir, "pmn_output"), pmn_corr_method, pmn_filter_by, pmn_corr_cutoff, pmn_pvalue_cutoff, pmn_q_value_cutoff, pmn_p_adjust_method)[1]

    res_path <- con_mln_int(
      gsea_file = curr_gsea, mpn_file = curr_mpn, ppn_file = curr_jaccard, pmn_file = curr_pmn,
      output_dir = file.path(output_dir, "mln_final"), visualize = visualize,
      layout_method = layout_method, node_colors = node_colors, node_shapes = node_shapes,
      base_node_size = base_node_size, plot_width = plot_width, plot_height = plot_height, plot_dpi = plot_dpi
    )
    if (!is.null(res_path)) final_outputs <- c(final_outputs, res_path)
  }

  message("\n--- Multi-Layered Network Pipeline Successfully Completed ---")
  return(final_outputs)
}
