#' Multi-layered Network Construction (End-to-End)
#'
#' @details
#' Orchestrates the construction of a multi-layered network by standardizing inputs,
#' executing the layer-building functions (PPN, MPN, PMN), and triggering integration.
#'
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
  # --- UPDATED DEFAULTS TO MATCH YOUR THEME ---
  node_colors = c("Microbe" = "#9AA374", "Pathway" = "#C1ABAD", "Metabolite" = "#4E7286"),
  node_shapes = c("Microbe" = "hexagon", "Pathway" = "dot", "Metabolite" = "diamond"),
  base_node_size = 6, plot_width = 12, plot_height = 10, plot_dpi = 600
) {
  format <- match.arg(format)
  ppn_da_method <- match.arg(ppn_da_method)
  ppn_map_database <- match.arg(ppn_map_database)
  ppn_jaccard_method <- match.arg(ppn_jaccard_method)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  message("--- Preprocessing Data (Format: ", format, ") ---")

  processed_contrib_file <- file.path(output_dir, "processed_contribution.csv")
  taxonomy_file_to_pass <- NULL

  if (format == "universal") {
    df <- read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE)
    required <- c("SampleID", "PathwayID", "TaxonID", "contribution")
    if (!all(required %in% colnames(df))) stop("Universal format requires columns: ", paste(required, collapse = ", "))

    colnames(df)[colnames(df) == "contribution"] <- "taxon_function_abun"
    colnames(df)[colnames(df) == "PathwayID"] <- "FunctionID"
    colnames(df)[colnames(df) == "TaxonID"] <- "FeatureID"
    df$TaxonID <- df$FeatureID
    write.csv(df, processed_contrib_file, row.names = FALSE)
  } else if (format == "picrust") {
    if (is.null(taxonomy_file)) stop("Taxonomy file is required for PICRUSt format.")

    contrib <- read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE)
    tax <- read_input_file(taxonomy_file, file_type = "csv", stringsAsFactors = FALSE)

    merged <- merge(contrib, tax, by = "FeatureID", all = FALSE)
    keep_cols <- c("SampleID", "PathwayID", "FeatureID", "TaxonID", "taxon_function_abun")
    if (!all(keep_cols %in% colnames(merged))) stop("Missing required columns. Check exact spelling.")

    final_df <- merged[, keep_cols]
    colnames(final_df)[colnames(final_df) == "PathwayID"] <- "FunctionID"
    write.csv(final_df, processed_contrib_file, row.names = FALSE)
  } else if (format == "humann") {
    df <- read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE, check.names = FALSE)
    col1_name <- colnames(df)[1]
    sample_cols <- colnames(df)[-1]
    long_df <- tidyr::pivot_longer(df, cols = dplyr::all_of(sample_cols), names_to = "SampleID", values_to = "taxon_function_abun")
    long_df <- long_df[grepl("\\|", long_df[[col1_name]]), ]
    split_data <- stringr::str_split_fixed(long_df[[col1_name]], "\\|", 2)
    long_df$FunctionID <- split_data[, 1]
    long_df$FeatureID <- split_data[, 2]
    long_df$TaxonID <- long_df$FeatureID
    final_df <- long_df[, c("SampleID", "FunctionID", "FeatureID", "TaxonID", "taxon_function_abun")]
    write.csv(final_df, processed_contrib_file, row.names = FALSE)
  }

  ppn_dir <- file.path(output_dir, "ppn_output")
  ppn_results <- con_ppn_int(
    gene_abun_file = gene_abun_file, metadata_file = metadata_file, map_file = map_file,
    output_dir = ppn_dir, ppn_da_method = ppn_da_method, ppn_map_database = ppn_map_database,
    ppn_rank_by = ppn_rank_by, ppn_p_adjust_method = ppn_p_adjust_method,
    ppn_pvalue_cutoff = ppn_pvalue_cutoff, ppn_jaccard_cutoff = ppn_jaccard_cutoff,
    ppn_jaccard_method = ppn_jaccard_method, comparisons_list = comparisons_list
  )

  gsea_files <- ppn_results$gsea_paths
  jaccard_files <- ppn_results$jaccard_paths
  if (length(gsea_files) == 0) stop("No significant GSEA results found.")

  final_networks <- c()
  for (i in seq_along(gsea_files)) {
    curr_gsea <- gsea_files[i]
    curr_jaccard <- jaccard_files[i]

    message("\nProcessing integration for: ", basename(curr_gsea))

    mpn_dir <- file.path(output_dir, "mpn_output")
    mpn_files <- con_mpn_int(
      path_con_file = processed_contrib_file, metadata_file = metadata_file,
      taxonomy_file = taxonomy_file_to_pass, output_dir = mpn_dir, mpn_filtering = mpn_filtering
    )

    if (length(mpn_files) == 0) next
    curr_mpn <- mpn_files[1]

    pmn_dir <- file.path(output_dir, "pmn_output")
    pmn_files <- con_pmn_int(
      path_abun_file = path_abun_file, met_con_file = met_con_file, gsea_file = curr_gsea,
      metadata_file = metadata_file, output_dir = pmn_dir, pmn_corr_method = pmn_corr_method,
      pmn_filter_by = pmn_filter_by, pmn_corr_cutoff = pmn_corr_cutoff,
      pmn_pvalue_cutoff = pmn_pvalue_cutoff, pmn_q_value_cutoff = pmn_q_value_cutoff,
      pmn_p_adjust_method = pmn_p_adjust_method
    )
    curr_pmn <- if (length(pmn_files) > 0) pmn_files[1] else NULL

    mln_dir <- file.path(output_dir, "mln_final")
    final_net <- con_mln_int(
      gsea_file = curr_gsea, mpn_file = curr_mpn, ppn_file = curr_jaccard, pmn_file = curr_pmn,
      output_dir = mln_dir, visualize = visualize, layout_method = layout_method,
      node_colors = node_colors, node_shapes = node_shapes, base_node_size = base_node_size,
      plot_width = plot_width, plot_height = plot_height, plot_dpi = plot_dpi
    )
    final_networks <- c(final_networks, final_net)
  }

  message("Multi-Layered Network Pipeline Complete.")
  final_networks
}
