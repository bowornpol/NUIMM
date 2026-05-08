#' Multi-layered Network Construction (End-to-End)
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
  node_colors = c("Microbe" = "#D55E00", "Pathway" = "#0072B2", "Metabolite" = "#009E73"),
  node_shapes = c("Microbe" = 24, "Pathway" = 21, "Metabolite" = 22),
  base_node_size = 6, plot_width = 12, plot_height = 10, plot_dpi = 600
) {
  format <- match.arg(format)
  ppn_da_method <- match.arg(ppn_da_method)
  ppn_map_database <- match.arg(ppn_map_database)
  ppn_jaccard_method <- match.arg(ppn_jaccard_method)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Install 'data.table' for rapid preprocessing.")

  message("--- Preprocessing Data (Format: ", format, ") ---")

  processed_contrib_file <- file.path(output_dir, "processed_contribution.csv")
  taxonomy_file_to_pass <- NULL

  # --- FAST DATA.TABLE PREPROCESSING ---
  if (format == "universal") {
    df <- data.table::as.data.table(read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE))
    required <- c("SampleID", "PathwayID", "TaxonID", "contribution")
    if (!all(required %in% colnames(df))) stop("Universal format requires columns: ", paste(required, collapse = ", "))

    data.table::setnames(df, old = c("contribution", "PathwayID", "TaxonID"), new = c("taxon_function_abun", "FunctionID", "FeatureID"), skip_absent = TRUE)
    df[, TaxonID := FeatureID]
    data.table::fwrite(df, processed_contrib_file)
  } else if (format == "picrust") {
    if (is.null(taxonomy_file)) stop("Taxonomy file is required for PICRUSt format.")
    contrib <- data.table::as.data.table(read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE))
    tax <- data.table::as.data.table(read_input_file(taxonomy_file, file_type = "csv", stringsAsFactors = FALSE))

    merged <- merge(contrib, tax, by = "FeatureID", all = FALSE)
    # FIX: Swapped .() for list()
    final_df <- merged[, list(SampleID, FunctionID = PathwayID, FeatureID, TaxonID, taxon_function_abun)]
    data.table::fwrite(final_df, processed_contrib_file)
  } else if (format == "humann") {
    df <- data.table::as.data.table(read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE, check.names = FALSE))
    col1_name <- colnames(df)[1]

    long_df <- data.table::melt(df, id.vars = col1_name, variable.name = "SampleID", value.name = "taxon_function_abun")
    long_df <- long_df[grepl("\\|", get(col1_name))]
    long_df[, c("FunctionID", "FeatureID") := data.table::tstrsplit(get(col1_name), "\\|", keep = 1:2)]
    long_df[, TaxonID := FeatureID]

    # FIX: Swapped .() for list()
    final_df <- long_df[, list(SampleID, FunctionID, FeatureID, TaxonID, taxon_function_abun)]
    data.table::fwrite(final_df, processed_contrib_file)
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

  if (length(gsea_files) == 0) stop("No significant GSEA results found. Try relaxing ppn_pvalue_cutoff.")

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

    if (length(mpn_files) == 0) {
      warning("MPN layer could not be generated. Skipping integration for this comparison.")
      next
    }
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

    if (is.na(curr_mpn) || !file.exists(curr_mpn)) {
      warning("MPN file invalid. Skipping.")
      next
    }

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
