#' Multi-layered Network Construction (End-to-End)
#'
#' Standardizes input (Universal, HUMAnN, PICRUSt), builds PPN, MPN, and PMN layers,
#' and integrates them into a final visualized network.
#'
#' @param gene_abun_file Character path to gene abundance file.
#' @param path_abun_file Character path to pathway abundance file.
#' @param path_con_file Character path to pathway contribution file.
#' @param met_con_file Character path to metabolite concentration file.
#' @param metadata_file Character path to metadata file (Must contain 'SampleID' and 'class').
#' @param taxonomy_file Character path to taxonomy file (Required for 'picrust').
#' @param map_file Character path to pathway-to-gene mapping file.
#' @param output_dir Character path to output directory.
#' @param format Input format options: "universal", "humann", "picrust".
#' @param da_method DA method options: "deseq2", "edger", "maaslin2", "simple".
#' @param map_database Database options: "kegg", "metacyc".
#' @param ppn_rank_by GSEA ranking options: "signed_log_pvalue", "log2foldchange", "pvalue".
#' @param ppn_p_adjust_method GSEA p-adj options: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none".
#' @param ppn_pvalue_cutoff Numeric p-value cutoff for GSEA.
#' @param ppn_jaccard_cutoff Numeric Jaccard cutoff for PPN edges.
#' @param ppn_min_gs_size Numeric min gene set size.
#' @param ppn_max_gs_size Numeric max gene set size.
#' @param ppn_eps Numeric GSEA epsilon.
#' @param mpn_filtering MPN filtering options: "unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%".
#' @param pmn_corr_method Correlation options: "spearman", "pearson", "kendall".
#' @param pmn_filter_by Filter options: "none", "p_value", "q_value".
#' @param pmn_corr_cutoff Numeric correlation cutoff.
#' @param pmn_pvalue_cutoff Numeric p-value cutoff.
#' @param pmn_q_value_cutoff Numeric q-value cutoff.
#' @param pmn_p_adjust_method PMN q-value adj options: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none".
#' @return Vector of final network file paths.
#' @export
con_mln <- function(
    gene_abun_file,
    path_abun_file,
    path_con_file,
    met_con_file,
    metadata_file,
    taxonomy_file = NULL,
    map_file,
    output_dir,
    format = c("universal", "humann", "picrust"),
    da_method = c("deseq2", "edger", "maaslin2", "simple"),
    map_database = c("kegg", "metacyc"),
    ppn_rank_by = c("signed_log_pvalue", "log2foldchange", "pvalue"),
    ppn_p_adjust_method = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none"),
    ppn_pvalue_cutoff = 0.05,
    ppn_jaccard_cutoff = 0.2,
    ppn_min_gs_size = 10,
    ppn_max_gs_size = 500,
    ppn_eps = 1e-10,
    mpn_filtering = c("unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%"),
    pmn_corr_method = c("spearman", "pearson", "kendall"),
    pmn_filter_by = c("none", "p_value", "q_value"),
    pmn_corr_cutoff = 0.3,
    pmn_pvalue_cutoff = 0.05,
    pmn_q_value_cutoff = 0.05,
    pmn_p_adjust_method = "fdr"
) {
  format <- match.arg(format)
  da_method <- match.arg(da_method)
  map_database <- match.arg(map_database)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  message("--- Preprocessing Data (Format: ", format, ") ---")

  processed_contrib_file <- file.path(output_dir, "processed_contribution.csv")

  if (format == "universal") {
    df <- read.csv(path_con_file, stringsAsFactors = FALSE)
    required <- c("SampleID", "PathwayID", "TaxonID", "contribution")
    if (!all(required %in% colnames(df))) stop("Universal format requires columns: ", paste(required, collapse=", "))

    colnames(df)[colnames(df) == "contribution"] <- "taxon_function_abun"
    colnames(df)[colnames(df) == "PathwayID"] <- "FunctionID"
    colnames(df)[colnames(df) == "TaxonID"] <- "FeatureID"
    df$TaxonID <- df$FeatureID
    write.csv(df, processed_contrib_file, row.names = FALSE)

    processed_taxonomy_file <- file.path(output_dir, "processed_taxonomy.csv")
    dummy_tax <- data.frame(FeatureID = unique(df$FeatureID), TaxonID = unique(df$FeatureID))
    write.csv(dummy_tax, processed_taxonomy_file, row.names = FALSE)
    taxonomy_file <- processed_taxonomy_file

  } else if (format == "picrust") {
    if (is.null(taxonomy_file)) stop("Taxonomy file is required for PICRUSt format.")
    contrib <- read.csv(path_con_file, stringsAsFactors = FALSE)
    tax <- read.csv(taxonomy_file, stringsAsFactors = FALSE)

    merged <- merge(contrib, tax, by = "FeatureID")
    final_df <- merged[, c("SampleID", "PathwayID", "FeatureID", "TaxonID", "taxon_function_abun")]
    colnames(final_df)[colnames(final_df) == "PathwayID"] <- "FunctionID"
    write.csv(final_df, processed_contrib_file, row.names = FALSE)

  } else if (format == "humann") {
    df <- read.csv(path_con_file, stringsAsFactors = FALSE, check.names = FALSE)
    sample_cols <- colnames(df)[-1]
    long_df <- tidyr::pivot_longer(df, cols = dplyr::all_of(sample_cols), names_to = "SampleID", values_to = "taxon_function_abun")

    col1_name <- colnames(df)[1]
    long_df <- long_df[grep("\\|", long_df[[col1_name]]), ]
    split_data <- stringr::str_split_fixed(long_df[[col1_name]], "\\|", 2)
    long_df$FunctionID <- split_data[, 1]
    long_df$FeatureID <- split_data[, 2]
    long_df$TaxonID <- long_df$FeatureID

    final_df <- long_df[, c("SampleID", "FunctionID", "FeatureID", "TaxonID", "taxon_function_abun")]
    write.csv(final_df, processed_contrib_file, row.names = FALSE)

    processed_taxonomy_file <- file.path(output_dir, "processed_taxonomy.csv")
    dummy_tax <- data.frame(FeatureID = unique(final_df$FeatureID), TaxonID = unique(final_df$FeatureID))
    write.csv(dummy_tax, processed_taxonomy_file, row.names = FALSE)
    taxonomy_file <- processed_taxonomy_file
  }

  ppn_dir <- file.path(output_dir, "ppn_output")
  ppn_results <- con_ppn_int(
    gene_abun_file = gene_abun_file,
    metadata_file = metadata_file,
    map_file = map_file,
    output_dir = ppn_dir,
    da_method = da_method,
    map_database = map_database,
    rank_by = ppn_rank_by,
    p_adjust_method = ppn_p_adjust_method,
    pvalue_cutoff = ppn_pvalue_cutoff,
    jaccard_cutoff = ppn_jaccard_cutoff,
    min_gs_size = ppn_min_gs_size,
    max_gs_size = ppn_max_gs_size,
    eps = ppn_eps
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
      path_con_file = processed_contrib_file,
      metadata_file = metadata_file,
      taxonomy_file = taxonomy_file,
      output_dir = mpn_dir,
      mpn_filtering = mpn_filtering
    )
    curr_mpn <- mpn_files[1]

    pmn_dir <- file.path(output_dir, "pmn_output")
    pmn_files <- con_pmn_int(
      path_abun_file = path_abun_file,
      met_con_file = met_con_file,
      gsea_file = curr_gsea,
      metadata_file = metadata_file,
      output_dir = pmn_dir,
      corr_method = pmn_corr_method,
      pmn_filter_by = pmn_filter_by,
      corr_cutoff = pmn_corr_cutoff,
      pvalue_cutoff = pmn_pvalue_cutoff,
      q_value_cutoff = pmn_q_value_cutoff,
      p_adjust_method = pmn_p_adjust_method
    )
    curr_pmn <- pmn_files[1]

    mln_dir <- file.path(output_dir, "mln_final")
    final_net <- con_mln_int(
      gsea_file = curr_gsea,
      mpn_file = curr_mpn,
      ppn_file = curr_jaccard,
      pmn_file = curr_pmn,
      output_dir = mln_dir
    )
    final_networks <- c(final_networks, final_net)
  }

  message("Multi-Layered Network Pipeline Complete.")
  return(final_networks)
}
