#' Internal Pathway-Metabolite Network Helper
#'
#' @details
#' Connects Pathways to Metabolites using vectorized, high-performance matrix
#' correlation instead of nested loops.
#'
#' @param path_abun_file Character path to pathway abundance file.
#' @param met_con_file Character path to metabolite concentration file.
#' @param gsea_file Character path to GSEA file.
#' @param metadata_file Character path to metadata file.
#' @param output_dir Character path to output directory.
#' @param pmn_corr_method Correlation options: "spearman", "pearson", "kendall".
#' @param pmn_filter_by Filter options: "none", "p_value", "q_value".
#' @param pmn_corr_cutoff Numeric correlation cutoff.
#' @param pmn_pvalue_cutoff Numeric p-value cutoff.
#' @param pmn_q_value_cutoff Numeric q-value cutoff.
#' @param pmn_p_adjust_method P-adj options: "fdr", "holm", "hochberg", etc.
#' @return Vector of PMN file paths.
#' @keywords internal
con_pmn_int <- function(
  path_abun_file, met_con_file, gsea_file, metadata_file, output_dir,
  pmn_corr_method = c("spearman", "pearson", "kendall"),
  pmn_filter_by = "p_value", pmn_corr_cutoff = 0.3, pmn_pvalue_cutoff = 0.05,
  pmn_q_value_cutoff = 0.05, pmn_p_adjust_method = "fdr"
) {
  pmn_corr_method <- match.arg(pmn_corr_method)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (!requireNamespace("WGCNA", quietly = TRUE)) stop("Install 'WGCNA' for fast matrix correlations.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Install 'data.table'.")

  path_abun <- read_input_file(path_abun_file, file_type = "csv", row.names = 1, check.names = FALSE)
  met_con <- read_input_file(met_con_file, file_type = "csv", row.names = 1, check.names = FALSE)
  meta <- read_input_file(metadata_file, file_type = "csv")
  gsea <- read_input_file(gsea_file, file_type = "csv")

  sig_paths <- gsea$ID
  path_abun <- path_abun[rownames(path_abun) %in% sig_paths, , drop = FALSE]

  common <- intersect(colnames(path_abun), rownames(met_con))
  common <- intersect(common, meta$SampleID)

  # Prepare matrices (Rows = Samples, Cols = Features)
  mat_x <- t(path_abun[, common, drop = FALSE])
  mat_y <- as.matrix(met_con[common, , drop = FALSE])

  message("  Computing fast matrix correlations...")

  # Fast vectorized correlation
  cor_res <- WGCNA::corAndPvalue(x = mat_x, y = mat_y, use = "pairwise.complete.obs", method = pmn_corr_method)

  # Fast matrix melt using data.table
  dt_cor <- data.table::as.data.table(as.table(cor_res$cor))
  data.table::setnames(dt_cor, c("V1", "V2", "N"), c("FunctionID", "MetaboliteID", "correlation"))

  dt_pval <- data.table::as.data.table(as.table(cor_res$p))
  data.table::setnames(dt_pval, c("V1", "V2", "N"), c("FunctionID", "MetaboliteID", "p_value"))

  results <- merge(dt_cor, dt_pval, by = c("FunctionID", "MetaboliteID"))
  results <- as.data.frame(results[!is.na(results$correlation), ])

  if (nrow(results) > 0) {
    results$q_value <- p.adjust(results$p_value, method = pmn_p_adjust_method)
    results <- results[abs(results$correlation) >= pmn_corr_cutoff, ]

    if (pmn_filter_by == "p_value") {
      results <- results[results$p_value <= pmn_pvalue_cutoff, ]
    } else if (pmn_filter_by == "q_value") {
      results <- results[results$q_value <= pmn_q_value_cutoff, ]
    }

    fname <- file.path(output_dir, "pmn_results.csv")
    write.csv(results, fname, row.names = FALSE)
    return(c(fname))
  } else {
    return(c())
  }
}
