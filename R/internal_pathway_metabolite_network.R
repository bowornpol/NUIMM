#' Internal Pathway-Metabolite Network Helper
#'
#' Correlates pathway abundance and metabolites.
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
#' @param pmn_p_adjust_method P-adj options: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none".
#' @return Vector of PMN file paths.
#' @keywords internal
con_pmn_int <- function(
    path_abun_file,
    met_con_file,
    gsea_file,
    metadata_file,
    output_dir,
    pmn_corr_method = c("spearman", "pearson", "kendall"),
    pmn_filter_by = "p_value",
    pmn_corr_cutoff = 0.3,
    pmn_pvalue_cutoff = 0.05,
    pmn_q_value_cutoff = 0.05,
    pmn_p_adjust_method = "fdr"
) {
  pmn_corr_method <- match.arg(pmn_corr_method)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  path_abun <- read.csv(path_abun_file, row.names = 1, check.names = FALSE)
  met_con <- read.csv(met_con_file, row.names = 1, check.names = FALSE)
  meta <- read.csv(metadata_file)
  gsea <- read.csv(gsea_file)

  sig_paths <- gsea$ID
  path_abun <- path_abun[rownames(path_abun) %in% sig_paths, ]

  common <- intersect(colnames(path_abun), rownames(met_con))
  common <- intersect(common, meta$SampleID)

  path_abun <- t(path_abun[, common])
  met_con <- met_con[common, ]
  results <- data.frame()

  for (path in colnames(path_abun)) {
    for (met in colnames(met_con)) {
      x <- as.numeric(path_abun[, path])
      y <- as.numeric(met_con[, met])
      test <- tryCatch(cor.test(x, y, method = pmn_corr_method), error = function(e) NULL)
      if (!is.null(test)) {
        results <- rbind(results, data.frame(FunctionID = path, MetaboliteID = met, correlation = test$estimate, p_value = test$p.value))
      }
    }
  }

  if(nrow(results) > 0) {
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
