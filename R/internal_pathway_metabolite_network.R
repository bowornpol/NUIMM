#' Internal Pathway-Metabolite Network Helper
#' @keywords internal
#' @noRd
#' @param pmn_mode Correlation mode: "delta" (default) computes paired deltas then correlates; "pooled" uses all samples from both groups.
#' @param pmn_filter_by Significance filter: "none", "pvalue", or "padjust".
#' @param pmn_corr_cutoff Minimum absolute correlation for pathway-metabolite edges.
#' @param pmn_pvalue_cutoff P-value cutoff for pathway-metabolite edges.
#' @param pmn_padjust_cutoff Adjusted p-value cutoff for pathway-metabolite edges.
#' @param pmn_padjust_method P-value adjustment method for pathway-metabolite correlations.
con_pmn_int <- function(
  path_abun_file, met_con_file, gsea_file, metadata_file, output_dir,
  pmn_corr_method = c("spearman", "pearson", "kendall"),
  pmn_mode = c("delta", "pooled"),
  pmn_filter_by = c("none", "pvalue", "padjust"), pmn_corr_cutoff = 0.3, pmn_pvalue_cutoff = 0.05,
  pmn_padjust_cutoff = 0.05, pmn_padjust_method = "fdr",
  comparisons_list = NULL
) {
  pmn_corr_method <- match.arg(pmn_corr_method)
  pmn_mode <- match.arg(pmn_mode)
  pmn_filter_by <- match.arg(pmn_filter_by)
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
  
  meta <- meta[meta$SampleID %in% common, ]
  
  # Extract dynamic groups
  if (is.null(comparisons_list)) {
    conditions <- sort(unique(meta$class))
    if (length(conditions) >= 2) {
      baseline_grp <- conditions[1]
      treatment_grp <- conditions[2]
    } else {
      baseline_grp <- conditions[1]
      treatment_grp <- conditions[1]
    }
  } else {
    baseline_grp <- comparisons_list[[1]][1]
    treatment_grp <- comparisons_list[[1]][2]
  }

  if (pmn_mode == "pooled") {
    # Use all samples for correlation
    message("  Computing pathway-metabolite correlations (Pooled mode).")
    cor_res <- WGCNA::corAndPvalue(x = mat_x, y = mat_y, use = "pairwise.complete.obs", method = pmn_corr_method)
    
    dt_cor <- data.table::as.data.table(as.table(cor_res$cor))
    data.table::setnames(dt_cor, c("V1", "V2", "N"), c("FunctionID", "MetaboliteID", "correlation"))
    dt_pval <- data.table::as.data.table(as.table(cor_res$p))
    data.table::setnames(dt_pval, c("V1", "V2", "N"), c("FunctionID", "MetaboliteID", "p_value"))
    
    results <- merge(dt_cor, dt_pval, by = c("FunctionID", "MetaboliteID"))
    n_initial <- nrow(dt_cor)
    
  } else if (pmn_mode == "delta") {
    message("  Computing pathway-metabolite correlations (Delta mode).")
    # Identify paired samples
    v1_meta <- meta[meta$class == baseline_grp, ]
    v3_meta <- meta[meta$class == treatment_grp, ]
    
    # Strip group suffix to match subject IDs
    v1_meta$Subj <- sub(paste0("_", baseline_grp, "$"), "", v1_meta$SampleID)
    v3_meta$Subj <- sub(paste0("_", treatment_grp, "$"), "", v3_meta$SampleID)
    
    paired_subjs <- intersect(v1_meta$Subj, v3_meta$Subj)
    if (length(paired_subjs) < 3) stop("Not enough paired subjects for delta correlation.")
    
    v1_samples <- v1_meta$SampleID[match(paired_subjs, v1_meta$Subj)]
    v3_samples <- v3_meta$SampleID[match(paired_subjs, v3_meta$Subj)]
    
    delta_x <- mat_x[v3_samples, , drop=FALSE] - mat_x[v1_samples, , drop=FALSE]
    delta_y <- mat_y[v3_samples, , drop=FALSE] - mat_y[v1_samples, , drop=FALSE]
    
    cor_res <- WGCNA::corAndPvalue(x = delta_x, y = delta_y, use = "pairwise.complete.obs", method = pmn_corr_method)
    
    dt_cor <- data.table::as.data.table(as.table(cor_res$cor))
    data.table::setnames(dt_cor, c("V1", "V2", "N"), c("FunctionID", "MetaboliteID", "correlation"))
    dt_pval <- data.table::as.data.table(as.table(cor_res$p))
    data.table::setnames(dt_pval, c("V1", "V2", "N"), c("FunctionID", "MetaboliteID", "p_value"))
    
    results <- merge(dt_cor, dt_pval, by = c("FunctionID", "MetaboliteID"))
    n_initial <- nrow(dt_cor)
  }

  results <- as.data.frame(results[!is.na(results$correlation), ])

  if (nrow(results) > 0) {
    results$p_adjust <- p.adjust(results$p_value, method = pmn_padjust_method)
    results <- results[abs(results$correlation) >= pmn_corr_cutoff, ]

    if (pmn_filter_by == "pvalue") {
      results <- results[results$p_value <= pmn_pvalue_cutoff, ]
    } else if (pmn_filter_by == "padjust") {
      results <- results[results$p_adjust <= pmn_padjust_cutoff, ]
    }

    message(sprintf("    Retained %d/%d significant pathway-metabolite correlations.", nrow(results), n_initial))

    fname <- file.path(output_dir, "pmn_results.csv")
    write.csv(results, fname, row.names = FALSE)
    return(c(fname))
  } else {
    return(c())
  }
}
