#' Internal Pathway-Metabolite Network Helper
#' @keywords internal
#' @noRd
#' @param pmn_mode Correlation mode: "delta" computes paired deltas; "pooled" uses all samples;
#'   "differential" tests for correlation differences between groups.
#' @param pmn_filter_by Significance filter: "none", "pvalue", or "padjust".
#' @param pmn_corr_cutoff Minimum absolute correlation for pathway-metabolite edges.
#' @param pmn_pvalue_cutoff P-value cutoff for pathway-metabolite edges.
#' @param pmn_padjust_cutoff Adjusted p-value cutoff for pathway-metabolite edges.
#' @param pmn_padjust_method P-value adjustment method for pathway-metabolite correlations.
#' @note In differential mode, the test method is automatically selected based on
#'   `pmn_corr_method`: **Pearson** uses the Fisher Z-transformation (exact under
#'   bivariate normality); **Spearman/Kendall** uses a permutation test (999
#'   permutations, assumption-free).
con_pmn_int <- function(
  path_abun_file, met_con_file, gsea_file, metadata_file, output_dir,
  pmn_corr_method = c("spearman", "pearson", "kendall"),
  pmn_mode = c("delta", "pooled", "differential"),
  pmn_filter_by = c("none", "pvalue", "padjust"), pmn_corr_cutoff = 0.3, pmn_pvalue_cutoff = 0.05,
  pmn_padjust_cutoff = 0.05, pmn_padjust_method = "fdr",
  comparisons_list = NULL, pmn_n_perm = 999
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

  # Warn about zero-variance features (they produce NA correlations)
  zv_path <- which(apply(mat_x, 2, function(col) var(col, na.rm = TRUE) == 0))
  zv_met <- which(apply(mat_y, 2, function(col) var(col, na.rm = TRUE) == 0))
  if (length(zv_path) > 0) message(sprintf("    Note: %d pathway(s) with zero variance will yield NA correlations.", length(zv_path)))
  if (length(zv_met) > 0) message(sprintf("    Note: %d metabolite(s) with zero variance will yield NA correlations.", length(zv_met)))
  
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
    message("  Computing pathway-metabolite correlations (Pooled mode).")
    
    # Filter to only the samples relevant to this comparison
    base_samples <- meta$SampleID[meta$class == baseline_grp]
    treat_samples <- meta$SampleID[meta$class == treatment_grp]
    pool_samples <- c(base_samples, treat_samples)
    
    mat_x_pool <- mat_x[pool_samples, , drop=FALSE]
    mat_y_pool <- mat_y[pool_samples, , drop=FALSE]
    
    cor_res <- WGCNA::corAndPvalue(x = mat_x_pool, y = mat_y_pool, use = "pairwise.complete.obs", method = pmn_corr_method)
    
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
    
    # Strip group suffix to match subject IDs (use fixed matching to avoid regex issues)
    v1_meta$Subj <- sub(paste0("_", baseline_grp), "", v1_meta$SampleID, fixed = TRUE)
    v3_meta$Subj <- sub(paste0("_", treatment_grp), "", v3_meta$SampleID, fixed = TRUE)
    
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
    
  } else if (pmn_mode == "differential") {
    message("  Computing pathway-metabolite correlations (Differential mode).")
    
    # Extract samples for each group
    base_samples <- meta$SampleID[meta$class == baseline_grp]
    treat_samples <- meta$SampleID[meta$class == treatment_grp]
    
    if (length(base_samples) < 4 || length(treat_samples) < 4) stop("Not enough samples for differential correlation (need at least 4 per group).")
    
    # Calculate correlation for baseline group
    mat_x_base <- mat_x[base_samples, , drop=FALSE]
    mat_y_base <- mat_y[base_samples, , drop=FALSE]
    cor_base <- WGCNA::corAndPvalue(x = mat_x_base, y = mat_y_base, use = "pairwise.complete.obs", method = pmn_corr_method)
    
    # Calculate correlation for treatment group
    mat_x_treat <- mat_x[treat_samples, , drop=FALSE]
    mat_y_treat <- mat_y[treat_samples, , drop=FALSE]
    cor_treat <- WGCNA::corAndPvalue(x = mat_x_treat, y = mat_y_treat, use = "pairwise.complete.obs", method = pmn_corr_method)
    
    bound_cor <- function(r) {
      r[is.na(r)] <- 0
      r[r > 0.999] <- 0.999
      r[r < -0.999] <- -0.999
      return(r)
    }
    
    r_b <- bound_cor(cor_base$cor)
    r_t <- bound_cor(cor_treat$cor)
    delta_r <- r_t - r_b
    
    if (pmn_corr_method == "pearson") {
      # --- Fisher Z-transformation (exact for Pearson under bivariate normality) ---
      message("    Using Fisher Z-test (exact for Pearson correlations).")
      n_base <- length(base_samples)
      n_treat <- length(treat_samples)
      
      z_base <- 0.5 * log((1 + r_b) / (1 - r_b))
      z_treat <- 0.5 * log((1 + r_t) / (1 - r_t))
      
      se_diff <- sqrt((1 / (n_base - 3)) + (1 / (n_treat - 3)))
      z_stat <- (z_treat - z_base) / se_diff
      
      p_diff <- 2 * (1 - pnorm(abs(z_stat)))
    } else {
      # --- Permutation test (assumption-free, correct for Spearman/Kendall) ---
      n_perm <- pmn_n_perm
      message(sprintf("    Using permutation test (%d permutations) for %s correlations.", n_perm, pmn_corr_method))
      
      all_samples <- c(base_samples, treat_samples)
      n_b <- length(base_samples)
      mat_x_all <- mat_x[all_samples, , drop=FALSE]
      mat_y_all <- mat_y[all_samples, , drop=FALSE]
      
      # Observed delta is already computed as delta_r
      # Build null distribution by permuting group labels
      null_delta <- array(0, dim = c(ncol(mat_x_all), ncol(mat_y_all), n_perm))
      
      for (perm_i in seq_len(n_perm)) {
        perm_idx <- sample(length(all_samples))
        perm_base <- all_samples[perm_idx[1:n_b]]
        perm_treat <- all_samples[perm_idx[(n_b+1):length(all_samples)]]
        
        cor_perm_b <- WGCNA::cor(x = mat_x_all[perm_base, , drop=FALSE],
                                 y = mat_y_all[perm_base, , drop=FALSE],
                                 use = "pairwise.complete.obs", method = pmn_corr_method)
        cor_perm_t <- WGCNA::cor(x = mat_x_all[perm_treat, , drop=FALSE],
                                 y = mat_y_all[perm_treat, , drop=FALSE],
                                 use = "pairwise.complete.obs", method = pmn_corr_method)
        
        cor_perm_b <- bound_cor(cor_perm_b)
        cor_perm_t <- bound_cor(cor_perm_t)
        null_delta[, , perm_i] <- cor_perm_t - cor_perm_b
      }
      
      # Two-tailed empirical p-value: proportion of |null_delta| >= |observed_delta|
      obs_abs <- abs(delta_r)
      p_diff <- matrix(0, nrow = nrow(delta_r), ncol = ncol(delta_r))
      for (i in seq_len(nrow(delta_r))) {
        for (j in seq_len(ncol(delta_r))) {
          p_diff[i, j] <- (sum(abs(null_delta[i, j, ]) >= obs_abs[i, j]) + 1) / (n_perm + 1)
        }
      }
      dimnames(p_diff) <- dimnames(delta_r)
    }
    
    dt_cor <- data.table::as.data.table(as.table(delta_r))
    data.table::setnames(dt_cor, c("V1", "V2", "N"), c("FunctionID", "MetaboliteID", "correlation"))
    
    dt_pval <- data.table::as.data.table(as.table(p_diff))
    data.table::setnames(dt_pval, c("V1", "V2", "N"), c("FunctionID", "MetaboliteID", "p_value"))
    
    results <- merge(dt_cor, dt_pval, by = c("FunctionID", "MetaboliteID"))
    n_initial <- nrow(dt_cor)
  }


  results <- as.data.frame(results[!is.na(results$correlation), ])

  if (nrow(results) > 0) {
    results$p_adjust <- p.adjust(results$p_value, method = pmn_padjust_method)
    results$direction <- ifelse(results$correlation > 0, "positive", "negative")
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
