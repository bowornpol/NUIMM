#' Internal Pathway-Metabolite Network Helper
#' @keywords internal
#' @noRd
#' @param pmn_method Engine: "correlation" (WGCNA) or "bdgraph" (Bayesian GCGM).
#' @param pmn_bdgraph_prior Path to pathway-metabolite prior CSV (PathwayID, MetaboliteID).
#' @param pmn_bdgraph_cutoff Minimum |delta P| for BDgraph edges. Default 0.5.
#' @param map_file Path to pathway-gene mapping (reused from GSEA for PPN priors).
#' @param pmn_mode Mode: "delta", "pooled", or "differential".
#' @param pmn_corr_method Correlation method (correlation engine only).
#' @param pmn_filter_by Significance filter (correlation engine only).
#' @param pmn_corr_cutoff Min correlation cutoff (correlation engine only).
#' @param pmn_pvalue_cutoff P-value cutoff (correlation engine only).
#' @param pmn_padjust_cutoff Adjusted p-value cutoff (correlation engine only).
#' @param pmn_padjust_method P-value adjustment method (correlation engine only).
#' @param pmn_n_perm Permutations for differential mode (correlation engine only).
con_pmn_int <- function(
  path_abun_file, met_con_file, gsea_file, metadata_file, output_dir,
  pmn_method = c("correlation", "bdgraph"),
  pmn_bdgraph_prior = NULL, pmn_bdgraph_cutoff = 0.5,
  pmn_bdgraph_iter = 5000, pmn_bdgraph_burnin = NULL,
  pmn_bdgraph_algorithm = c("bdmcmc", "rjmcmc"),
  pmn_bdgraph_method = c("gcgm", "ggm"),
  pmn_bdgraph_jump = 1, pmn_bdgraph_cores = 1, ppn_bdgraph_min_shared = 1,
  map_file = NULL,
  pmn_corr_method = c("spearman", "pearson", "kendall"),
  pmn_mode = c("delta", "pooled", "differential"),
  pmn_filter_by = c("none", "pvalue", "padjust"), pmn_corr_cutoff = 0.3, pmn_pvalue_cutoff = 0.05,
  pmn_padjust_cutoff = 0.05, pmn_padjust_method = "fdr",
  comparisons_list = NULL, pmn_n_perm = 999
) {
  pmn_method <- match.arg(pmn_method)
  pmn_bdgraph_algorithm <- match.arg(pmn_bdgraph_algorithm)
  pmn_bdgraph_method <- match.arg(pmn_bdgraph_method)
  pmn_corr_method <- match.arg(pmn_corr_method)
  pmn_mode <- match.arg(pmn_mode)
  pmn_filter_by <- match.arg(pmn_filter_by)
  if (is.null(pmn_bdgraph_burnin)) pmn_bdgraph_burnin <- round(pmn_bdgraph_iter / 2)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Load shared data
  path_abun <- read_input_file(path_abun_file, file_type = "csv", row.names = 1, check.names = FALSE)
  met_con <- read_input_file(met_con_file, file_type = "csv", row.names = 1, check.names = FALSE)
  meta <- read_input_file(metadata_file, file_type = "csv")
  gsea <- read_input_file(gsea_file, file_type = "csv")

  sig_paths <- gsea$ID
  path_abun <- path_abun[rownames(path_abun) %in% sig_paths, , drop = FALSE]

  common <- intersect(colnames(path_abun), rownames(met_con))
  common <- intersect(common, meta$SampleID)

  mat_x <- t(path_abun[, common, drop = FALSE])
  mat_y <- as.matrix(met_con[common, , drop = FALSE])
  meta <- meta[meta$SampleID %in% common, ]

  # Extract groups
  if (is.null(comparisons_list)) {
    conditions <- sort(unique(meta$class))
    baseline_grp <- conditions[1]
    treatment_grp <- if (length(conditions) >= 2) conditions[2] else conditions[1]
  } else {
    baseline_grp <- comparisons_list[[1]][1]
    treatment_grp <- comparisons_list[[1]][2]
  }
  comp_suffix <- paste0("_", baseline_grp, "_vs_", treatment_grp)

  # ============================================================
  # BDgraph Engine
  # ============================================================
  if (pmn_method == "bdgraph") {
    if (!requireNamespace("BDgraph", quietly = TRUE)) stop("Install 'BDgraph' for Bayesian structural learning.")
    message("  Constructing PMN + PPN via BDgraph (Bayesian GCGM, mode: ", pmn_mode, ").")

    pathway_names <- colnames(mat_x)
    metabolite_names <- colnames(mat_y)
    all_features <- c(pathway_names, metabolite_names)
    p <- length(all_features)
    n_path <- length(pathway_names)

    # --- Filter zero-variance features (prevents GCGM crash) ---
    var_x <- apply(mat_x, 2, var, na.rm = TRUE)
    var_y <- apply(mat_y, 2, var, na.rm = TRUE)
    zv_path <- which(is.na(var_x) | var_x == 0)
    zv_met  <- which(is.na(var_y) | var_y == 0)
    if (length(zv_path) > 0) {
      message(sprintf("    Removing %d zero-variance pathway(s): %s", length(zv_path), paste(names(zv_path), collapse = ", ")))
      mat_x <- mat_x[, -zv_path, drop = FALSE]
      pathway_names <- colnames(mat_x)
    }
    if (length(zv_met) > 0) {
      message(sprintf("    Removing %d zero-variance metabolite(s): %s", length(zv_met), paste(names(zv_met), collapse = ", ")))
      mat_y <- mat_y[, -zv_met, drop = FALSE]
      metabolite_names <- colnames(mat_y)
    }
    all_features <- c(pathway_names, metabolite_names)
    p <- length(all_features)
    n_path <- length(pathway_names)
    if (p < 3) {
      message("    Not enough features after zero-variance filtering. Skipping BDgraph.")
      return(list(pmn_path = NULL, ppn_path = NULL))
    }

    # --- Build Prior Matrix ---
    g_prior <- matrix(0.5, p, p)
    colnames(g_prior) <- rownames(g_prior) <- all_features
    diag(g_prior) <- 0  # No self-loops

    # PPN priors: pathways sharing genes get 0.9
    if (!is.null(map_file) && file.exists(map_file)) {
      map_raw <- read_input_file(map_file, file_type = "csv", header = FALSE, stringsAsFactors = FALSE)
      if (ncol(map_raw) > 2) {
        long_df <- tidyr::pivot_longer(map_raw, cols = -c(1, 2), names_to = "drop", values_to = "gene")
        gene_sets <- split(long_df$gene[!is.na(long_df$gene) & long_df$gene != ""], long_df[[1]][!is.na(long_df$gene) & long_df$gene != ""])
      } else {
        gene_sets <- split(map_raw[[2]], map_raw[[1]])
      }
      sig_gene_sets <- gene_sets[names(gene_sets) %in% pathway_names]
      pw_names <- names(sig_gene_sets)
      if (length(pw_names) > 1) {
        for (ii in 1:(length(pw_names) - 1)) {
          for (jj in (ii + 1):length(pw_names)) {
            shared <- length(intersect(sig_gene_sets[[ii]], sig_gene_sets[[jj]]))
            if (shared >= ppn_bdgraph_min_shared) {
              g_prior[pw_names[ii], pw_names[jj]] <- 0.9
              g_prior[pw_names[jj], pw_names[ii]] <- 0.9
            }
          }
        }
      }
      message(sprintf("    Injected PPN priors from map_file (%d pathways).", length(sig_gene_sets)))
    }

    # PMN priors: mapped pathway-metabolite pairs get 0.9
    if (!is.null(pmn_bdgraph_prior) && file.exists(pmn_bdgraph_prior)) {
      prior_map <- read.csv(pmn_bdgraph_prior, stringsAsFactors = FALSE)
      n_injected <- 0
      for (row_i in seq_len(nrow(prior_map))) {
        pi <- prior_map$PathwayID[row_i]
        mi <- prior_map$MetaboliteID[row_i]
        if (pi %in% all_features && mi %in% all_features) {
          g_prior[pi, mi] <- 0.9
          g_prior[mi, pi] <- 0.9
          n_injected <- n_injected + 1
        }
      }
      message(sprintf("    Injected PMN priors: %d pathway-metabolite pairs.", n_injected))
    }

    # --- Helper: run BDgraph on a sample subset ---
    run_bdgraph_subset <- function(samples) {
      joint <- cbind(mat_x[samples, , drop = FALSE], mat_y[samples, , drop = FALSE])
      fit <- BDgraph::bdgraph(data = joint, method = pmn_bdgraph_method, algorithm = pmn_bdgraph_algorithm,
                              iter = pmn_bdgraph_iter, burnin = pmn_bdgraph_burnin,
                              jump = pmn_bdgraph_jump, cores = pmn_bdgraph_cores,
                              g.prior = g_prior, verbose = FALSE)
      return(fit)
    }

    get_R <- function(fit) {
      K <- fit$K_hat
      if (is.null(K)) stop("BDgraph fit does not contain K_hat. MCMC may not have converged.")
      d <- diag(K)
      if (any(d <= 0)) {
        message("    Warning: K_hat has non-positive diagonal entries. Clamping to 1e-10.")
        d[d <= 0] <- 1e-10
      }
      D <- diag(1 / sqrt(d))
      R <- - (D %*% K %*% D)
      diag(R) <- 1
      return(R)
    }

    base_samples <- meta$SampleID[meta$class == baseline_grp]
    treat_samples <- meta$SampleID[meta$class == treatment_grp]

    if (pmn_mode == "differential") {
      message(sprintf("    Running BD-MCMC on %s (%d samples)...", baseline_grp, length(base_samples)))
      fit_base <- run_bdgraph_subset(base_samples)
      message(sprintf("    Running BD-MCMC on %s (%d samples)...", treatment_grp, length(treat_samples)))
      fit_treat <- run_bdgraph_subset(treat_samples)
      
      P_treat <- BDgraph::plinks(fit_treat); P_treat <- P_treat + t(P_treat)
      P_base  <- BDgraph::plinks(fit_base);  P_base  <- P_base  + t(P_base)
      P_mat <- P_treat - P_base
      R_mat <- get_R(fit_treat) - get_R(fit_base)
      select_mat <- abs(P_mat) >= pmn_bdgraph_cutoff

    } else if (pmn_mode == "delta") {
      v1_meta <- meta[meta$class == baseline_grp, ]
      v3_meta <- meta[meta$class == treatment_grp, ]
      v1_meta$Subj <- strip_group_suffix(v1_meta$SampleID, baseline_grp)
      v3_meta$Subj <- strip_group_suffix(v3_meta$SampleID, treatment_grp)
      paired_subjs <- intersect(v1_meta$Subj, v3_meta$Subj)
      if (length(paired_subjs) < 3) stop("Not enough paired subjects for delta BDgraph.")
      v1_s <- v1_meta$SampleID[match(paired_subjs, v1_meta$Subj)]
      v3_s <- v3_meta$SampleID[match(paired_subjs, v3_meta$Subj)]
      delta_mat <- cbind(mat_x[v3_s, , drop = FALSE] - mat_x[v1_s, , drop = FALSE],
                         mat_y[v3_s, , drop = FALSE] - mat_y[v1_s, , drop = FALSE])
      message(sprintf("    Running BD-MCMC on delta matrix (%d subjects)...", length(paired_subjs)))
      fit_delta <- BDgraph::bdgraph(data = delta_mat, method = pmn_bdgraph_method, algorithm = pmn_bdgraph_algorithm,
                                    iter = pmn_bdgraph_iter, burnin = pmn_bdgraph_burnin,
                                    jump = pmn_bdgraph_jump, cores = pmn_bdgraph_cores,
                                    g.prior = g_prior, verbose = FALSE)
      P_mat <- BDgraph::plinks(fit_delta); P_mat <- P_mat + t(P_mat)
      R_mat <- get_R(fit_delta)
      select_mat <- P_mat >= pmn_bdgraph_cutoff

    } else {
      pool_samples <- c(base_samples, treat_samples)
      message(sprintf("    Running BD-MCMC on pooled samples (%d)...", length(pool_samples)))
      fit_pool <- run_bdgraph_subset(pool_samples)
      P_mat <- BDgraph::plinks(fit_pool); P_mat <- P_mat + t(P_mat)
      R_mat <- get_R(fit_pool)
      select_mat <- P_mat >= pmn_bdgraph_cutoff
    }

    # --- Extract Edges ---
    colnames(P_mat) <- rownames(P_mat) <- all_features
    colnames(R_mat) <- rownames(R_mat) <- all_features
    pmn_list <- vector("list", p * (p - 1) / 2)
    ppn_list <- vector("list", p * (p - 1) / 2)
    pmn_k <- 0L; ppn_k <- 0L

    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        if (select_mat[i, j]) {
          ni <- all_features[i]; nj <- all_features[j]
          is_path_i <- ni %in% pathway_names; is_path_j <- nj %in% pathway_names
          is_met_i <- ni %in% metabolite_names; is_met_j <- nj %in% metabolite_names

          effect <- R_mat[i, j]
          prob <- P_mat[i, j]
          dir <- ifelse(effect > 0, "positive", "negative")

          if (is_path_i && is_met_j) {
            pmn_k <- pmn_k + 1L
            pmn_list[[pmn_k]] <- data.frame(FunctionID = ni, MetaboliteID = nj, correlation = effect,
                                            p_value = prob, direction = dir, stringsAsFactors = FALSE)
          } else if (is_met_i && is_path_j) {
            pmn_k <- pmn_k + 1L
            pmn_list[[pmn_k]] <- data.frame(FunctionID = nj, MetaboliteID = ni, correlation = effect,
                                            p_value = prob, direction = dir, stringsAsFactors = FALSE)
          } else if (is_path_i && is_path_j) {
            ppn_k <- ppn_k + 1L
            ppn_list[[ppn_k]] <- data.frame(FunctionID_1 = ni, FunctionID_2 = nj, jaccard_index = effect,
                                             delta_P = prob, direction = dir, stringsAsFactors = FALSE)
          }
        }
      }
    }
    pmn_edges <- if (pmn_k > 0) do.call(rbind, pmn_list[1:pmn_k]) else data.frame()
    ppn_edges <- if (ppn_k > 0) do.call(rbind, ppn_list[1:ppn_k]) else data.frame()

    message(sprintf("    BDgraph extracted %d PMN edges and %d PPN edges (|delta P| >= %.2f).", nrow(pmn_edges), nrow(ppn_edges), pmn_bdgraph_cutoff))

    pmn_path <- file.path(output_dir, paste0("pmn_bdgraph", comp_suffix, ".csv"))
    ppn_path <- file.path(output_dir, paste0("ppn_bdgraph", comp_suffix, ".csv"))
    write.csv(pmn_edges, pmn_path, row.names = FALSE)
    write.csv(ppn_edges, ppn_path, row.names = FALSE)

    return(list(pmn_path = pmn_path, ppn_path = ppn_path))
  }

  # ============================================================
  # Correlation Engine (Original Code — Unchanged)
  # ============================================================
  if (!requireNamespace("WGCNA", quietly = TRUE)) stop("Install 'WGCNA' for fast matrix correlations.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Install 'data.table'.")

  zv_path <- which(apply(mat_x, 2, function(col) var(col, na.rm = TRUE) == 0))
  zv_met <- which(apply(mat_y, 2, function(col) var(col, na.rm = TRUE) == 0))
  if (length(zv_path) > 0) message(sprintf("    Note: %d pathway(s) with zero variance.", length(zv_path)))
  if (length(zv_met) > 0) message(sprintf("    Note: %d metabolite(s) with zero variance.", length(zv_met)))

  if (pmn_mode == "pooled") {
    message("  Computing pathway-metabolite correlations (Pooled mode).")
    base_samples <- meta$SampleID[meta$class == baseline_grp]
    treat_samples <- meta$SampleID[meta$class == treatment_grp]
    pool_samples <- c(base_samples, treat_samples)
    cor_res <- WGCNA::corAndPvalue(x = mat_x[pool_samples, , drop=FALSE], y = mat_y[pool_samples, , drop=FALSE], use = "pairwise.complete.obs", method = pmn_corr_method)
    dt_cor <- data.table::as.data.table(as.table(cor_res$cor))
    data.table::setnames(dt_cor, c("V1", "V2", "N"), c("FunctionID", "MetaboliteID", "correlation"))
    dt_pval <- data.table::as.data.table(as.table(cor_res$p))
    data.table::setnames(dt_pval, c("V1", "V2", "N"), c("FunctionID", "MetaboliteID", "p_value"))
    results <- merge(dt_cor, dt_pval, by = c("FunctionID", "MetaboliteID"))
    n_initial <- nrow(dt_cor)

  } else if (pmn_mode == "delta") {
    message("  Computing pathway-metabolite correlations (Delta mode).")
    v1_meta <- meta[meta$class == baseline_grp, ]
    v3_meta <- meta[meta$class == treatment_grp, ]
    v1_meta$Subj <- strip_group_suffix(v1_meta$SampleID, baseline_grp)
    v3_meta$Subj <- strip_group_suffix(v3_meta$SampleID, treatment_grp)
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
    base_samples <- meta$SampleID[meta$class == baseline_grp]
    treat_samples <- meta$SampleID[meta$class == treatment_grp]
    if (length(base_samples) < 4 || length(treat_samples) < 4) stop("Not enough samples for differential correlation.")
    cor_base <- WGCNA::corAndPvalue(x = mat_x[base_samples, , drop=FALSE], y = mat_y[base_samples, , drop=FALSE], use = "pairwise.complete.obs", method = pmn_corr_method)
    cor_treat <- WGCNA::corAndPvalue(x = mat_x[treat_samples, , drop=FALSE], y = mat_y[treat_samples, , drop=FALSE], use = "pairwise.complete.obs", method = pmn_corr_method)
    bound_cor <- function(r) { r[is.na(r)] <- 0; r[r > 0.999] <- 0.999; r[r < -0.999] <- -0.999; return(r) }
    r_b <- bound_cor(cor_base$cor); r_t <- bound_cor(cor_treat$cor)
    delta_r <- r_t - r_b

    if (pmn_corr_method == "pearson") {
      message("    Using Fisher Z-test.")
      n_base <- length(base_samples); n_treat <- length(treat_samples)
      z_base <- 0.5 * log((1 + r_b) / (1 - r_b)); z_treat <- 0.5 * log((1 + r_t) / (1 - r_t))
      se_diff <- sqrt((1 / (n_base - 3)) + (1 / (n_treat - 3)))
      z_stat <- (z_treat - z_base) / se_diff
      p_diff <- 2 * (1 - pnorm(abs(z_stat)))
    } else {
      n_perm <- pmn_n_perm
      message(sprintf("    Using permutation test (%d permutations) for %s.", n_perm, pmn_corr_method))
      all_samples <- c(base_samples, treat_samples); n_b <- length(base_samples)
      mat_x_all <- mat_x[all_samples, , drop=FALSE]; mat_y_all <- mat_y[all_samples, , drop=FALSE]
      null_delta <- array(0, dim = c(ncol(mat_x_all), ncol(mat_y_all), n_perm))
      for (perm_i in seq_len(n_perm)) {
        perm_idx <- sample(length(all_samples))
        perm_base <- all_samples[perm_idx[1:n_b]]; perm_treat <- all_samples[perm_idx[(n_b+1):length(all_samples)]]
        cor_perm_b <- bound_cor(WGCNA::cor(x = mat_x_all[perm_base, , drop=FALSE], y = mat_y_all[perm_base, , drop=FALSE], use = "pairwise.complete.obs", method = pmn_corr_method))
        cor_perm_t <- bound_cor(WGCNA::cor(x = mat_x_all[perm_treat, , drop=FALSE], y = mat_y_all[perm_treat, , drop=FALSE], use = "pairwise.complete.obs", method = pmn_corr_method))
        null_delta[, , perm_i] <- cor_perm_t - cor_perm_b
      }
      obs_abs <- abs(delta_r); p_diff <- matrix(0, nrow = nrow(delta_r), ncol = ncol(delta_r))
      for (k in seq_len(n_perm)) { p_diff <- p_diff + (abs(null_delta[, , k]) >= obs_abs) }
      p_diff <- (p_diff + 1) / (n_perm + 1); dimnames(p_diff) <- dimnames(delta_r)
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
    if (pmn_filter_by == "pvalue") results <- results[results$p_value <= pmn_pvalue_cutoff, ]
    else if (pmn_filter_by == "padjust") results <- results[results$p_adjust <= pmn_padjust_cutoff, ]
    message(sprintf("    Retained %d/%d significant pathway-metabolite correlations.", nrow(results), n_initial))
    fname <- file.path(output_dir, paste0("pmn_results", comp_suffix, ".csv"))
    write.csv(results, fname, row.names = FALSE)
    return(c(fname))
  } else {
    return(c())
  }
}
