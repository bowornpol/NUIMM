utils::globalVariables(c("relative_contribution", "FunctionID", "taxon_function_abun", "total_abun", "TaxonID"))

#' Internal Microbe-Pathway Network Helper
#' @keywords internal
#' @noRd
#' @param path_con_file Path to contribution file.
#' @param metadata_file Path to sample metadata CSV.
#' @param taxonomy_file Path to taxonomy mapping file (optional).
#' @param output_dir Path to output directory.
#' @param mpn_filtering Filtering method: "unfiltered", "mean", "median", or "topN%".
#' @param mpn_mode Mode: "delta" computes paired changes; "pooled" combines all samples; "differential" computes unpaired cross-sectional differences via Wilcoxon.
#' @param mpn_filter_by Significance filter for delta/differential modes: "pvalue" or "padjust".
#' @param mpn_pvalue_cutoff P-value cutoff for significance filtering (default 0.05).
#' @param mpn_padjust_cutoff Adjusted p-value cutoff for significance filtering (default 0.05).
#' @param comparisons_list Optional list of pairwise group comparisons. Required for "delta" mode.

con_mpn_int <- function(
  path_con_file, metadata_file, taxonomy_file = NULL, output_dir,
  mpn_filtering = "top10%",
  mpn_mode = c("delta", "pooled", "differential"),
  mpn_filter_by = c("pvalue", "padjust"),
  mpn_pvalue_cutoff = 0.05,
  mpn_padjust_cutoff = 0.05,
  comparisons_list = NULL
) {
  mpn_mode <- match.arg(mpn_mode)
  mpn_filter_by <- match.arg(mpn_filter_by)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Read Data efficiently using data.table to prevent out of memory aborts
  if (requireNamespace("data.table", quietly = TRUE)) {
    contrib <- data.table::fread(path_con_file, stringsAsFactors = FALSE)
  } else {
    contrib <- read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE)
  }
  meta <- read_input_file(metadata_file, file_type = "csv", stringsAsFactors = FALSE)

  # Merge Metadata Efficiently (Prevents 10GB RAM spike from merge copy)
  contrib$class <- meta$class[match(contrib$SampleID, meta$SampleID)]
  merged <- contrib

  if (!is.null(taxonomy_file)) {
    taxonomy <- read_input_file(taxonomy_file, file_type = "csv", stringsAsFactors = FALSE)
    merged <- merge(merged, taxonomy, by = "FeatureID")
  }

  output_paths <- c()

  # --- Pooled Mode: Aggregate samples without statistical testing ---
  if (mpn_mode == "pooled") {
    message("  Analyzing Microbe-Pathway layer (Pooled mode).")

    if (is.null(comparisons_list)) {
      conditions <- sort(unique(merged$class))
      baseline_grp <- conditions[1]
      treatment_grp <- if (length(conditions) >= 2) conditions[2] else conditions[1]
    } else {
      baseline_grp <- comparisons_list[[1]][1]
      treatment_grp <- comparisons_list[[1]][2]
    }
    comp_suffix <- paste0("_", baseline_grp, "_vs_", treatment_grp)

    # Filter to only the samples relevant to this comparison
    pool_data <- merged[merged$class %in% c(baseline_grp, treatment_grp), ]

    # Compute relative contribution across all samples
    res <- pool_data |>
      dplyr::group_by(FunctionID, TaxonID) |>
      dplyr::summarise(taxon_function_abun = sum(taxon_function_abun), .groups = "drop") |>
      dplyr::group_by(FunctionID) |>
      dplyr::mutate(total_abun = sum(taxon_function_abun)) |>
      dplyr::mutate(relative_contribution = ifelse(total_abun == 0, 0, taxon_function_abun / total_abun)) |>
      dplyr::ungroup() |>
      as.data.frame()

    # Filtering Logic
    res <- .mpn_apply_filter(res, mpn_filtering)

    if (nrow(res) > 0) {
      message(sprintf("    Retained %d microbe-pathway associations (%d taxa).", nrow(res), length(unique(res$TaxonID))))
      fname <- file.path(output_dir, paste0("mpn_pooled", comp_suffix, ".csv"))
      write.csv(res, fname, row.names = FALSE)
      output_paths <- c(output_paths, fname)
    }

  # --- Delta Mode: Paired change in contribution per subject ---
  } else if (mpn_mode == "delta") {
    message("  Analyzing Microbe-Pathway layer (Delta mode).")

    # --- Dynamic Group Extraction ---
    if (is.null(comparisons_list)) {
      conditions <- sort(unique(merged$class))
      if (length(conditions) >= 2) {
        baseline_grp <- conditions[1]
        treatment_grp <- conditions[2]
      } else stop("Delta mode requires at least two groups.")
    } else {
      baseline_grp <- comparisons_list[[1]][1]
      treatment_grp <- comparisons_list[[1]][2]
    }
    comp_suffix <- paste0("_", baseline_grp, "_vs_", treatment_grp)

    # Aggregate contributions per sample
    per_sample <- merged |>
      dplyr::group_by(SampleID, FunctionID, TaxonID, class) |>
      dplyr::summarise(taxon_function_abun = sum(taxon_function_abun), .groups = "drop") |>
      as.data.frame()

    # Normalize to relative contribution per sample to account for sequencing depth
    per_sample <- per_sample |>
      dplyr::group_by(SampleID, FunctionID) |>
      dplyr::mutate(
        sample_pathway_total = sum(taxon_function_abun),
        taxon_function_abun = ifelse(sample_pathway_total == 0, 0, taxon_function_abun / sample_pathway_total)
      ) |>
      dplyr::select(-sample_pathway_total) |>
      dplyr::ungroup() |>
      as.data.frame()
    message("    Applied within-sample relative contribution normalization.")

    # Identify paired subjects
    base_samps <- unique(per_sample$SampleID[per_sample$class == baseline_grp])
    treat_samps <- unique(per_sample$SampleID[per_sample$class == treatment_grp])

    base_subj <- strip_group_suffix(base_samps, baseline_grp)
    treat_subj <- strip_group_suffix(treat_samps, treatment_grp)

    paired_subjs <- intersect(base_subj, treat_subj)
    if (length(paired_subjs) < 3) stop("Not enough paired subjects for delta MPN (need at least 3).")

    message(sprintf("    Identified %d paired subjects (%s vs %s).", length(paired_subjs), baseline_grp, treatment_grp))

    # Calculate paired deltas
    message("    Computing paired subject deltas.")
    
    # Filter to paired subjects
    # Strip group suffix from each sample to derive SubjectID
    per_sample$SubjectID <- strip_group_suffix(
      strip_group_suffix(per_sample$SampleID, treatment_grp),
      baseline_grp
    )
    dt <- per_sample[per_sample$SubjectID %in% paired_subjs & per_sample$class %in% c(baseline_grp, treatment_grp), ]
    
    # Ensure all combinations are present to handle missing values as zeros
    unique_pairs <- unique(dt[, c("FunctionID", "TaxonID")])
    
    # Pre-allocate using tidyr::complete for extreme memory efficiency
    # This avoids the out of memory crash caused by expand.grid on massive datasets
    all_combs <- dt |>
      tidyr::complete(
        tidyr::nesting(FunctionID, TaxonID),
        SubjectID = paired_subjs,
        class = c(baseline_grp, treatment_grp),
        fill = list(taxon_function_abun = 0)
      ) |>
      as.data.frame()
      
    # Pivot wider to get base and treat in separate columns
    wide_combs <- tidyr::pivot_wider(
      all_combs, 
      id_cols = c("FunctionID", "TaxonID", "SubjectID"),
      names_from = "class",
      values_from = "taxon_function_abun",
      values_fill = 0
    )
    
    # Rename for logic compatibility
    names(wide_combs)[names(wide_combs) == baseline_grp] <- "base_val"
    names(wide_combs)[names(wide_combs) == treatment_grp] <- "treat_val"
    
    wide_combs$delta <- wide_combs$treat_val - wide_combs$base_val
    all_combs <- wide_combs
    
    message("    Executing Wilcoxon signed-rank tests.")
    
    # Wilcoxon wrapper
    fast_wilcox <- function(d) {
      if (all(d == 0)) return(NA_real_)
      tryCatch(wilcox.test(d, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
    }
    
    # Compute stats per pair using dplyr
    res <- all_combs |>
      dplyr::group_by(FunctionID, TaxonID) |>
      dplyr::summarise(
        median_delta = stats::median(delta),
        p_value = fast_wilcox(delta),
        n_subjects = length(paired_subjs),
        .groups = "drop"
      ) |>
      as.data.frame()
    
    # Filter out NAs (pairs with all 0 delta or test failure)
    res <- res[!is.na(res$p_value), ]

    if (nrow(res) > 0) {
      res$p_adjust <- stats::p.adjust(res$p_value, method = "fdr")

      # Filter significant changes using the user-selected criterion
      if (mpn_filter_by == "padjust") {
        res_sig <- res[res$p_adjust < mpn_padjust_cutoff, ]
        filter_label <- "q"
        cutoff_val <- mpn_padjust_cutoff
      } else {
        res_sig <- res[res$p_value < mpn_pvalue_cutoff, ]
        filter_label <- "p"
        cutoff_val <- mpn_pvalue_cutoff
      }

      if (nrow(res_sig) > 0) {
        # Compute actual relative contribution from pooled abundance (both timepoints)
        pooled_contrib <- dt[dt$class %in% c(baseline_grp, treatment_grp), ] |>
          dplyr::group_by(FunctionID, TaxonID) |>
          dplyr::summarise(taxon_function_abun = sum(taxon_function_abun), .groups = "drop") |>
          dplyr::group_by(FunctionID) |>
          dplyr::mutate(total_abun = sum(taxon_function_abun)) |>
          dplyr::mutate(relative_contribution = ifelse(total_abun == 0, 0, taxon_function_abun / total_abun)) |>
          dplyr::ungroup() |>
          as.data.frame()

        # Merge actual relative contribution into significant results
        res_sig <- merge(res_sig, pooled_contrib[, c("FunctionID", "TaxonID", "taxon_function_abun", "total_abun", "relative_contribution")],
                         by = c("FunctionID", "TaxonID"), all.x = TRUE)
        res_sig$relative_contribution[is.na(res_sig$relative_contribution)] <- 0

        # Apply the same topN% / mean / median filter on magnitude
        res_sig <- .mpn_apply_filter(res_sig, mpn_filtering)
      }

      if (nrow(res_sig) > 0) {
        message(sprintf("    Retained %d significant delta associations (%s < %g, %d taxa).", nrow(res_sig), filter_label, cutoff_val, length(unique(res_sig$TaxonID))))
        fname <- file.path(output_dir, paste0("mpn_delta", comp_suffix, ".csv"))
        write.csv(res_sig, fname, row.names = FALSE)
        output_paths <- c(output_paths, fname)
      } else {
        message(sprintf("    No significant delta associations detected (%s < %g).", filter_label, cutoff_val))
      }
    } else {
      message("    No non-zero deltas computed.")
    }
    
  # --- Differential Mode: Unpaired change for cross-sectional data ---
  } else if (mpn_mode == "differential") {
    message("  Analyzing Microbe-Pathway layer (Differential mode).")
    
    if (is.null(comparisons_list)) {
      conditions <- sort(unique(merged$class))
      if (length(conditions) >= 2) {
        baseline_grp <- conditions[1]
        treatment_grp <- conditions[2]
      } else stop("Differential mode requires at least two groups.")
    } else {
      baseline_grp <- comparisons_list[[1]][1]
      treatment_grp <- comparisons_list[[1]][2]
    }
    comp_suffix <- paste0("_", baseline_grp, "_vs_", treatment_grp)
    
    per_sample <- merged |>
      dplyr::group_by(SampleID, FunctionID, TaxonID, class) |>
      dplyr::summarise(taxon_function_abun = sum(taxon_function_abun), .groups = "drop") |>
      dplyr::group_by(SampleID, FunctionID) |>
      dplyr::mutate(
        sample_pathway_total = sum(taxon_function_abun),
        taxon_function_abun = ifelse(sample_pathway_total == 0, 0, taxon_function_abun / sample_pathway_total)
      ) |>
      dplyr::select(-sample_pathway_total) |>
      dplyr::ungroup() |>
      as.data.frame()
      
    unique_pairs <- unique(per_sample[, c("FunctionID", "TaxonID")])
    all_samps <- unique(merged[merged$class %in% c(baseline_grp, treatment_grp), c("SampleID", "class")])
    
    # Use tidyr::complete to implicitly fill missing sample-pathway-taxon combinations with 0
    # This avoids a massive expand.grid and merge that causes Out-Of-Memory crashes
    dt <- per_sample[per_sample$class %in% c(baseline_grp, treatment_grp), ]
    
    dt <- dt |>
      tidyr::complete(
        tidyr::nesting(FunctionID, TaxonID),
        SampleID = all_samps$SampleID,
        fill = list(taxon_function_abun = 0)
      ) |>
      as.data.frame()
      
    # Restore the class column for the newly completed rows
    dt$class <- all_samps$class[match(dt$SampleID, all_samps$SampleID)]
    
    fast_unpaired_wilcox <- function(val, cls) {
      if (length(unique(cls)) < 2) return(NA_real_)
      tryCatch(wilcox.test(val ~ cls, exact = FALSE)$p.value, error = function(e) NA_real_)
    }
    
    message("    Executing Unpaired Wilcoxon rank-sum tests.")
    res <- dt |>
      dplyr::group_by(FunctionID, TaxonID) |>
      dplyr::summarise(
        median_base = stats::median(taxon_function_abun[class == baseline_grp], na.rm=TRUE),
        median_treat = stats::median(taxon_function_abun[class == treatment_grp], na.rm=TRUE),
        p_value = fast_unpaired_wilcox(taxon_function_abun, class),
        .groups = "drop"
      ) |>
      dplyr::mutate(median_diff = median_treat - median_base) |>
      as.data.frame()
      
    res <- res[!is.na(res$p_value), ]
    if (nrow(res) > 0) {
      res$p_adjust <- stats::p.adjust(res$p_value, method = "fdr")

      # Filter significant changes using the user-selected criterion
      if (mpn_filter_by == "padjust") {
        res_sig <- res[res$p_adjust < mpn_padjust_cutoff, ]
        filter_label <- "q"
        cutoff_val <- mpn_padjust_cutoff
      } else {
        res_sig <- res[res$p_value < mpn_pvalue_cutoff, ]
        filter_label <- "p"
        cutoff_val <- mpn_pvalue_cutoff
      }

      if (nrow(res_sig) > 0) {
        # Compute actual relative contribution from pooled abundance (both groups)
        pooled_contrib <- dt[dt$class %in% c(baseline_grp, treatment_grp), ] |>
          dplyr::group_by(FunctionID, TaxonID) |>
          dplyr::summarise(taxon_function_abun = sum(taxon_function_abun), .groups = "drop") |>
          dplyr::group_by(FunctionID) |>
          dplyr::mutate(total_abun = sum(taxon_function_abun)) |>
          dplyr::mutate(relative_contribution = ifelse(total_abun == 0, 0, taxon_function_abun / total_abun)) |>
          dplyr::ungroup() |>
          as.data.frame()

        res_sig <- merge(res_sig, pooled_contrib[, c("FunctionID", "TaxonID", "taxon_function_abun", "total_abun", "relative_contribution")],
                         by = c("FunctionID", "TaxonID"), all.x = TRUE)
        res_sig$relative_contribution[is.na(res_sig$relative_contribution)] <- 0
        res_sig <- .mpn_apply_filter(res_sig, mpn_filtering)
      }
      
      if (nrow(res_sig) > 0) {
        message(sprintf("    Retained %d significant differential associations (%s < %g, %d taxa).", nrow(res_sig), filter_label, cutoff_val, length(unique(res_sig$TaxonID))))
        fname <- file.path(output_dir, paste0("mpn_differential", comp_suffix, ".csv"))
        write.csv(res_sig, fname, row.names = FALSE)
        output_paths <- c(output_paths, fname)
      } else message(sprintf("    No significant differential associations detected (%s < %g).", filter_label, cutoff_val))
    } else message("    No differential statistics computed.")
  }

  return(output_paths)
}


#' Internal helper: Apply MPN filtering logic
#' @keywords internal
#' @noRd
.mpn_apply_filter <- function(res, mpn_filtering) {
  if (mpn_filtering != "unfiltered") {
    if (mpn_filtering %in% c("mean", "median")) {
      FUN_used <- if (mpn_filtering == "mean") mean else median
      thresh <- aggregate(relative_contribution ~ FunctionID, res, FUN_used)
      colnames(thresh)[2] <- "threshold"
      res <- merge(res, thresh, by = "FunctionID")
      res <- res[res$relative_contribution >= res$threshold, ]
      res$threshold <- NULL
    } else if (grepl("top", mpn_filtering)) {
      perc <- as.numeric(gsub("top|%", "", mpn_filtering)) / 100
      res <- res |>
        dplyr::group_by(FunctionID) |>
        dplyr::arrange(dplyr::desc(relative_contribution)) |>
        dplyr::slice_head(prop = perc) |>
        dplyr::ungroup() |>
        as.data.frame()
    }
  }
  return(res)
}
