utils::globalVariables(c("relative_contribution", "FunctionID", "taxon_function_abun", "total_abun", "TaxonID"))

#' Internal Microbe-Pathway Network Helper
#' @keywords internal
#' @param path_con_file Path to contribution file.
#' @param metadata_file Path to sample metadata CSV.
#' @param taxonomy_file Path to taxonomy mapping file (optional).
#' @param output_dir Path to output directory.
#' @param mpn_filtering Filtering method: "unfiltered", "mean", "median", or "topN%".
#' @param mpn_mode Mode: "delta" (default) computes paired changes (Treatment - Baseline) in per-sample relative contribution and tests significance via Wilcoxon; "pooled" combines all samples from both groups and computes relative contribution without statistical testing.
#' @param comparisons_list Optional list of pairwise group comparisons. Required for "delta" mode.

con_mpn_int <- function(
  path_con_file, metadata_file, taxonomy_file = NULL, output_dir,
  mpn_filtering = "top10%",
  mpn_mode = c("delta", "pooled"),
  comparisons_list = NULL
) {
  mpn_mode <- match.arg(mpn_mode)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Read Data
  contrib <- read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE)
  meta <- read_input_file(metadata_file, file_type = "csv", stringsAsFactors = FALSE)

  # Merge Metadata
  merged <- merge(contrib, meta, by = "SampleID")

  if (!is.null(taxonomy_file)) {
    taxonomy <- read_input_file(taxonomy_file, file_type = "csv", stringsAsFactors = FALSE)
    merged <- merge(merged, taxonomy, by = "FeatureID")
  }

  output_paths <- c()

  # --- Pooled Mode: Aggregate samples without statistical testing ---
  if (mpn_mode == "pooled") {
    message("  Analyzing Microbe-Pathway layer (Pooled mode).")

    # Compute relative contribution across all samples
    res <- merged |>
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
      fname <- file.path(output_dir, "mpn_pooled.csv")
      write.csv(res, fname, row.names = FALSE)
      output_paths <- c(output_paths, fname)
    }

  # --- Delta Mode: Paired change in contribution per subject ---
  } else if (mpn_mode == "delta") {
    message("  Analyzing Microbe-Pathway layer (Delta mode).")

    # Force 'unfiltered' mode because Wilcoxon test handles significance
    if (mpn_filtering != "unfiltered") {
      message(sprintf("    Note: mpn_filtering='%s' superseded by Wilcoxon test.", mpn_filtering))
      mpn_filtering <- "unfiltered"
    }

    # --- Dynamic Group Extraction ---
    if (is.null(comparisons_list)) {
      conditions <- sort(unique(merged$class))
      if (length(conditions) >= 2) {
        baseline_grp <- conditions[1]
        treatment_grp <- conditions[2]
      } else {
        stop("Delta mode requires at least two groups.")
      }
    } else {
      baseline_grp <- comparisons_list[[1]][1]
      treatment_grp <- comparisons_list[[1]][2]
    }

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

    base_subj <- sub(paste0("_", baseline_grp, "$"), "", base_samps)
    treat_subj <- sub(paste0("_", treatment_grp, "$"), "", treat_samps)

    paired_subjs <- intersect(base_subj, treat_subj)
    if (length(paired_subjs) < 3) stop("Not enough paired subjects for delta MPN (need at least 3).")

    message(sprintf("    Identified %d paired subjects (%s vs %s).", length(paired_subjs), baseline_grp, treatment_grp))

    # Calculate paired deltas
    message("    Computing paired subject deltas.")
    
    # Filter to paired subjects
    per_sample$SubjectID <- sub(paste0("_(", baseline_grp, "|", treatment_grp, ")$"), "", per_sample$SampleID)
    dt <- per_sample[per_sample$SubjectID %in% paired_subjs & per_sample$class %in% c(baseline_grp, treatment_grp), ]
    
    # Ensure all combinations are present to handle missing values as zeros
    unique_pairs <- unique(dt[, c("FunctionID", "TaxonID")])
    
    # Pre-allocate using expand.grid
    all_combs <- expand.grid(
      PairIdx = seq_len(nrow(unique_pairs)),
      SubjectID = paired_subjs,
      stringsAsFactors = FALSE
    )
    all_combs$FunctionID <- unique_pairs$FunctionID[all_combs$PairIdx]
    all_combs$TaxonID <- unique_pairs$TaxonID[all_combs$PairIdx]
    all_combs$PairIdx <- NULL
    
    # Merge Baseline data
    base_dt <- dt[dt$class == baseline_grp, c("FunctionID", "TaxonID", "SubjectID", "taxon_function_abun")]
    names(base_dt)[4] <- "base_val"
    all_combs <- merge(all_combs, base_dt, by = c("FunctionID", "TaxonID", "SubjectID"), all.x = TRUE)
    
    # Merge Treatment data
    treat_dt <- dt[dt$class == treatment_grp, c("FunctionID", "TaxonID", "SubjectID", "taxon_function_abun")]
    names(treat_dt)[4] <- "treat_val"
    all_combs <- merge(all_combs, treat_dt, by = c("FunctionID", "TaxonID", "SubjectID"), all.x = TRUE)
    
    # Fill NAs with 0
    all_combs$base_val[is.na(all_combs$base_val)] <- 0
    all_combs$treat_val[is.na(all_combs$treat_val)] <- 0
    
    all_combs$delta <- all_combs$treat_val - all_combs$base_val
    
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
        mean_delta = mean(delta),
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

      # Filter significant changes using unadjusted p-value < 0.05
      res_sig <- res[res$p_value < 0.05, ]

      # Also compute relative contribution magnitude for filtering
      res_sig$taxon_function_abun <- abs(res_sig$mean_delta)
      res_sig$total_abun <- ave(res_sig$taxon_function_abun, res_sig$FunctionID, FUN = sum)
      res_sig$relative_contribution <- ifelse(res_sig$total_abun == 0, 0, res_sig$taxon_function_abun / res_sig$total_abun)

      # Apply the same topN% / mean / median filter on magnitude
      res_sig <- .mpn_apply_filter(res_sig, mpn_filtering)

      if (nrow(res_sig) > 0) {
        message(sprintf("    Retained %d significant delta associations (p < 0.05, %d taxa).", nrow(res_sig), length(unique(res_sig$TaxonID))))
        fname <- file.path(output_dir, "mpn_delta.csv")
        write.csv(res_sig, fname, row.names = FALSE)
        output_paths <- c(output_paths, fname)
      } else {
        message("    No significant delta associations detected (p < 0.05).")
      }
    } else {
      message("    No non-zero deltas computed.")
    }
  }

  return(output_paths)
}


#' Internal helper: Apply MPN filtering logic
#' @keywords internal
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
