#' Pathway-metabolite network construction
#'
#' This function calculates correlations between pathway abundances and metabolite
#' concentrations, optionally filtering results based on GSEA findings and
#' statistical significance. It supports both Spearman and Pearson correlation
#' methods and various p-value/q-value adjustment options.
#'
#' The function performs the following key steps:
#' 1. Loads pathway abundance, metabolite concentration, and optional GSEA results
#'    and metadata files.
#' 2. Prepares data by aligning samples across datasets and converting pathway
#'    abundances to relative values.
#' 3. Optionally filters pathways based on identifiers found in a provided GSEA
#'    results file.
#' 4. Determines groups for correlation analysis based on metadata, with priority
#'    given to a specific group implied by the GSEA results filename if provided.
#' 5. Performs pairwise correlations between pathways and metabolites for each
#'    identified group using `stats::cor.test`.
#' 6. Calculates adjusted p-values (q-values) for all correlations within a group.
#' 7. Filters correlation results based on absolute correlation coefficient,
#'    p-value, or q-value cutoffs.
#' 8. Saves the filtered correlation results as CSV files for each processed group.
#'
#' @param pathway_abundance_file A character string specifying the path to the
#'   CSV file containing pathway abundance data. Expected: Pathways as rows,
#'   samples as columns, with the first column being pathway IDs.
#' @param metabolite_concentration_file A character string specifying the path
#'   to the CSV file containing metabolite concentration data. Expected: Samples
#'   as rows, metabolites as columns, with the first column being sample IDs.
#' @param gsea_results_file An optional character string specifying the path
#'   to a GSEA results CSV file (e.g., from `construct_pathway_pathway_network`).
#'   If provided, pathways will be filtered to those present in this file, and
#'   the filename might influence group processing. Set to `NULL` if not used.
#' @param metadata_file An optional character string specifying the path to the
#'   CSV file containing sample metadata. Must include 'SampleID' and 'class'
#'   columns if provided. If `NULL`, all samples are treated as a single 'overall' group.
#' @param output_file A character string specifying the path to the directory
#'   where the output CSV files (correlation results) will be saved.
#'   The directory will be created if it does not exist.
#' @param correlation_method A character string specifying the correlation method.
#'   Must be "spearman" or "pearson".
#' @param filter_by A character string specifying how to filter the correlation
#'   results. Must be "none", "p_value", or "q_value".
#' @param corr_cutoff A numeric value specifying the absolute correlation
#'   coefficient cutoff. Only correlations with an absolute value greater than
#'   or equal to this cutoff will be kept.
#' @param p_value_cutoff A numeric value specifying the p-value cutoff. Required
#'   if `filter_by` is "p_value".
#' @param q_value_cutoff A numeric value specifying the q-value (adjusted p-value)
#'   cutoff. Required if `filter_by` is "q_value".
#' @param q_adjust_method A character string specifying the method for q-value
#'   (p-value) adjustment. Must be "bonferroni" or "fdr" (False Discovery Rate).
#'   Used if `filter_by` is "q_value".
#' @return The function's primary output is CSV files saved
#'   to the specified `output_file` directory, containing the filtered pathway-metabolite
#'   correlation results for each processed group.
#' @importFrom dplyr %>% filter mutate across select
#' @importFrom tidyr gather
#' @importFrom stats cor.test p.adjust
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom stringr str_match
#' @importFrom tools file_path_sans_ext
#' @export
construct_pathway_metabolite_network <- function(
  pathway_abundance_file,
  metabolite_concentration_file,
  gsea_results_file,
  metadata_file,
  output_file,
  correlation_method = c("spearman", "pearson"),
  filter_by = c("none", "p_value", "q_value"),
  corr_cutoff,
  p_value_cutoff,
  q_value_cutoff,
  q_adjust_method = c("bonferroni", "fdr")
) {
  # Validate inputs
  correlation_method <- match.arg(correlation_method)
  filter_by <- match.arg(filter_by)
  q_adjust_method <- match.arg(q_adjust_method)

  message("Starting pathway-metabolite network construction.")
  message("Correlation method: ", correlation_method)
  message("Filtering results by: ", filter_by)

  if (filter_by == "p_value" && is.null(p_value_cutoff)) {
    stop("Error: 'p_value_cutoff' must be specified if 'filter_by' is 'p_value'.")
  }
  if (filter_by == "q_value" && is.null(q_value_cutoff)) {
    stop("Error: 'q_value_cutoff' must be specified if 'filter_by' is 'q_value'.")
  }

  message("Absolute correlation coefficient cutoff: ", corr_cutoff)
  if (filter_by == "p_value") message("P-value cutoff: ", p_value_cutoff)
  if (filter_by == "q_value") message("Q-value cutoff: ", q_value_cutoff, " (using ", q_adjust_method, " correction)")

  # Create output directory if it doesn't exist
  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
    message("Created output directory: ", output_file)
  } else {
    message("Output directory already exists: ", output_file)
  }

  # Determine gsea_suffix for filename and group processing priority
  gsea_suffix <- NULL
  gsea_target_group_from_filename <- NULL

  if (!is.null(gsea_results_file)) {
    # Check if the GSEA file exists before trying to read its name
    if (file.exists(gsea_results_file)) {
      gsea_basename <- tools::file_path_sans_ext(basename(gsea_results_file))

      # Account for additional parameters after the group names
      match_result <- stringr::str_match(gsea_basename, "^gsea_results_([^_]+)_vs_([^_]+).*$")

      if (!is.na(match_result[1,1])) {
        gsea_source_group <- match_result[1,2]
        gsea_target_group_from_filename <- match_result[1,3]

        gsea_suffix <- gsea_target_group_from_filename
        message("Derived GSEA filename suffix: '", gsea_suffix, "' (from GSEA target group)")

      } else {
        message("  Warning: GSEA results file name '", basename(gsea_results_file), "' did not strictly match the 'gsea_results_SOURCE_vs_TARGET*.csv' pattern. No GSEA-specific group filtering will be applied based on filename.")
      }
    } else {
      message("  Warning: GSEA results file '", gsea_results_file, "' not found. No GSEA-specific group filtering will be applied based on filename.")
    }
  } else {
    message("No GSEA results file provided. Processing will default to groups from metadata or overall if no metadata.")
  }

  # 1. Load data
  message("1. Loading pathway abundance, metabolite concentration, and GSEA results...")

  # Load Pathway Abundance (Expected: Pathways as Rows, Samples as Columns)
  pathway_abun_absolute <- tryCatch(
    read.csv(pathway_abundance_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading pathway abundance file '", pathway_abundance_file, "': ", e$message, sep = ""))
    }
  )

  # Load Metabolite Concentration (Expected: Samples as Rows, Metabolites as Columns)
  metabolite_conc_absolute <- tryCatch(
    read.csv(metabolite_concentration_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading metabolite concentration file '", metabolite_concentration_file, "': ", e$message, sep = ""))
    }
  )

  message("  Loaded pathway abundance. Dimensions: ", paste(dim(pathway_abun_absolute), collapse = "x"))
  message("  Loaded metabolite concentration. Dimensions: ", paste(dim(metabolite_conc_absolute), collapse = "x"))

  # 2. Load metadata if provided
  metadata <- NULL
  if (!is.null(metadata_file)) {
    message("2. Loading sample metadata from: ", metadata_file)
    metadata <- tryCatch(
      read.csv(metadata_file, header = TRUE, stringsAsFactors = FALSE),
      error = function(e) {
        warning(paste("Warning: Error loading metadata file '", metadata_file, "': ", e$message, ". Proceeding without group-specific analysis.", sep = ""))
        return(NULL)
      }
    )
    if (!is.null(metadata)) {
      message("  Loaded metadata. Dimensions: ", paste(dim(metadata), collapse = "x"))
      if (!"SampleID" %in% colnames(metadata)) {
        warning("Metadata file missing 'SampleID' column. Proceeding without group-specific analysis.")
        metadata <- NULL
      } else if (!"class" %in% colnames(metadata)) {
        warning("Metadata file missing 'class' column. Proceeding without group-specific analysis.")
        metadata <- NULL
      } else {
        metadata$SampleID <- as.character(metadata$SampleID)
        metadata$class <- as.factor(metadata$class)
        rownames(metadata) <- metadata$SampleID
      }
    }
  } else {
    message("2. No metadata file provided. Performing correlation on overall dataset (single group).")
  }

  # 3. Data preparation: Ensure SampleIDs align
  message("3. Preparing data for correlation...")

  # Get common samples across all datasets
  common_samples <- intersect(colnames(pathway_abun_absolute), rownames(metabolite_conc_absolute))
  if (!is.null(metadata)) {
    common_samples <- intersect(common_samples, rownames(metadata))
  }

  if (length(common_samples) == 0) {
    stop("No common samples found between pathway abundance, metabolite concentration, and metadata (if provided). Please check SampleIDs for consistency across all input files.")
  }
  message("  Number of common samples: ", length(common_samples))

  # Filter abundance data to common samples
  pathway_abun_filtered_by_samples <- pathway_abun_absolute[, common_samples, drop = FALSE]

  # metabolite_conc_filtered_by_samples will be Samples x Metabolites
  metabolite_conc_filtered_by_samples <- metabolite_conc_absolute[common_samples, , drop = FALSE]

  # Filter metadata to common samples OR create an "overall" group if no metadata
  if (!is.null(metadata)) {
    metadata_filtered <- metadata[common_samples, , drop = FALSE]
  } else {
    metadata_filtered <- data.frame(SampleID = common_samples, class = "overall", row.names = common_samples)
    metadata_filtered$class <- as.factor(metadata_filtered$class)
    message("  No metadata provided, processing all samples as a single 'overall' group.")
  }

  # Convert pathway abundance to relative values (SAMPLE-WISE normalization)
  message("  Converting pathway abundance to relative values (sample-wise) using dplyr syntax...")
  pathway_abun_temp_tibble <- pathway_abun_filtered_by_samples %>%
    tibble::rownames_to_column(var = "FunctionID")

  relative_abundance_tibble <- pathway_abun_temp_tibble %>%
    mutate(across(-FunctionID, ~ {
      col_sum <- sum(., na.rm = TRUE)
      if (col_sum == 0) {
        0
      } else {
        . / col_sum
      }
    }, .names = "{col}"))

  pathway_abun_relative <- relative_abundance_tibble %>%
    tibble::column_to_rownames(var = "FunctionID")

  pathway_abun_relative[is.na(pathway_abun_filtered_by_samples)] <- NA
  message("  Pathway abundance converted to relative values.")

  # Filter normalized pathways based on GSEA results
  pathway_abun_for_correlation <- pathway_abun_relative
  if (!is.null(gsea_results_file) && file.exists(gsea_results_file)) { # Check file existence again before loading
    message("  Filtering normalized pathways based on GSEA results...")
    gsea_results_df <- tryCatch(
      read.csv(gsea_results_file, header = TRUE, stringsAsFactors = FALSE),
      error = function(e) {
        stop(paste("Error loading GSEA results file '", gsea_results_file, "': ", e$message, sep = ""))
      }
    )

    pathway_id_col <- NULL
    if ("ID" %in% colnames(gsea_results_df)) {
      pathway_id_col <- "ID"
    } else {
      stop("GSEA results file must contain an 'ID' column with pathway identifiers.")
    }

    pathways_from_gsea <- unique(gsea_results_df[[pathway_id_col]])

    initial_normalized_pathway_count <- nrow(pathway_abun_relative)
    pathway_abun_for_correlation <- pathway_abun_relative[rownames(pathway_abun_relative) %in% pathways_from_gsea, , drop = FALSE]

    if (nrow(pathway_abun_for_correlation) == 0) {
      stop("No pathways from GSEA results were found in the normalized pathway abundance file. Please check pathway identifiers in both files.")
    }
    message("  Filtered normalized pathway abundance to ", nrow(pathway_abun_for_correlation), " pathways (from ", initial_normalized_pathway_count, " initial) based on GSEA results.")
  } else {
    message("  No GSEA results file provided or file not found. Proceeding with all normalized pathways for correlation.")
  }

  # Ensure all data are numeric for correlation (convert data frames to numeric)
  pathway_abun_for_correlation[] <- lapply(pathway_abun_for_correlation, as.numeric)
  metabolite_conc_filtered_by_samples[] <- lapply(metabolite_conc_filtered_by_samples, as.numeric)

  # Determine groups to process based on metadata AND GSEA filename priority
  all_metadata_groups <- unique(metadata_filtered$class)
  groups_to_process <- all_metadata_groups # Default: process all groups from metadata

  if (!is.null(gsea_target_group_from_filename)) { # If a target group was successfully extracted from GSEA filename
    if (gsea_target_group_from_filename %in% all_metadata_groups) {
      # If the GSEA target group exists in metadata, process ONLY that group
      groups_to_process <- gsea_target_group_from_filename
      message("  Restricting processing to group: '", gsea_target_group_from_filename, "' as implied by GSEA results filename.")
    } else {
      warning("  GSEA target group '", gsea_target_group_from_filename, "' (from filename) not found in metadata groups. Processing all groups from metadata.")
      # Fallback to all groups from metadata if the target group from filename isn't in metadata
      groups_to_process <- all_metadata_groups
    }
  } else {
    message("  Processing all groups from metadata (no specific GSEA target group identified from filename).")
  }

  if (length(groups_to_process) == 0) {
    stop("No valid groups to process after filtering based on GSEA filename and metadata. Please check group names consistency.")
  }
  message("Final groups to process: ", paste(groups_to_process, collapse = ", "))

  all_correlation_results <- list()

  # 4. Loop over each (now potentially filtered) group to perform correlations
  message("4. Performing correlations for each selected group...")
  for (current_group in groups_to_process) {
    message("  Processing group: '", current_group, "'")

    # Get samples belonging to the current group
    samples_in_group <- rownames(metadata_filtered[metadata_filtered$class == current_group, , drop = FALSE])

    if (length(samples_in_group) < 3) {
      warning("  Not enough samples (less than 3) in group '", current_group, "' to perform meaningful correlations. Skipping.")
      next
    }

    # Subset data for the current group
    path_data_group <- pathway_abun_for_correlation[, samples_in_group, drop = FALSE]
    met_data_group <- metabolite_conc_filtered_by_samples[samples_in_group, , drop = FALSE]

    # Transpose pathway data for correlation: cor.test expects vectors or Samples as Rows
    path_data_group <- t(path_data_group)

    # Get column names for pathways and metabolites
    pathway_names <- colnames(path_data_group)
    metabolite_names <- colnames(met_data_group)

    # Initialize empty list to store correlation results in long format
    results_list_for_group <- list()

    message("    Calculating correlations using cor.test for each pair...")
    # Loop through each pathway and each metabolite to get individual correlation tests
    for (path_name in pathway_names) {
      for (met_name in metabolite_names) {
        vec_path <- path_data_group[, path_name]
        vec_met <- met_data_group[, met_name]

        # Calculate number of complete observations for this specific pair
        valid_obs_count <- sum(complete.cases(vec_path, vec_met))

        # Only perform cor.test if there are at least 3 valid pairs
        if (valid_obs_count >= 3) {
          cor_args <- list(x = vec_path, y = vec_met, method = correlation_method)
          if (correlation_method == "spearman") {
            cor_args$exact <- FALSE # For large N, exact p-value is computationally intensive for Spearman
          }

          cor_test_result <- tryCatch({
            do.call(stats::cor.test, cor_args)
          }, error = function(e) {
            return(NULL) # Return NULL if error occurs
          })

          if (!is.null(cor_test_result)) {
            results_list_for_group[[length(results_list_for_group) + 1]] <- data.frame(
              FunctionID = path_name,
              MetaboliteID = met_name,
              Correlation = cor_test_result$estimate,
              P_value = cor_test_result$p.value,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    if (length(results_list_for_group) == 0) {
      message("    No valid pathway-metabolite correlations found for group '", current_group, "'. Skipping output for this group.")
      next # Skip to next group if no results
    }

    # Combine all results for the current group into a single data frame
    combined_results <- do.call(rbind, results_list_for_group)

    # Calculate Q-values using the specified adjustment method for ALL p-values in this group
    combined_results$Q_value <- stats::p.adjust(combined_results$P_value, method = q_adjust_method)

    # 5. Apply filtering based on user choice (now including Q_value)
    message("  Applying filters for group '", current_group, "'...")

    # Filter by absolute correlation coefficient first
    combined_results_filtered <- combined_results %>%
      filter(abs(Correlation) >= corr_cutoff)

    if (filter_by == "p_value") {
      combined_results_filtered <- combined_results_filtered %>%
        filter(P_value <= p_value_cutoff)
      message("  Filtered by p-value <= ", p_value_cutoff)
    } else if (filter_by == "q_value") {
      combined_results_filtered <- combined_results_filtered %>%
        filter(Q_value <= q_value_cutoff)
      message("  Filtered by q-value <= ", q_value_cutoff, " (", q_adjust_method, " correction)")
    }

    if (nrow(combined_results_filtered) == 0) {
      message("  No significant correlations found after filtering for group '", current_group, "'. Skipping output for this group.")
      next
    }

    combined_results_filtered$Group <- current_group
    all_correlation_results[[current_group]] <- combined_results_filtered

    # --- Save results with specific filename based on parameters ---

    # Determine the group name for the filename
    # Prioritize gsea_suffix if it was successfully extracted and applies to this group
    group_name_for_file <- current_group
    if (!is.null(gsea_suffix) && current_group == gsea_suffix) {
      group_name_for_file <- gsea_suffix
    }
    # If current_group is "overall" (meaning no metadata was provided), and gsea_suffix exists,
    # then use gsea_suffix in the filename. Otherwise, "overall" is correct.
    if (current_group == "overall" && !is.null(gsea_suffix)) {
      group_name_for_file <- gsea_suffix
    }

    # Format cutoffs for filename (keeping decimals)
    # Use formatC to avoid scientific notation for very small numbers, use "f" for fixed
    corr_cutoff_fname <- formatC(corr_cutoff, format = "f", digits = 2)

    # Construct the simplified output filename
    output_filename <- paste0(
      "pathway_metabolite_network_",
      group_name_for_file, "_",
      correlation_method, "_",
      corr_cutoff_fname,
      ".csv"
    )

    output_filepath <- file.path(output_file, output_filename)
    write.csv(combined_results_filtered, output_filepath, row.names = FALSE)
    message("  Saved results for group '", current_group, "' to: ", output_filepath)
  }

  message("Pathway-metabolite network construction complete.")
  return(invisible(NULL))
}
