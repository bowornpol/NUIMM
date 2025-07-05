#' Multi-layered network construction
#'
#' This function integrates various network layers (Microbe-Pathway, Pathway-Pathway,
#' and Pathway-Metabolite) into a single comprehensive network. It uses GSEA
#' results as a central filter to ensure that only relevant pathways and their
#' associated connections are included in the final integrated network.
#'
#' @details
#' The function performs the following steps:
#' 1. Creates an output directory for the combined network file.
#' 2. Parses the GSEA results filename to identify a target group, which can
#'    influence the naming of the output file.
#' 3. Loads the GSEA results file to identify a definitive set of pathways
#'    to be integrated across all layers.
#' 4. Iteratively loads and processes each network layer (Microbe-Pathway,
#'    Pathway-Pathway, Pathway-Metabolite), filtering edges to include only
#'    connections involving the identified GSEA pathways.
#' 5. Standardizes the column names for each layer's edges (Feature1, Feature2,
#'    Edge_Score, Edge_Type).
#' 6. Combines all collected and filtered network edges into a single data frame.
#' 7. Saves the final integrated multi-layered network to a CSV file in the
#'    specified output directory, with a filename dynamically generated based
#'    on the GSEA target group or a default 'overall' suffix.
#'
#' @param gsea_results_file A character string specifying the path to the
#'   GSEA results file (e.g., from `construct_pathway_pathway_network`).
#'   This file is crucial for defining the set of pathways to be included
#'   in the multi-layered network. Must contain an 'ID' column for pathways.
#' @param microbe_pathway_file A character string specifying the path to the
#'   Microbe-Pathway network file (e.g., from `construct_microbe_pathway_network`).
#'   Expected columns: 'TaxonID', 'FunctionID', 'relative_contribution'.
#'   If the file is not found or invalid, this layer will be skipped.
#' @param pathway_jaccard_file A character string specifying the path to the
#'   Pathway-Pathway Jaccard index file (e.g., from `construct_pathway_pathway_network`).
#'   Expected columns: 'FunctionID_1', 'FunctionID_2', 'jaccard_index'.
#'   If the file is not found or invalid, this layer will be skipped.
#' @param pathway_metabolite_file A character string specifying the path to the
#'   Pathway-Metabolite network file (e.g., from `construct_pathway_metabolite_network`).
#'   Expected columns: 'FunctionID', 'MetaboliteID', 'Correlation'.
#'   If the file is not found or invalid, this layer will be skipped.
#' @param output_directory A character string specifying the path to the directory
#'   where the final integrated multi-layered network CSV file will be saved.
#'   The directory will be created if it does not exist.
#' @param file_type A character string indicating the type of input files.
#'   Must be "csv" (for comma-separated) or "tsv" (for tab-separated).
#' @return The function's primary output is a single CSV file
#'   saved to the specified `output_directory`, containing the combined and
#'   filtered edges from all integrated network layers.
#' @export
construct_multi_layered_network <- function(
  gsea_results_file,
  microbe_pathway_file,
  pathway_jaccard_file,
  pathway_metabolite_file,
  output_directory, # Renamed from output_file to clarify it's a directory
  file_type = c("csv", "tsv")
) {
  file_type <- match.arg(file_type)
  message("Starting multi-layered network construction.")

  # Create output directory if it doesn't exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    message("Created output directory: ", output_directory)
  } else {
    message("Output directory already exists: ", output_directory)
  }

  # Determine gsea_suffix and target groups for processing
  gsea_suffix <- NULL
  gsea_target_group_from_filename <- NULL

  if (!is.null(gsea_results_file)) {
    # Check if the GSEA file actually exists before trying to parse its name
    if (file.exists(gsea_results_file)) {
      gsea_basename <- tools::file_path_sans_ext(basename(gsea_results_file))

      # Account for additional parameters after the group names
      match_result <- stringr::str_match(gsea_basename, "^gsea_results_([^_]+)_vs_([^_]+).*$")

      if (!is.na(match_result[1,1])) { # If the full pattern matches successfully
        gsea_source_group <- match_result[1,2]
        gsea_target_group_from_filename <- match_result[1,3]

        gsea_suffix <- gsea_target_group_from_filename # Use target group as suffix for filename
        message("Derived GSEA filename suffix: '", gsea_suffix, "' (from GSEA target group)")

      } else {
        message("  Warning: GSEA results file name '", basename(gsea_results_file), "' did not strictly match the 'gsea_results_SOURCE_vs_TARGET*.csv' pattern for precise group identification. Output filename will use 'overall' suffix.")
      }
    } else {
      message("  Warning: GSEA results file '", gsea_results_file, "' not found. Output filename will use 'overall' suffix.")
    }
  } else {
    message("No GSEA results file provided, defaulting to processing without GSEA-specific group identification.")
  }

  # 1. Load GSEA results to identify all pathways to be included
  message("\n1. Identifying pathways from GSEA results for filtering...")
  if (!file.exists(gsea_results_file)) {
    stop("GSEA results file not found: '", gsea_results_file, "'. Cannot proceed.")
  }

  gsea_df <- tryCatch(
    read_input_file(gsea_results_file, file_type = file_type, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error reading GSEA results file '", gsea_results_file, "': ", e$message, sep = ""))
    }
  )

  if (!"ID" %in% colnames(gsea_df)) {
    stop("GSEA results file '", gsea_df, "' must contain an 'ID' column for pathways.")
  }

  pathways_to_integrate_set <- unique(as.character(gsea_df$ID))

  if (length(pathways_to_integrate_set) == 0) {
    stop("No pathways identified from GSEA results '", gsea_results_file, "'. Cannot construct multi-layered network.")
  }
  message("    Identified ", length(pathways_to_integrate_set), " unique pathways from GSEA results for integration.")

  all_network_edges <- list()

  # 2. Load and process Microbe-Pathway Network edges
  message("\n2. Processing Microbe-Pathway network edges...")
  if (!file.exists(microbe_pathway_file)) {
    message("    Microbe-Pathway network file not found: '", microbe_pathway_file, "'. Skipping this layer.")
  } else {
    mp_df <- tryCatch(
      read_input_file(microbe_pathway_file, file_type = file_type, stringsAsFactors = FALSE),
      error = function(e) {
        warning(paste("Error reading Microbe-Pathway file '", basename(microbe_pathway_file), "': ", e$message, ". Skipping this layer.", sep = ""))
        return(NULL)
      }
    )

    if (!is.null(mp_df) && all(c("TaxonID", "FunctionID", "relative_contribution") %in% colnames(mp_df))) {
      # Filter by pathways identified from GSEA files (FunctionID must be a GSEA pathway)
      mp_filtered <- dplyr::mutate(
        dplyr::select(
          dplyr::filter(mp_df, FunctionID %in% pathways_to_integrate_set),
          Feature1 = TaxonID,
          Feature2 = FunctionID,
          Edge_Score = relative_contribution
        ),
        Edge_Type = "Microbe-Pathway"
      )
      if (nrow(mp_filtered) > 0) {
        all_network_edges[[length(all_network_edges) + 1]] <- mp_filtered
        message("    Added ", nrow(mp_filtered), " microbe-pathway edges from ", basename(microbe_pathway_file), ".")
      } else {
        message("    No relevant microbe-pathway edges found in ", basename(microbe_pathway_file), " after filtering by GSEA pathways.")
      }
    } else {
      warning("    Microbe-Pathway file '", basename(microbe_pathway_file), "' missing required columns ('TaxonID', 'FunctionID', 'relative_contribution') or is empty. Skipping this layer.")
    }
  }

  # 3. Load and process Pathway-Pathway Network edges
  message("\n3. Processing Pathway-Pathway network edges...")
  if (!file.exists(pathway_jaccard_file)) {
    message("    Pathway-Pathway network file not found: '", pathway_jaccard_file, "'. Skipping this layer.")
  } else {
    pp_df <- tryCatch(
      read_input_file(pathway_jaccard_file, file_type = file_type, stringsAsFactors = FALSE),
      error = function(e) {
        warning(paste("Error reading Pathway-Pathway file '", basename(pathway_jaccard_file), "': ", e$message, ". Skipping this layer.", sep = ""))
        return(NULL)
      }
    )

    if (!is.null(pp_df) && all(c("FunctionID_1", "FunctionID_2", "jaccard_index") %in% colnames(pp_df))) {
      # Filter by pathways identified from GSEA files (both FunctionID_1 AND FunctionID_2 must be GSEA pathways)
      pp_filtered <- dplyr::mutate(
        dplyr::select(
          dplyr::filter(pp_df, FunctionID_1 %in% pathways_to_integrate_set & FunctionID_2 %in% pathways_to_integrate_set),
          Feature1 = FunctionID_1,
          Feature2 = FunctionID_2,
          Edge_Score = jaccard_index
        ),
        Edge_Type = "Pathway-Pathway"
      )
      if (nrow(pp_filtered) > 0) {
        all_network_edges[[length(all_network_edges) + 1]] <- pp_filtered
        message("    Added ", nrow(pp_filtered), " pathway-pathway edges from ", basename(pathway_jaccard_file), ".")
      } else {
        message("    No relevant pathway-pathway edges found in ", basename(pathway_jaccard_file), " after filtering by GSEA pathways.")
      }
    } else {
      warning("    Pathway-Pathway file '", basename(pathway_jaccard_file), "' missing required columns ('FunctionID_1', 'FunctionID_2', 'jaccard_index') or is empty. Skipping this layer.")
    }
  }

  # 4. Load and process Pathway-Metabolite Network edges
  message("\n4. Processing Pathway-Metabolite network edges...")
  if (!file.exists(pathway_metabolite_file)) {
    message("    Pathway-Metabolite network file not found: '", pathway_metabolite_file, "'. Skipping this layer.")
  } else {
    pm_df <- tryCatch(
      read_input_file(pathway_metabolite_file, file_type = file_type, stringsAsFactors = FALSE),
      error = function(e) {
        warning(paste("Error reading Pathway-Metabolite file '", basename(pathway_metabolite_file), "': ", e$message, ". Skipping this layer.", sep = ""))
        return(NULL)
      }
    )

    if (!is.null(pm_df) && all(c("FunctionID", "MetaboliteID", "Correlation") %in% colnames(pm_df))) {
      # Filter by pathways identified from GSEA files (FunctionID must be a GSEA pathway)
      pm_filtered <- dplyr::mutate(
        dplyr::select(
          dplyr::filter(pm_df, FunctionID %in% pathways_to_integrate_set),
          Feature1 = FunctionID,
          Feature2 = MetaboliteID,
          Edge_Score = Correlation
        ),
        Edge_Type = "Pathway-Metabolite"
      )
      if (nrow(pm_filtered) > 0) {
        all_network_edges[[length(all_network_edges) + 1]] <- pm_filtered
        message("    Added ", nrow(pm_filtered), " pathway-metabolite edges from ", basename(pathway_metabolite_file), ".")
      } else {
        message("    No relevant pathway-metabolite edges found in ", basename(pathway_metabolite_file), " after filtering by GSEA pathways.")
      }
    } else {
      warning("    Pathway-Metabolite file '", basename(pathway_metabolite_file), "' missing required columns ('FunctionID', 'MetaboliteID', 'Correlation') or is empty. Skipping this layer.")
    }
  }

  # 5. Combine all network edges into a single data frame
  message("\n5. Combining all network edges...")
  if (length(all_network_edges) == 0) {
    stop("No network edges were collected from any layer after filtering by GSEA pathways. Please check input files and GSEA results.")
  }

  # Define the required output columns
  final_cols <- c("Feature1", "Feature2", "Edge_Score", "Edge_Type")

  final_network_df <- dplyr::select(dplyr::bind_rows(all_network_edges), dplyr::all_of(final_cols)) # Ensure only the specified columns are kept and in order

  message("    Total integrated edges: ", nrow(final_network_df))

  # 6. Save the final integrated network with a dynamic filename
  dynamic_output_filename <- if (!is.null(gsea_suffix)) {
    paste0("multi_layered_network_", gsea_suffix, ".csv")
  } else {
    # Fallback if GSEA suffix could not be determined
    paste0("multi_layered_network_overall.csv")
  }

  final_output_filepath <- file.path(output_directory, dynamic_output_filename)
  message("\n6. Saving the final integrated multi-layered network to: ", final_output_filepath)
  write.csv(final_network_df, final_output_filepath, row.names = FALSE)
  message("Multi-layered network construction complete.")
  return(invisible(NULL))
}
