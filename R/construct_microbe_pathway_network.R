#' Microbe-pathway network construction
#'
#' This function processes microbial contribution, metadata, and taxonomy data
#' to construct a microbe-pathway network, calculating relative contributions
#' of taxa to specific functions. It supports various filtering options
#' to focus on the most significant contributions.
#'
#' @details
#' The function performs several key steps:
#' 1. Loads input CSV/TSV files (contribution, metadata, taxonomy).
#' 2. Handles missing 'class' column in metadata by assigning a default 'all' class.
#' 3. Merges the datasets based on SampleID and FeatureID.
#' 4. Aggregates taxon-function abundance per class.
#' 5. Calculates the relative contribution of each taxon to a function.
#' 6. Applies optional filtering (mean, median, top percentages) based on
#'    relative contribution.
#' 7. Saves the processed network data as a CSV file for each class.
#'
#' @param contrib_file A character string specifying the path to the contribution data file.
#'   Expected columns include 'SampleID', 'FeatureID', 'FunctionID', and 'taxon_function_abun'.
#' @param metadata_file A character string specifying the path to the sample metadata file.
#'   Expected columns include 'SampleID' and optionally 'class'.
#' @param taxonomy_file A character string specifying the path to the taxonomy data file.
#'   Expected columns include 'FeatureID' and 'TaxonID'.
#' @param output_file A character string specifying the path to the directory
#'   where the output CSV files will be saved. The directory will be created
#'   if it does not exist.
#' @param file_type A character string indicating the type of input files.
#'   Must be "csv" (for comma-separated) or "tsv" (for tab-separated).
#' @param filtering A character string specifying the filtering method to apply.
#'   Must be one of "unfiltered", "mean", "median", "top10%", "top25%", "top50%",
#'   or "top75%". "unfiltered" means no filtering is applied.
#' @return The function’s primary output is CSV files saved to the specified \code{output_file} directory, one for each unique class.
#' @export
con_mpn <- function(
  contrib_file,
  metadata_file,
  taxonomy_file,
  output_file,
  file_type = c("csv", "tsv"),
  filtering = c("unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%")
) {
  file_type <- match.arg(file_type)
  filtering <- match.arg(filtering)
  message("Calling internal microbe-pathway network construction.")

  # Call the internal function. It will save files and return their paths.
  # We ignore the returned paths here as per the public function's contract.
  con_mpn_int( # No leading underscore
    contrib_file = contrib_file,
    metadata_file = metadata_file,
    taxonomy_file = taxonomy_file,
    output_file = output_file,
    file_type = file_type,
    filtering = filtering
  )

  message("Microbe-pathway network construction complete (files saved).")
  return(invisible(NULL))
}
