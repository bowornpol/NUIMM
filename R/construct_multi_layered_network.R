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
#'    edge_score, edge_type).
#' 6. Combines all collected and filtered network edges into a single data frame.
#' 7. Saves the final integrated multi-layered network to a CSV file in the
#'    specified output directory, with a filename dynamically generated based
#'    on the GSEA target group or a default 'overall' suffix.
#'
#' @param gsea_results_file A character string specifying the path to the
#'   GSEA results file (from `construct_pathway_pathway_network`).
#'   This file is crucial for defining the set of pathways to be included
#'   in the multi-layered network. Must contain an 'ID' column for pathways.
#' @param microbe_pathway_file A character string specifying the path to the
#'   Microbe-Pathway network file (from `construct_microbe_pathway_network`).
#'   Expected columns: 'TaxonID', 'FunctionID', 'relative_contribution'.
#'   If the file is not found or invalid, this layer will be skipped.
#' @param pathway_jaccard_file A character string specifying the path to the
#'   Pathway-Pathway Jaccard index file (from `construct_pathway_pathway_network`).
#'   Expected columns: 'FunctionID_1', 'FunctionID_2', 'jaccard_index'.
#'   If the file is not found or invalid, this layer will be skipped.
#' @param pathway_metabolite_file A character string specifying the path to the
#'   Pathway-Metabolite network file (from `construct_pathway_metabolite_network`).
#'   Expected columns: 'FunctionID', 'MetaboliteID', 'correlation'.
#'   If the file is not found or invalid, this layer will be skipped.
#' @param output_directory A character string specifying the path to the directory
#'   where the final integrated multi-layered network CSV file will be saved.
#'   The directory will be created if it does not exist.
#' @param file_type A character string indicating the type of input files.
#'   Must be "csv" (for comma-separated) or "tsv" (for tab-separated).
#' @return The function's primary output is a single CSV file saved to the specified \code{output_directory}, containing the combined and filtered edges from all integrated network layers.
#' @export
con.mln <- function(
  gsea_results_file,
  microbe_pathway_file,
  pathway_jaccard_file,
  pathway_metabolite_file,
  output_directory,
  file_type = c("csv", "tsv")
) {
  file_type <- match.arg(file_type)
  message("Calling internal multi-layered network construction.")

  # Call the internal function. It will save files and return their paths.
  # We ignore the returned path here as per the public function's contract.
  con.mln.int( # No leading underscore
    gsea_results_file = gsea_results_file,
    microbe_pathway_file = microbe_pathway_file,
    pathway_jaccard_file = pathway_jaccard_file,
    pathway_metabolite_file = pathway_metabolite_file,
    output_directory = output_directory,
    file_type = file_type
  )

  message("Multi-layered network construction complete (file saved).")
  return(invisible(NULL))
}
