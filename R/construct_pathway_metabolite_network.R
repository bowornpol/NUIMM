#' Pathway-metabolite network construction
#'
#' This function calculates correlations between pathway abundances and metabolite
#' concentrations, optionally filtering results based on GSEA findings and
#' statistical significance. It supports both Spearman and Pearson correlation
#' methods and various p-value/q-value adjustment options.
#'
#' @details
#' The edge weights are transformed from `Edge_Score` as follows:
#' - If `Edge_Score` is `NA`, weight is `Inf`.
#' - If `Edge_Score` is less than 1, weight is `1 / Edge_Score`.
#' - If `Edge_Score` is exactly 1, weight is `1 / (Edge_Score + 0.1)` to avoid `1/1=1` which might not be ideal for shortest path.
#' - If `Edge_Score` is greater than 1, weight is `Edge_Score`.
#' This transformation aims to represent stronger connections (higher `Edge_Score`)
#' as shorter paths (lower weights).
#'
#' @param pathway_abundance_file A character string specifying the path to the
#'   pathway abundance data file. Expected: Pathways as rows,
#'   samples as columns, with the first column being pathway IDs.
#' @param metabolite_concentration_file A character string specifying the path
#'   to the metabolite concentration data file. Expected: Samples
#'   as rows, metabolites as columns, with the first column being sample IDs.
#' @param gsea_results_file An optional character string specifying the path
#'   to a GSEA results file (e.g., from `construct_pathway_pathway_network`).
#'   If provided, pathways will be filtered to those present in this file, and
#'   the filename might influence group processing. Set to `NULL` if not used.
#' @param metadata_file An optional character string specifying the path to the
#'   sample metadata file. Must include 'SampleID' and 'class'
#'   columns if provided. If `NULL`, all samples are treated as a single 'overall' group.
#' @param output_file A character string specifying the path to the directory
#'   where the output CSV files (correlation results) will be saved.
#'   The directory will be created if it does not exist.
#' @param file_type A character string indicating the type of input files.
#'   Must be "csv" (for comma-separated) or "tsv" (for tab-separated).
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
#' @return Invisible \code{NULL}. The function's primary output is CSV files saved to the specified \code{output_file} directory, containing the filtered pathway-metabolite correlation results for each processed group.
#' @export
construct_pathway_metabolite_network <- function(
  pathway_abundance_file,
  metabolite_concentration_file,
  gsea_results_file,
  metadata_file,
  output_file,
  file_type = c("csv", "tsv"),
  correlation_method = c("spearman", "pearson"),
  filter_by = c("none", "p_value", "q_value"),
  corr_cutoff,
  p_value_cutoff,
  q_value_cutoff,
  q_adjust_method = c("bonferroni", "fdr")
) {
  file_type <- match.arg(file_type)
  correlation_method <- match.arg(correlation_method)
  filter_by <- match.arg(filter_by)
  q_adjust_method <- match.arg(q_adjust_method)
  message("Calling internal pathway-metabolite network construction.")

  # Call the internal function. It will save files and return their paths.
  # We ignore the returned paths here as per the public function's contract.
  construct_pathway_metabolite_network_internal( # No leading underscore
    pathway_abundance_file = pathway_abundance_file,
    metabolite_concentration_file = metabolite_concentration_file,
    gsea_results_file = gsea_results_file,
    metadata_file = metadata_file,
    output_file = output_file,
    file_type = file_type,
    correlation_method = correlation_method,
    filter_by = filter_by,
    corr_cutoff = corr_cutoff,
    p_value_cutoff = p_value_cutoff,
    q_value_cutoff = q_value_cutoff,
    q_adjust_method = q_adjust_method
  )

  message("Pathway-metabolite network construction complete (files saved).")
  return(invisible(NULL))
}
