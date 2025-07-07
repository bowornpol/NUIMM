#' Pathway-pathway network construction
#'
#' This function performs differential expression analysis (using DESeq2) and
#' Gene Set Enrichment Analysis (GSEA) to identify enriched pathways and
#' subsequently computes Jaccard indices to quantify overlap between core
#' enrichment genes of significant pathways. It allows for flexible p-value
#' adjustment and gene ranking methods.
#'
#' @details
#' The overall workflow includes:
#' 1. Loading gene abundance, sample metadata, and pathway-to-gene mapping files.
#' 2. Preparing data for DESeq2, including sample alignment and rounding counts.
#' 3. Performing DESeq2 differential expression analysis for all pairwise
#'    comparisons between conditions defined in the metadata.
#' 4. Reshaping the pathway mapping file into the `clusterProfiler` TERM2GENE format.
#' 5. Running GSEA for each pairwise comparison, ranking genes by either
#'    signed log10 p-value or log2FoldChange.
#' 6. Saving the GSEA results for each comparison.
#' 7. Calculating Jaccard indices between core enrichment gene sets of significant
#'    pathways within each comparison, saving these overlap results.
#'
#' @param abundance_file A character string specifying the path to the
#'   gene abundance data file. The first column should be Gene IDs,
#'   and subsequent columns should be sample counts (integers).
#' @param metadata_file A character string specifying the path to the sample metadata file.
#'   Must include 'SampleID' and 'class' columns.
#' @param map_file A character string specifying the path to the
#'   pathway-to-gene mapping file. Expected to be a two-column CSV/TSV
#'   where the first column 'FunctionID' and all other columns containing the corresponding Gene IDs.
#' @param output_file A character string specifying the path to the directory
#'   where output CSV files (GSEA results and Jaccard indices) will be saved.
#'   The directory will be created if it does not exist.
#' @param file_type A character string indicating the type of input files.
#'   Must be "csv" (for comma-separated) or "tsv" (for tab-separated).
#' @param pvalueCutoff A numeric value specifying the adjusted p-value cutoff
#'   for determining significance in GSEA results.
#' @param pAdjustMethod A character string specifying the method for p-value
#'   adjustment. Must be one of "fdr" (False Discovery Rate), "holm", "hochberg",
#'   "hommel", "bonferroni", "BH", "BY", or "none".
#' @param rank_by A character string specifying the method to rank genes for GSEA.
#'   Must be either "signed_log_pvalue" (sign of log2FoldChange * -log10(p-value))
#'   or "log2FoldChange".
#' @return The function's primary output is CSV files saved to the specified \code{output_file} directory, containing GSEA results and pathway Jaccard indices for each comparison.
#' @export
con_ppn <- function(
  abundance_file,
  metadata_file,
  map_file,
  output_file,
  file_type = c("csv", "tsv"),
  pvalueCutoff,
  pAdjustMethod = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none"),
  rank_by = c("signed_log_pvalue", "log2FoldChange")
) {
  file_type <- match.arg(file_type)
  pAdjustMethod <- match.arg(pAdjustMethod)
  rank_by <- match.arg(rank_by)
  message("Calling internal pathway-pathway network construction.")

  # Call the internal function. It will save files and return their paths.
  # We ignore the returned paths here as per the public function's contract.
  con_ppn_int( # No leading underscore
    abundance_file = abundance_file,
    metadata_file = metadata_file,
    map_file = map_file,
    output_file = output_file,
    file_type = file_type,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    rank_by = rank_by
  )

  message("Pathway-pathway network construction complete (files saved).")
  return(invisible(NULL))
}
