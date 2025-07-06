#' Multi-layered network construction (end-to-end)
#'
#' This comprehensive function orchestrates the construction of all individual
#' network layers (Microbe-Pathway, Pathway-Pathway, Pathway-Metabolite) and
#' then integrates them into a single, cohesive multi-layered network. It
#' combines the functionality of `construct_microbe_pathway_network`,
#' `construct_pathway_pathway_network`, `construct_pathway_metabolite_network`,
#' and `construct_multi_layered_network` into one convenient, end-to-end call.
#'
#' @details
#' The pipeline executes in the following sequential order to respect data dependencies:
#' 1.  **Pathway-Pathway Network Construction (`construct_pathway_pathway_network`):**
#'     Performs differential expression analysis and Gene Set Enrichment Analysis (GSEA)
#'     to identify enriched pathways and subsequently computes Jaccard indices
#'     to quantify overlap between core enrichment genes of significant pathways.
#'     (This step must run first as its GSEA results are a dependency for other layers).
#' 2.  **Microbe-Pathway Network Construction (`construct_microbe_pathway_network`):**
#'     Calculates the relative contributions of microbial taxa to specific functions.
#' 3.  **Pathway-Metabolite Network Construction (`construct_pathway_metabolite_network`):**
#'     Computes correlations between pathway abundances and metabolite concentrations,
#'     optionally filtering results based on GSEA findings and statistical significance.
#' 4.  **Multi-Layered Network Construction (`construct_multi_layered_network`):**
#'     Integrates the various network layers into a single comprehensive network,
#'     using each GSEA result as a filter for generating a distinct multi-layered network.
#'
#' @param abundance_file A character string specifying the path to the
#'   gene abundance data file. The first column should be gene IDs,
#'   and subsequent columns should be sample counts (integers).
#' @param metadata_file A character string specifying the path to the sample metadata file.
#'   Must include 'SampleID' and 'class' columns.
#' @param map_file A character string specifying the path to the
#'   pathway-to-gene mapping file (e.g., KEGG_pathways_to_KO.csv).
#' @param contrib_file A character string specifying the path to the contribution data file
#'   for microbe-pathway network construction. Expected columns include 'SampleID', 'FeatureID',
#'   'FunctionID', and 'taxon_function_abun'.
#' @param taxonomy_file A character string specifying the path to the taxonomy data file
#'   for microbe-pathway network construction. Expected columns include 'FeatureID' and 'TaxonID'.
#' @param pathway_abundance_file A character string specifying the path to the
#'   pathway abundance data file for pathway-metabolite network construction.
#'   Expected: Pathways as rows, samples as columns, with the first column being pathway IDs.
#' @param metabolite_concentration_file A character string specifying the path
#'   to the metabolite concentration data file for pathway-metabolite network construction.
#'   Expected: Samples as rows, metabolites as columns, with the first column being sample IDs.
#' @param base_output_directory A character string specifying the base directory
#'   where all output files will be saved. Subdirectories for each network type
#'   will be created within this directory.
#' @param file_type A character string indicating the type of input files for all steps.
#'   Must be "csv" (for comma-separated) or "tsv" (for tab-separated).
#' @param ppn_pvalueCutoff A numeric value specifying the adjusted p-value cutoff
#'   for determining significance in GSEA results within `construct_pathway_pathway_network`.
#' @param ppn_pAdjustMethod A character string specifying the method for p-value
#'   adjustment within `construct_pathway_pathway_network`.
#'   Must be one of "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none".
#' @param ppn_rank_by A character string specifying the method to rank genes for GSEA
#'   within `construct_pathway_pathway_network`.
#'   Must be either "signed_log_pvalue" or "log2FoldChange".
#' @param mpn_filtering A character string specifying the filtering method to apply
#'   for microbe-pathway network construction. Must be one of "unfiltered", "mean",
#'   "median", "top10%", "top25%", "top50%", or "top75%". "unfiltered" means no filtering is applied.
#' @param pmn_correlation_method A character string specifying the correlation method
#'   for pathway-metabolite network construction. Must be "spearman" or "pearson".
#' @param pmn_filter_by A character string specifying how to filter the correlation
#'   results for pathway-metabolite network construction. Must be "none", "p_value", or "q_value".
#' @param pmn_corr_cutoff A numeric value specifying the absolute correlation
#'   coefficient cutoff for pathway-metabolite network construction.
#' @param pmn_p_value_cutoff A numeric value specifying the p-value cutoff for
#'   pathway-metabolite network construction. Required if `pmn_filter_by` is "p_value".
#' @param pmn_q_value_cutoff A numeric value specifying the q-value (adjusted p-value)
#'   cutoff for pathway-metabolite network construction. Required if `pmn_filter_by` is "q_value".
#' @param pmn_q_adjust_method A character string specifying the method for q-value
#'   (p-value) adjustment for pathway-metabolite network construction.
#'   Must be "bonferroni" or "fdr". Used if `pmn_filter_by` is "q_value".
#' @return A character vector of paths to the final integrated multi-layered network CSV files.
#' @export
construct_multi_layered_network_full <- function(
  abundance_file,
  metadata_file,
  map_file,
  contrib_file,
  taxonomy_file,
  pathway_abundance_file,
  metabolite_concentration_file,
  base_output_directory,
  file_type = c("csv", "tsv"),
  ppn_pvalueCutoff,
  ppn_pAdjustMethod = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none"),
  ppn_rank_by = c("signed_log_pvalue", "log2FoldChange"),
  mpn_filtering = c("unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%"),
  pmn_correlation_method = c("spearman", "pearson"),
  pmn_filter_by = c("none", "p_value", "q_value"),
  pmn_corr_cutoff,
  pmn_p_value_cutoff = NULL,
  pmn_q_value_cutoff = NULL,
  pmn_q_adjust_method = NULL
) {
  file_type <- match.arg(file_type)
  ppn_pAdjustMethod <- match.arg(ppn_pAdjustMethod)
  ppn_rank_by <- match.arg(ppn_rank_by)
  mpn_filtering <- match.arg(mpn_filtering)
  pmn_correlation_method <- match.arg(pmn_correlation_method)
  pmn_filter_by <- match.arg(pmn_filter_by)
  if (!is.null(pmn_q_adjust_method)) {
    pmn_q_adjust_method <- match.arg(pmn_q_adjust_method)
  }

  # Create base output directory if it doesn't exist
  if (!dir.exists(base_output_directory)) {
    dir.create(base_output_directory, recursive = TRUE)
    message("Created base output directory: ", base_output_directory)
  }

  # Define output directories for each step
  ppn_output_dir <- file.path(base_output_directory, "pathway_pathway_network_output")
  mpn_output_dir <- file.path(base_output_directory, "microbe_pathway_network_output")
  pmn_output_dir <- file.path(base_output_directory, "pathway_metabolite_network_output")
  multi_layered_output_dir <- file.path(base_output_directory, "multi_layered_network_final")

  # Ensure all sub-directories exist
  if (!dir.exists(ppn_output_dir)) dir.create(ppn_output_dir, recursive = TRUE)
  if (!dir.exists(mpn_output_dir)) dir.create(mpn_output_dir, recursive = TRUE)
  if (!dir.exists(pmn_output_dir)) dir.create(pmn_output_dir, recursive = TRUE)
  if (!dir.exists(multi_layered_output_dir)) dir.create(multi_layered_output_dir, recursive = TRUE)


  # --- Step 1: Construct Pathway-Pathway Network ---
  message("\n--- Step 1/4: Constructing Pathway-Pathway Network ---")
  ppn_results <- construct_pathway_pathway_network_internal(
    abundance_file = abundance_file,
    metadata_file = metadata_file,
    map_file = map_file,
    output_file = ppn_output_dir,
    file_type = file_type,
    pvalueCutoff = ppn_pvalueCutoff,
    pAdjustMethod = ppn_pAdjustMethod,
    rank_by = ppn_rank_by
  )

  gsea_results_paths <- ppn_results$gsea_paths
  jaccard_results_paths <- ppn_results$jaccard_paths

  if (is.null(gsea_results_paths) || length(gsea_results_paths) == 0) {
    stop("Pathway-Pathway Network construction did not generate any GSEA results. Cannot proceed with multi-layered network construction.")
  }

  # Initialize a list to store paths of all generated multi-layered networks
  all_integrated_network_paths <- list()

  # Loop through each GSEA comparison result
  for (i in seq_along(gsea_results_paths)) {
    current_gsea_file <- gsea_results_paths[i]
    current_jaccard_file <- jaccard_results_paths[i] # Assuming Jaccard files align with GSEA files

    message(paste0("\n--- Processing GSEA comparison: ", basename(current_gsea_file), " ---"))

    # --- Step 2: Construct Microbe-Pathway Network ---
    # This step is run for each GSEA comparison to ensure relevant filtering is applied if needed.
    # Note: The current MPN function doesn't take GSEA as a direct filter.
    # If it needs to be filtered by GSEA, that logic would need to be in internal_microbe_pathway_network.R
    # or a specific MPN file should be chosen here based on the GSEA comparison.
    # For now, we assume MPN output is general or chosen outside this loop if GSEA-specific.
    message("\n--- Step 2/4: Constructing Microbe-Pathway Network ---")
    mpn_output_paths <- construct_microbe_pathway_network_internal(
      contrib_file = contrib_file,
      metadata_file = metadata_file,
      taxonomy_file = taxonomy_file,
      output_file = mpn_output_dir,
      file_type = file_type,
      filtering = mpn_filtering
    )

    if (is.null(mpn_output_paths) || length(mpn_output_paths) == 0) {
      warning("Microbe-Pathway Network construction did not generate any output files. This layer might be empty in the multi-layered network.")
      selected_microbe_pathway_file <- NULL
    } else {
      # Select the microbe-pathway file. If MPN also produces multiple files
      # (e.g., by class), you might need more complex logic here to select the
      # appropriate one for the current GSEA comparison.
      # For simplicity, if MPN output is general for all comparisons, pick the first.
      selected_microbe_pathway_file <- mpn_output_paths[1]
      message("Selected Microbe-Pathway file for integration: ", selected_microbe_pathway_file)
    }

    # --- Step 3: Construct Pathway-Metabolite Network ---
    message("\n--- Step 3/4: Constructing Pathway-Metabolite Network ---")
    pmn_output_paths <- construct_pathway_metabolite_network_internal(
      pathway_abundance_file = pathway_abundance_file,
      metabolite_concentration_file = metabolite_concentration_file,
      gsea_results_file = current_gsea_file, # Use the current GSEA file for filtering
      metadata_file = metadata_file,
      output_file = pmn_output_dir,
      file_type = file_type,
      correlation_method = pmn_correlation_method,
      filter_by = pmn_filter_by,
      corr_cutoff = pmn_corr_cutoff,
      p_value_cutoff = pmn_p_value_cutoff,
      q_value_cutoff = pmn_q_value_cutoff,
      q_adjust_method = pmn_q_adjust_method
    )

    if (is.null(pmn_output_paths) || length(pmn_output_paths) == 0) {
      warning("Pathway-Metabolite Network construction did not generate any output files for this GSEA comparison. This layer might be empty in the multi-layered network.")
      selected_pathway_metabolite_file <- NULL
    } else {
      # Select the pathway-metabolite file. It's likely that pmn_output_paths
      # will contain only one file that corresponds to the current GSEA comparison
      # due to the 'gsea_results_file' parameter. If it produces multiple (e.g., per group),
      # you might need to select the appropriate one.
      # Assuming only one relevant file is produced or the first is sufficient.
      selected_pathway_metabolite_file <- pmn_output_paths[1] # Use the first PMN file that was just generated for this GSEA
      message("Selected Pathway-Metabolite file for integration: ", selected_pathway_metabolite_file)
    }


    # --- Step 4: Construct Multi-Layered Network (Integrates all previous layers) ---
    message("\n--- Step 4/4: Constructing Multi-Layered Network for current GSEA comparison ---")
    final_multi_layered_network_path <- construct_multi_layered_network_internal( # Call internal version
      gsea_results_file = current_gsea_file, # Use the current GSEA file
      microbe_pathway_file = selected_microbe_pathway_file, # Use selected MPN file from Step 2
      pathway_jaccard_file = current_jaccard_file, # Use the current Jaccard file
      pathway_metabolite_file = selected_pathway_metabolite_file, # Use selected PMN file from Step 3
      output_directory = multi_layered_output_dir,
      file_type = file_type
    )

    if (!is.null(final_multi_layered_network_path)) {
      all_integrated_network_paths <- c(all_integrated_network_paths, final_multi_layered_network_path)
    }
  }

  message("\nMulti-layered network construction complete for all GSEA comparisons.")
  # Return all generated multi-layered network paths
  return(invisible(unlist(all_integrated_network_paths)))
}
