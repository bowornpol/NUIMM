#' Multi-layered network construction (end-to-end integration)
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
#' 4.  **Multi-Layered Network Integration (`construct_multi_layered_network`):**
#'     Combines the outputs from the previous three steps into a single, integrated
#'     network, filtering based on the GSEA results from Step 1.
#'
#' All intermediate and final output files are saved within organized subdirectories
#' of the specified `base_output_directory`.
#'
#' @param abundance_file A character string specifying the path to the gene
#'   abundance data file (input for Pathway-Pathway Network).
#' @param metadata_file A character string specifying the path to the sample
#'   metadata file (input for Pathway-Pathway, Microbe-Pathway, Pathway-Metabolite Networks).
#' @param map_file A character string specifying the path to the pathway-to-gene
#'   mapping file (input for Pathway-Pathway Network).
#' @param contrib_file A character string specifying the path to the microbial
#'   contribution data file (input for Microbe-Pathway Network).
#' @param taxonomy_file A character string specifying the path to the taxonomy
#'   data file (input for Microbe-Pathway Network).
#' @param pathway_abundance_file A character string specifying the path to the
#'   pathway abundance data file (input for Pathway-Metabolite Network).
#' @param metabolite_concentration_file A character string specifying the path
#'   to the metabolite concentration data file (input for Pathway-Metabolite Network).
#' @param base_output_directory A character string specifying the root directory
#'   where all intermediate and final output files will be saved. Subdirectories
#'   will be created within this path (e.g., "pathway_pathway_results",
#'   "microbe_pathway_results", "pathway_metabolite_results",
#'   "multi_layered_network_final").
#' @param file_type A character string indicating the type of all input files.
#'   Must be "csv" (for comma-separated) or "tsv" (for tab-separated).
#' @param ppn_pvalueCutoff A numeric value: adjusted p-value cutoff for GSEA
#'   in `construct_pathway_pathway_network`.
#' @param ppn_pAdjustMethod A character string: p-value adjustment method for GSEA
#'   in `construct_pathway_pathway_network`.
#' @param ppn_rank_by A character string: gene ranking method for GSEA
#'   in `construct_pathway_pathway_network`.
#' @param mpn_filtering A character string: filtering method for
#'   `construct_microbe_pathway_network`.
#' @param pmn_correlation_method A character string: correlation method for
#'   `construct_pathway_metabolite_network`.
#' @param pmn_filter_by A character string: how to filter correlations in
#'   `construct_pathway_metabolite_network`.
#' @param pmn_corr_cutoff A numeric value: absolute correlation coefficient cutoff
#'   for `construct_pathway_metabolite_network`.
#' @param pmn_p_value_cutoff A numeric value: p-value cutoff for
#'   `construct_pathway_metabolite_network` (required if `pmn_filter_by` is "p_value").
#' @param pmn_q_value_cutoff A numeric value: q-value cutoff for
#'   `construct_pathway_metabolite_network` (required if `pmn_filter_by` is "q_value").
#' @param pmn_q_adjust_method A character string: q-value adjustment method for
#'   `construct_pathway_metabolite_network`.
#' @return A character string representing the path to the final integrated
#'   multi-layered network CSV file.
#' @export
construct_multi_layered_network_full <- function(
  # Input files
  abundance_file, metadata_file, map_file,
  contrib_file, taxonomy_file,
  pathway_abundance_file, metabolite_concentration_file,

  # Output directory for all results
  base_output_directory,

  # Common file type for all inputs
  file_type = c("csv", "tsv"),

  # Parameters for construct_pathway_pathway_network
  ppn_pvalueCutoff,
  ppn_pAdjustMethod = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none"),
  ppn_rank_by = c("signed_log_pvalue", "log2FoldChange"),

  # Parameters for construct_microbe_pathway_network
  mpn_filtering = c("unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%"),

  # Parameters for construct_pathway_metabolite_network
  pmn_correlation_method = c("spearman", "pearson"),
  pmn_filter_by = c("none", "p_value", "q_value"),
  pmn_corr_cutoff,
  pmn_p_value_cutoff = NULL,
  pmn_q_value_cutoff = NULL,
  pmn_q_adjust_method = c("bonferroni", "fdr")
) {
  # Validate common file_type argument
  file_type <- match.arg(file_type)
  ppn_pAdjustMethod <- match.arg(ppn_pAdjustMethod)
  ppn_rank_by <- match.arg(ppn_rank_by)
  mpn_filtering <- match.arg(mpn_filtering)
  pmn_correlation_method <- match.arg(pmn_correlation_method)
  pmn_filter_by <- match.arg(pmn_filter_by)
  pmn_q_adjust_method <- match.arg(pmn_q_adjust_method)


  message("Starting full multi-layered network analysis pipeline.")
  message("All intermediate and final outputs will be saved under: ", base_output_directory)
  message("Input file type for all files: ", file_type)

  # Define subdirectories for outputs
  ppn_output_dir <- file.path(base_output_directory, "pathway_pathway_results")
  mpn_output_dir <- file.path(base_output_directory, "microbe_pathway_results")
  pmn_output_dir <- file.path(base_output_directory, "pathway_metabolite_results")
  mln_output_dir <- file.path(base_output_directory, "multi_layered_network_final")

  # Ensure base output directory exists
  if (!dir.exists(base_output_directory)) {
    dir.create(base_output_directory, recursive = TRUE)
    message("Created base output directory: ", base_output_directory)
  } else {
    message("Output directory already exists: ", base_output_directory)
  }

  # --- Step 1: Construct Pathway-Pathway Network (Dependency for PMN and MLN) ---
  message("\n--- Step 1/4: Constructing Pathway-Pathway Network ---")
  ppn_results <- construct_pathway_pathway_network_internal( # Call internal version
    abundance_file = abundance_file,
    metadata_file = metadata_file,
    map_file = map_file,
    output_file = ppn_output_dir,
    file_type = file_type,
    pvalueCutoff = ppn_pvalueCutoff,
    pAdjustMethod = ppn_pAdjustMethod,
    rank_by = ppn_rank_by
  )

  # Ensure GSEA results were generated and pick one for downstream use
  if (is.null(ppn_results$gsea_paths) || length(ppn_results$gsea_paths) == 0) {
    stop("Pathway-Pathway Network construction failed to generate any GSEA results. Cannot proceed with multi-layered integration.")
  }
  # For simplicity, we'll use the first GSEA results file generated.
  # If specific GSEA comparison is needed, the user should call functions separately.
  selected_gsea_file <- ppn_results$gsea_paths[1]
  message("Selected GSEA results file for integration: ", selected_gsea_file)

  # Ensure Jaccard results were generated
  if (is.null(ppn_results$jaccard_paths) || length(ppn_results$jaccard_paths) == 0) {
    warning("Pathway-Pathway Network construction did not generate any Jaccard index files. This layer might be empty in the multi-layered network.")
    selected_jaccard_file <- NULL # Set to NULL if no Jaccard file
  } else {
    selected_jaccard_file <- ppn_results$jaccard_paths[1] # Use the first Jaccard file
    message("Selected Pathway-Pathway Jaccard file for integration: ", selected_jaccard_file)
  }


  # --- Step 2: Construct Microbe-Pathway Network ---
  message("\n--- Step 2/4: Constructing Microbe-Pathway Network ---")
  mpn_output_paths <- construct_microbe_pathway_network_internal( # Call internal version
    contrib_file = contrib_file,
    metadata_file = metadata_file,
    taxonomy_file = taxonomy_file,
    output_file = mpn_output_dir,
    file_type = file_type,
    filtering = mpn_filtering
  )

  if (is.null(mpn_output_paths) || length(mpn_output_paths) == 0) {
    warning("Microbe-Pathway Network construction did not generate any output files. This layer might be empty in the multi-layered network.")
    selected_microbe_pathway_file <- NULL # Set to NULL if no MPN file
  } else {
    selected_microbe_pathway_file <- mpn_output_paths[1] # Use the first MPN file
    message("Selected Microbe-Pathway file for integration: ", selected_microbe_pathway_file)
  }


  # --- Step 3: Construct Pathway-Metabolite Network (Depends on GSEA results from Step 1) ---
  message("\n--- Step 3/4: Constructing Pathway-Metabolite Network ---")
  pmn_output_paths <- construct_pathway_metabolite_network_internal( # Call internal version
    pathway_abundance_file = pathway_abundance_file,
    metabolite_concentration_file = metabolite_concentration_file,
    gsea_results_file = selected_gsea_file, # Use the selected GSEA file from Step 1
    metadata_file = metadata_file, # Pass metadata to PMN for group filtering if applicable
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
    warning("Pathway-Metabolite Network construction did not generate any output files. This layer might be empty in the multi-layered network.")
    selected_pathway_metabolite_file <- NULL # Set to NULL if no PMN file
  } else {
    selected_pathway_metabolite_file <- pmn_output_paths[1] # Use the first PMN file
    message("Selected Pathway-Metabolite file for integration: ", selected_pathway_metabolite_file)
  }


  # --- Step 4: Construct Multi-Layered Network (Integrates all previous layers) ---
  message("\n--- Step 4/4: Constructing Multi-Layered Network ---")
  final_multi_layered_network_path <- construct_multi_layered_network_internal( # Call internal version
    gsea_results_file = selected_gsea_file, # Use the selected GSEA file from Step 1
    microbe_pathway_file = selected_microbe_pathway_file, # Use selected MPN file from Step 2
    pathway_jaccard_file = selected_jaccard_file, # Use selected Jaccard file from Step 1
    pathway_metabolite_file = selected_pathway_metabolite_file, # Use selected PMN file from Step 3
    output_directory = mln_output_dir,
    file_type = file_type
  )

  message("\nFull multi-layered network analysis pipeline complete.")
  message("Final integrated multi-layered network saved to: ", final_multi_layered_network_path)

  return(final_multi_layered_network_path)
}
