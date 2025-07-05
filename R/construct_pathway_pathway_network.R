#' Pathway-pathway network construction
#'
#' This function performs differential expression analysis (using DESeq2) and
#' Gene Set Enrichment Analysis (GSEA) to identify enriched pathways and
#' subsequently computes Jaccard indices to quantify overlap between core
#' enrichment genes of significant pathways. It allows for flexible p-value
#' adjustment and gene ranking methods.
#'
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
#' @param abundance_file A character string specifying the path to the CSV file
#'   containing gene abundance data. The first column should be gene IDs,
#'   and subsequent columns should be sample counts (integers).
#' @param metadata_file A character string specifying the path to the CSV file
#'   containing sample metadata. Must include 'SampleID' and 'class' columns.
#' @param map_file A character string specifying the path to the CSV file
#'   containing pathway-to-gene mapping. Expected to be a two-column CSV
#'   where the first column is Pathway ID and the second is Gene ID.
#' @param output_file A character string specifying the path to the directory
#'   where output CSV files (GSEA results and Jaccard indices) will be saved.
#'   The directory will be created if it does not exist.
#' @param pvalueCutoff A numeric value specifying the adjusted p-value cutoff
#'   for determining significance in GSEA results.
#' @param pAdjustMethod A character string specifying the method for p-value
#'   adjustment. Must be one of "fdr" (False Discovery Rate), "holm", "hochberg",
#'   "hommel", "bonferroni", "BH", "BY", or "none".
#' @param rank_by A character string specifying the method to rank genes for GSEA.
#'   Must be either "signed_log_pvalue" (sign of log2FoldChange * -log10(p-value))
#'   or "log2FoldChange".
#' @return The function's primary output is CSV files saved
#'   to the specified `output_file` directory, containing GSEA results and
#'   pathway Jaccard indices for each comparison.
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom clusterProfiler GSEA
#' @importFrom dplyr %>% filter select all_of distinct group_by arrange mutate row_number n ungroup
#' @importFrom tidyr gather
#' @export
construct_pathway_pathway_network <- function(
  abundance_file,
  metadata_file,
  map_file,
  output_file,
  pvalueCutoff,
  pAdjustMethod = c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none"),
  rank_by = c("signed_log_pvalue", "log2FoldChange")
) {
  # Validate inputs
  pAdjustMethod <- match.arg(pAdjustMethod)
  rank_by <- match.arg(rank_by) # Validate new parameter

  message("Starting pathway-pathway network construction.")
  message("Using p-value cutoff: ", pvalueCutoff, " and p-adjustment method: ", pAdjustMethod, ".")
  message("Genes will be ranked by: ", rank_by, ".") # Status message for ranking method

  # Create output directory if it doesn't exist
  if (!dir.exists(output_file)) {
    dir.create(output_file, recursive = TRUE)
    message("Created output directory: ", output_file)
  } else {
    message("Output directory already exists: ", output_file)
  }

  # 1. Load abundance data from provided file
  message("1. Loading abundance data from: ", abundance_file)
  gene_abundance <- tryCatch(
    read.csv(abundance_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading abundance file '", abundance_file, "': ", e$message, sep = ""))
    }
  )
  message("    Loaded abundance data. Dimensions: ", paste(dim(gene_abundance), collapse = "x"))

  # 2. Load sample metadata
  message("2. Loading sample metadata from: ", metadata_file)
  metadata <- tryCatch(
    read.csv(metadata_file, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) {
      stop(paste("Error loading metadata file '", metadata_file, "': ", e$message, sep = ""))
    }
  )
  message("    Loaded metadata. Dimensions: ", paste(dim(metadata), collapse = "x"))

  # Ensure the correct SampleID and class columns are present
  if (!"SampleID" %in% colnames(metadata)) {
    stop("The 'sample_metadata.csv' file must contain a 'SampleID' column.")
  }
  if (!"class" %in% colnames(metadata)) {
    stop("The 'sample_metadata.csv' file must contain a 'class' column for group definition.")
  }

  # 3. Filter metadata samples present in abundance data
  message("3. Filtering metadata to match samples in abundance data...")
  sample_ids <- colnames(gene_abundance)
  initial_metadata_rows <- nrow(metadata)
  metadata <- metadata %>% filter(SampleID %in% sample_ids)
  if (nrow(metadata) == 0) {
    stop("No matching samples found between abundance data and metadata after filtering.")
  }
  message("    Filtered metadata. Kept ", nrow(metadata), " out of ", initial_metadata_rows, " samples.")

  # 4. Set 'condition' factor from 'class' column
  message("4. Setting 'condition' factor from 'class' column and aligning data...")
  metadata$condition <- as.factor(metadata$class)
  rownames(metadata) <- metadata$SampleID # Use SampleID for row names

  # 5. Round abundance counts for DESeq2 compatibility and ensure sample order
  message("5. Rounding abundance counts and aligning sample order for DESeq2...")
  gene_abundance_rounded <- round(gene_abundance)
  gene_abundance_rounded <- gene_abundance_rounded[, rownames(metadata), drop = FALSE] # drop=FALSE to handle single sample case
  message("    Abundance data ready for DESeq2. Final dimensions: ", paste(dim(gene_abundance_rounded), collapse = "x"))

  # 6. Create DESeq2 dataset and run differential expression analysis
  message("6. Creating DESeq2 dataset and running differential expression analysis...")
  dds <- tryCatch(
    DESeqDataSetFromMatrix(countData = gene_abundance_rounded,
                           colData = metadata,
                           design = ~ condition),
    error = function(e) {
      stop(paste("Error creating DESeq2 dataset: ", e$message, ". Check abundance data (must be integers) and metadata consistency.", sep = ""))
    }
  )

  dds <- tryCatch(
    DESeq(dds),
    error = function(e) {
      stop(paste("Error running DESeq2 analysis: ", e$message, ". This might happen if groups have zero variance, or too few samples.", sep = ""))
    }
  )
  message("    DESeq2 analysis complete.")

  # 7. Get all pairwise condition comparisons
  conditions <- levels(metadata$condition)
  if (length(conditions) < 2) {
    stop("Less than two unique conditions found in metadata 'class' column. Cannot perform pairwise comparisons.")
  }
  comparisons <- combn(conditions, 2, simplify = FALSE)
  message("7. Identified ", length(comparisons), " pairwise comparisons: ",
          paste(sapply(comparisons, function(x) paste(x, collapse = " vs ")), collapse = ", "))

  # 8. Load pathway-to-gene mapping and reshape into TERM2GENE format
  message("8. Loading pathway-to-gene mapping from: ", map_file)
  map_raw <- tryCatch(
    read.csv(map_file, header = FALSE, fill = TRUE, stringsAsFactors = FALSE, skip = 1),
    error = function(e) {
      stop(paste("Error loading map file '", map_file, "': ", e$message, ". Ensure it's a valid CSV.", sep = ""))
    }
  )

  TERM2GENE <- map_raw %>%
    gather(key = "temp_col", value = "gene", -V1) %>% # V1 becomes 'term', others become 'gene'
    select(term = V1, gene) %>%
    filter(gene != "") %>% # Remove empty gene entries
    distinct() # Ensure unique pathway-gene pairs

  if (nrow(TERM2GENE) == 0) {
    stop("No valid pathway-gene mappings found after processing '", map_file, "'. Check file format.")
  }
  message("    Processed ", nrow(TERM2GENE), " unique pathway-gene mappings.")

  gsea_results_list <- list()

  # 9. Loop over each pairwise comparison to run GSEA
  message("9. Running GSEA for each pairwise comparison...")
  for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    cond1 <- comp[1]
    cond2 <- comp[2]
    comparison_name <- paste0(cond1, "_vs_", cond2)
    message("    Processing comparison (", i, "/", length(comparisons), "): ", cond2, " vs ", cond1)

    # Get DESeq2 results for contrast cond2 vs cond1
    res <- results(dds, contrast = c("condition", cond2, cond1))

    # Prepare ranked gene list
    ranked_df <- as.data.frame(res[, c("log2FoldChange", "pvalue")])
    # Filter out NA values for relevant columns before ranking
    ranked_df <- ranked_df[!is.na(ranked_df$log2FoldChange) & !is.na(ranked_df$pvalue), ]

    if (nrow(ranked_df) == 0) {
      warning("No valid log2FoldChange or pvalue for comparison ", comparison_name, ". Skipping GSEA.")
      next
    }

    # Apply chosen ranking method
    if (rank_by == "signed_log_pvalue") {
      # Handle cases where pvalue might be 0, leading to -log10(0) = Inf
      min_pvalue_for_log <- min(ranked_df$pvalue[ranked_df$pvalue > 0], na.rm = TRUE) / 2
      ranked_df$pvalue[ranked_df$pvalue == 0] <- min_pvalue_for_log
      ranked_df$rank <- sign(ranked_df$log2FoldChange) * -log10(ranked_df$pvalue)
      message("    Ranking by signed -log10(p-value) for ", comparison_name, ".")
    } else if (rank_by == "log2FoldChange") {
      ranked_df$rank <- ranked_df$log2FoldChange
      message("    Ranking by log2FoldChange for ", comparison_name, ".")
    }

    ranked_df <- ranked_df[order(ranked_df$rank, decreasing = TRUE), ]
    geneList <- setNames(ranked_df$rank, rownames(ranked_df))

    # Run GSEA using clusterProfiler
    gsea_res <- tryCatch(
      GSEA(geneList = geneList,
           TERM2GENE = TERM2GENE,
           pvalueCutoff = pvalueCutoff,
           pAdjustMethod = pAdjustMethod,
           seed = TRUE,
           verbose = FALSE),
      error = function(e) {
        warning(paste("GSEA failed for comparison ", comparison_name, ": ", e$message, ". Skipping.", sep = ""))
        return(NULL)
      }
    )

    if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
      message("    No significant GSEA results found for ", comparison_name, ".")
      next
    }

    # Save GSEA results dataframe
    gsea_df <- as.data.frame(gsea_res)
    key <- comparison_name
    gsea_results_list[[key]] <- gsea_df

    # --- Construct full output path for GSEA results with parameters ---

    # Format parameters for filename
    pvalueCutoff_fname <- as.character(pvalueCutoff) # Convert numeric to string
    pAdjustMethod_fname <- pAdjustMethod
    rank_by_fname <- rank_by

    gsea_output_filename <- paste0(
      "gsea_results_", key, "_",
      pvalueCutoff_fname, "_",
      pAdjustMethod_fname, "_",
      rank_by_fname,
      ".csv"
    )

    gsea_output_path <- file.path(output_file, gsea_output_filename)
    write.csv(gsea_df, gsea_output_path, row.names = FALSE)
    message("    Saved GSEA results for ", comparison_name, " to: ", gsea_output_path)
  }

  # 10. Compute Jaccard indices between pathways within each comparison's GSEA results
  message("10. Computing Jaccard indices for overlapping pathways...")
  jaccard_results_list <- list()

  if (length(gsea_results_list) == 0) {
    message("    No GSEA results to compute Jaccard indices. Skipping.")
  } else {
    for (key in names(gsea_results_list)) {
      message("    Calculating Jaccard index for comparison: ", key)
      gsea_df <- gsea_results_list[[key]]

      # Filter for pathways with core enrichment genes (i.e., not empty or NA)
      gsea_df_filtered <- gsea_df %>% filter(!is.na(core_enrichment) & core_enrichment != "")
      gene_sets <- strsplit(as.character(gsea_df_filtered$core_enrichment), "/")
      gene_sets <- lapply(gene_sets, function(x) unique(na.omit(x))) # Ensure unique genes per set and remove NAs

      n <- length(gene_sets)
      res_list <- list()

      if (n < 2) { # Ensure there are at least two pathways to compare
        message(paste("    Less than two significant pathways with core enrichment for comparison '", key, "'. Skipping Jaccard index calculation for this comparison.", sep = ""))
        next
      }

      for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
          genes_i <- gene_sets[[i]]
          genes_j <- gene_sets[[j]]

          # Skip if either gene set is empty
          if (length(genes_i) == 0 || length(genes_j) == 0) {
            next
          }

          intersection <- length(intersect(genes_i, genes_j))
          union <- length(union(genes_i, genes_j))
          jaccard <- ifelse(union == 0, 0, intersection / union)

          if (jaccard > 0) {
            res_list[[length(res_list) + 1]] <- data.frame(
              FunctionID_1 = gsea_df_filtered$ID[i], # Renamed from pathway_1
              FunctionID_2 = gsea_df_filtered$ID[j], # Renamed from pathway_2
              jaccard_index = jaccard,
              comparison = key,
              stringsAsFactors = FALSE
            )
          }
        }
      }

      if (length(res_list) > 0) {
        jaccard_df <- do.call(rbind, res_list)
        # Construct full output path for Jaccard results
        jaccard_output_path <- file.path(output_file, paste0("pathway_jaccard_", key, ".csv"))
        write.csv(jaccard_df, jaccard_output_path, row.names = FALSE)
        message("    Saved Jaccard index results for ", key, " to: ", jaccard_output_path)
      } else {
        message("    No Jaccard indices > 0 found for comparison '", key, "'. Skipping saving file.")
      }
    }
  }

  message("Pathway-pathway network construction complete.")
  return(invisible(NULL))
}
