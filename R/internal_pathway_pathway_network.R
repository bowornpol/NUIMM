utils::globalVariables(c("coef", "pval", "feature", "ID", "core_enrichment", "gene", "term"))

#' Internal Pathway-Pathway Network Helper
#'
#' @details
#' Performs differential abundance analysis, runs Gene Set Enrichment Analysis (GSEA),
#' and calculates Jaccard similarity between significant pathways to build
#' the pathway-pathway sub-network layer.
#'
#' @param gene_abun_file Path to gene/KO abundance table.
#' @param metadata_file Path to sample metadata CSV.
#' @param map_file Path to pathway-to-gene mapping file.
#' @param output_dir Path to output directory.
#' @param ppn_da_method Differential abundance method: "deseq2", "edger", "maaslin2", or "simple".
#' @param ppn_map_database Pathway database: "kegg", "metacyc", or "custom".
#' @param ppn_rank_by Gene ranking metric: "signed_log_pvalue", "log2foldchange", or "pvalue".
#' @param ppn_p_adjust_method P-value adjustment method for GSEA.
#' @param ppn_pvalue_cutoff P-value cutoff for GSEA significance.
#' @param ppn_jaccard_cutoff Minimum Jaccard index to retain edges.
#' @param ppn_jaccard_method Jaccard source: "gsea_core" or "map_file".
#' @param comparisons_list Optional list of pairwise group comparisons.
#' @return List with `gsea_paths` and `jaccard_paths` character vectors.
#' @keywords internal
con_ppn_int <- function(
  gene_abun_file, metadata_file, map_file, output_dir,
  ppn_da_method = c("deseq2", "edger", "maaslin2", "simple"),
  ppn_map_database = c("kegg", "metacyc", "custom"),
  ppn_rank_by = c("signed_log_pvalue", "log2foldchange", "pvalue"),
  ppn_p_adjust_method = "fdr", ppn_pvalue_cutoff = 0.05,
  ppn_jaccard_cutoff = 0.2, ppn_jaccard_method = c("gsea_core", "map_file"),
  comparisons_list = NULL
) {
  ppn_da_method <- match.arg(ppn_da_method)
  ppn_map_database <- match.arg(ppn_map_database)
  ppn_rank_by <- match.arg(ppn_rank_by)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  abun <- read_input_file(gene_abun_file, file_type = "csv", row.names = 1, check.names = FALSE)
  meta <- read_input_file(metadata_file, file_type = "csv", stringsAsFactors = FALSE)

  common_samples <- intersect(colnames(abun), meta$SampleID)
  if (length(common_samples) == 0) stop("No matching samples found!")

  abun <- abun[, common_samples, drop = FALSE]
  meta <- meta[meta$SampleID %in% common_samples, , drop = FALSE]
  rownames(meta) <- meta$SampleID

  if (is.null(comparisons_list)) {
    conditions <- sort(unique(meta$class))
    comparisons <- combn(conditions, 2, simplify = FALSE)
  } else {
    comparisons <- comparisons_list
  }

  map_path <- NULL
  if (ppn_map_database == "kegg") {
    map_path <- system.file("extdata", "kegg_mapping.csv", package = "NUIMM")
  } else if (ppn_map_database == "metacyc") {
    map_path <- system.file("extdata", "metacyc_mapping.csv", package = "NUIMM")
  } else if (ppn_map_database == "custom") {
    map_path <- map_file
  }
  
  # Fallback logic if the internal mapping file is missing
  if (is.null(map_path) || map_path == "") {
    if (!is.null(map_file) && file.exists(map_file)) {
      message("Warning: Internal database file for ", ppn_map_database, " not found. Falling back to the uploaded map_file.")
      map_path <- map_file
    } else {
      stop(paste("Mapping file for database", ppn_map_database, "not found. Please use 'custom' and upload a map file."))
    }
  }

  map_raw <- read_input_file(map_path, file_type = "csv", header = FALSE, stringsAsFactors = FALSE)
  TERM2GENE <- data.frame(term = map_raw[, 1], gene = map_raw[, 2])

  term_groups <- NULL
  if (ppn_jaccard_method == "map_file") term_groups <- split(TERM2GENE$gene, TERM2GENE$term)

  gsea_paths <- c()
  jaccard_paths <- c()

  for (comp in comparisons) {
    cond1 <- comp[1]
    cond2 <- comp[2]
    comp_name <- paste0(cond1, "_vs_", cond2)
    message("  Processing Pathway-Pathway Layer for comparison: ", comp_name)

    sub_meta <- meta[meta$class %in% c(cond1, cond2), ]
    sub_abun <- abun[, sub_meta$SampleID, drop = FALSE]
    res_df <- NULL

    if (ppn_da_method == "deseq2") {
      sub_abun_int <- round(sub_abun)
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = sub_abun_int, colData = sub_meta, design = ~class)
      dds <- DESeq2::DESeq(dds)
      res <- DESeq2::results(dds, contrast = c("class", cond2, cond1))
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
    } else if (ppn_da_method == "edger") {
      y <- edgeR::DGEList(counts = sub_abun, group = sub_meta$class)
      y <- edgeR::calcNormFactors(y)
      y <- edgeR::estimateDisp(y)
      et <- edgeR::exactTest(y, pair = c(cond1, cond2))
      res_df <- et$table
      res_df$log2FoldChange <- res_df$logFC
      res_df$pvalue <- res_df$PValue
      res_df$gene <- rownames(res_df)
    } else if (ppn_da_method == "maaslin2") {
      maaslin_out <- file.path(output_dir, paste0("maaslin_results_", comp_name))
      fit_data <- Maaslin2::Maaslin2(
        input_data = t(sub_abun), input_metadata = sub_meta, output = maaslin_out, # Fixed: Maaslin2 requires transposed abun
        fixed_effects = "class", reference = c("class", cond1)
      )
      res_df <- fit_data$results
      res_df <- dplyr::rename(res_df, log2FoldChange = coef, pvalue = pval, gene = feature)
    } else if (ppn_da_method == "simple") {
      grp1_vals <- sub_abun[, sub_meta$class == cond1, drop = FALSE]
      grp2_vals <- sub_abun[, sub_meta$class == cond2, drop = FALSE]

      if (!requireNamespace("matrixTests", quietly = TRUE)) stop("Install 'matrixTests' for fast Wilcoxon testing.")

      # Fast vectorized Wilcoxon
      wt_res <- matrixTests::row_wilcoxon_twosample(grp1_vals, grp2_vals)
      pvals <- wt_res$pvalue

      fc <- (rowMeans(grp2_vals) + 1e-6) / (rowMeans(grp1_vals) + 1e-6)
      res_df <- data.frame(gene = rownames(sub_abun), log2FoldChange = log2(fc), pvalue = pvals)
    }

    res_df <- res_df[!is.na(res_df$pvalue) & !is.na(res_df$log2FoldChange), ]

    if (ppn_rank_by == "signed_log_pvalue") {
      res_df$rank <- sign(res_df$log2FoldChange) * -log10(res_df$pvalue + 1e-300)
    } else if (ppn_rank_by == "log2foldchange") {
      res_df$rank <- res_df$log2FoldChange
    } else {
      res_df$rank <- -log10(res_df$pvalue + 1e-300)
    }

    res_df <- res_df[order(res_df$rank, decreasing = TRUE), ]
    gene_list <- setNames(res_df$rank, res_df$gene)

    gsea_res <- tryCatch(
      {
        clusterProfiler::GSEA(gene_list, TERM2GENE = TERM2GENE, pvalueCutoff = ppn_pvalue_cutoff, pAdjustMethod = ppn_p_adjust_method, verbose = FALSE)
      },
      error = function(e) {
        NULL
      }
    )

    if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
      gsea_out <- as.data.frame(gsea_res)
      message(sprintf("    Identified %d significant pathways via Gene Set Enrichment Analysis.", nrow(gsea_out)))
      fpath <- file.path(output_dir, paste0("gsea_results_", comp_name, ".csv"))
      write.csv(gsea_out, fpath, row.names = FALSE)
      gsea_paths <- c(gsea_paths, fpath)

      sig_paths <- gsea_out$ID
      jaccard_res <- data.frame()

      if (length(sig_paths) > 1) {
        combos <- combn(sig_paths, 2, simplify = FALSE)
        if (ppn_jaccard_method == "gsea_core") {
          genes_list <- strsplit(gsea_out$core_enrichment, "/")
          names(genes_list) <- sig_paths
          for (cb in combos) {
            u <- length(union(genes_list[[cb[1]]], genes_list[[cb[2]]]))
            idx <- ifelse(u > 0, length(intersect(genes_list[[cb[1]]], genes_list[[cb[2]]])) / u, 0)
            if (idx >= ppn_jaccard_cutoff) jaccard_res <- rbind(jaccard_res, data.frame(FunctionID_1 = cb[1], FunctionID_2 = cb[2], jaccard_index = idx))
          }
        } else if (ppn_jaccard_method == "map_file") {
          for (cb in combos) {
            g1 <- if (cb[1] %in% names(term_groups)) term_groups[[cb[1]]] else character(0)
            g2 <- if (cb[2] %in% names(term_groups)) term_groups[[cb[2]]] else character(0)
            u <- length(union(g1, g2))
            idx <- ifelse(u > 0, length(intersect(g1, g2)) / u, 0)
            if (idx >= ppn_jaccard_cutoff) jaccard_res <- rbind(jaccard_res, data.frame(FunctionID_1 = cb[1], FunctionID_2 = cb[2], jaccard_index = idx))
          }
        }
      }
      message(sprintf("    Constructed pathway-pathway network with %d edges based on Jaccard similarity.", nrow(jaccard_res)))
      jpath <- file.path(output_dir, paste0("pathway_jaccard_", comp_name, ".csv"))
      write.csv(jaccard_res, jpath, row.names = FALSE)
      jaccard_paths <- c(jaccard_paths, jpath)
    }
  }
  return(list(gsea_paths = gsea_paths, jaccard_paths = jaccard_paths))
}
