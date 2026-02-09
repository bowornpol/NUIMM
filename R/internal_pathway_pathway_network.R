#' Internal Pathway-Pathway Network Helper (Modified for Manual Reference)
#'
#' Runs DA, GSEA, and Jaccard index calculation.
#'
#' @param gene_abun_file Character path to gene abundance file.
#' @param metadata_file Character path to metadata file.
#' @param map_file Character path to mapping file.
#' @param output_dir Character path to output directory.
#' @param ppn_da_method DA method options: "deseq2", "edger", "maaslin2", "simple".
#' @param ppn_map_database Database options: "kegg", "metacyc", "custom".
#' @param ppn_rank_by GSEA ranking options: "signed_log_pvalue", "log2foldchange", "pvalue".
#' @param ppn_p_adjust_method GSEA p-adj options: "fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none".
#' @param ppn_pvalue_cutoff Numeric p-value cutoff.
#' @param ppn_jaccard_cutoff Numeric Jaccard cutoff.
#' @param comparisons_list List of character vectors (e.g., list(c("Healthy", "Disease"))). First element is Reference.
#' @return List of GSEA and Jaccard file paths.
#' @keywords internal
con_ppn_int <- function(
    gene_abun_file,
    metadata_file,
    map_file,
    output_dir,
    ppn_da_method = c("deseq2", "edger", "maaslin2", "simple"),
    ppn_map_database = c("kegg", "metacyc", "custom"),
    ppn_rank_by = c("signed_log_pvalue", "log2foldchange", "pvalue"),
    ppn_p_adjust_method = "fdr",
    ppn_pvalue_cutoff = 0.05,
    ppn_jaccard_cutoff = 0.2,
    comparisons_list = NULL
) {
  ppn_da_method <- match.arg(ppn_da_method)
  ppn_map_database <- match.arg(ppn_map_database)
  ppn_rank_by <- match.arg(ppn_rank_by)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 1. Load Data
  abun <- read.csv(gene_abun_file, row.names = 1, check.names = FALSE)
  meta <- read.csv(metadata_file, stringsAsFactors = FALSE)

  # 2. Match Samples
  common_samples <- intersect(colnames(abun), meta$SampleID)
  if(length(common_samples) == 0) stop("No matching samples found! Check SampleIDs.")

  abun <- abun[, common_samples]
  meta <- meta[meta$SampleID %in% common_samples, ]
  rownames(meta) <- meta$SampleID

  # 3. Define Comparisons
  if (is.null(comparisons_list)) {
    # Default: Sort alphabetically if user provided nothing
    conditions <- sort(unique(meta$class))
    comparisons <- combn(conditions, 2, simplify = FALSE)
  } else {
    # Use the manual list provided by the user
    comparisons <- comparisons_list
  }

  gsea_paths <- c()
  jaccard_paths <- c()

  for (comp in comparisons) {
    cond1 <- comp[1] # Reference Group
    cond2 <- comp[2] # Target Group
    comp_name <- paste0(cond1, "_vs_", cond2)
    message("Analyzing: ", comp_name)
    message("   -> Reference Group: ", cond1)
    message("   -> Target Group:    ", cond2)

    sub_meta <- meta[meta$class %in% c(cond1, cond2), ]
    sub_abun <- abun[, sub_meta$SampleID]
    res_df <- NULL

    # --- Differential Abundance ---
    if (ppn_da_method == "deseq2") {
      sub_abun_int <- round(sub_abun)
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = sub_abun_int, colData = sub_meta, design = ~ class)
      dds <- DESeq2::DESeq(dds)
      res <- DESeq2::results(dds, contrast = c("class", cond1, cond2))
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)

    } else if (ppn_da_method == "edger") {
      y <- edgeR::DGEList(counts = sub_abun, group = sub_meta$class)
      y <- edgeR::calcNormFactors(y)
      y <- edgeR::estimateDisp(y)
      et <- edgeR::exactTest(y, pair = c(cond2, cond1))
      res_df <- et$table
      res_df$log2FoldChange <- res_df$logFC
      res_df$pvalue <- res_df$PValue
      res_df$gene <- rownames(res_df)

    } else if (ppn_da_method == "maaslin2") {
      # FIXED: Dynamic folder name to prevent overwriting
      maaslin_out <- file.path(output_dir, paste0("maaslin_results_", comp_name))

      fit_data <- Maaslin2::Maaslin2(
        input_data      = sub_abun,
        input_metadata  = sub_meta,
        output          = maaslin_out,
        fixed_effects   = "class",
        reference       = c("class", cond1) # FIXED: Uses cond1 as Reference
      )
      res_df <- fit_data$results
      res_df <- dplyr::rename(res_df, log2FoldChange = coef, pvalue = pval, gene = feature)

    } else if (ppn_da_method == "simple") {
      grp1_vals <- sub_abun[, sub_meta$class == cond1]
      grp2_vals <- sub_abun[, sub_meta$class == cond2]
      pvals <- apply(sub_abun, 1, function(x) {
        g1 <- x[sub_meta$class == cond1]; g2 <- x[sub_meta$class == cond2]
        wt <- tryCatch(wilcox.test(g1, g2), error = function(e) NA)
        if(is.list(wt)) wt$p.value else NA
      })
      # Calculate Fold Change: Target (cond2) / Reference (cond1)
      fc <- (rowMeans(grp2_vals) + 1e-6) / (rowMeans(grp1_vals) + 1e-6)
      log2fc <- log2(fc)
      res_df <- data.frame(gene = rownames(sub_abun), log2FoldChange = log2fc, pvalue = pvals)
    }

    # --- GSEA Analysis ---
    map_raw <- read.csv(map_file, header = FALSE, stringsAsFactors = FALSE)
    TERM2GENE <- data.frame(term = map_raw[,1], gene = map_raw[,2])
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

    gsea_res <- tryCatch({
      clusterProfiler::GSEA(gene_list, TERM2GENE = TERM2GENE, pvalueCutoff = ppn_pvalue_cutoff, pAdjustMethod = ppn_p_adjust_method, verbose = FALSE)
    }, error = function(e) NULL)

    if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
      gsea_out <- as.data.frame(gsea_res)
      fname <- paste0("gsea_results_", comp_name, ".csv")
      fpath <- file.path(output_dir, fname)
      write.csv(gsea_out, fpath, row.names = FALSE)
      gsea_paths <- c(gsea_paths, fpath)

      # --- Jaccard Index ---
      sig_paths <- gsea_out$ID
      genes_list <- strsplit(gsea_out$core_enrichment, "/")
      names(genes_list) <- sig_paths
      jaccard_res <- data.frame()

      if (length(genes_list) > 1) {
        combos <- combn(names(genes_list), 2, simplify = FALSE)
        for (cb in combos) {
          u <- length(union(genes_list[[cb[1]]], genes_list[[cb[2]]]))
          i <- length(intersect(genes_list[[cb[1]]], genes_list[[cb[2]]]))
          idx <- ifelse(u > 0, i/u, 0)
          if (idx >= ppn_jaccard_cutoff) {
            jaccard_res <- rbind(jaccard_res, data.frame(FunctionID_1 = cb[1], FunctionID_2 = cb[2], jaccard_index = idx))
          }
        }
      }
      jname <- paste0("pathway_jaccard_", comp_name, ".csv")
      jpath <- file.path(output_dir, jname)
      write.csv(jaccard_res, jpath, row.names = FALSE)
      jaccard_paths <- c(jaccard_paths, jpath)
    }
  }
  return(list(gsea_paths = gsea_paths, jaccard_paths = jaccard_paths))
}
