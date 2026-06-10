utils::globalVariables(c("coef", "pval", "feature", "ID", "core_enrichment", "gene", "term",
                        "lfc", "taxon", "diff_abn"))

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
#' @param ppn_da_method Differential abundance method: "deseq2", "edger", "maaslin2", "maaslin3", "aldex2", "ancombc", or "wilcoxon".
#' @param ppn_map_database Pathway database: "kegg", "metacyc", or "custom".
#' @param ppn_rank_by Gene ranking metric: "signed_log_pvalue", "log2foldchange", or "pvalue".
#' @param ppn_min_gs_size Minimum number of genes in a pathway to be considered.
#' @param ppn_max_gs_size Maximum number of genes in a pathway to be considered.
#' @param ppn_exponent Weight of each step in enrichment score calculation.
#' @param ppn_eps Boundary for calculating the p-value.
#' @param ppn_n_perm Number of permutations.
#' @param ppn_seed Logical or numeric. If TRUE or numeric, sets seed for reproducibility.
#' @param ppn_padjust_method P-value adjustment method for GSEA.
#' @param ppn_filter_by Significance filter for GSEA: "none", "pvalue", or "padjust".
#' @param ppn_pvalue_cutoff P-value cutoff for GSEA significance.
#' @param ppn_padjust_cutoff Adjusted p-value cutoff for GSEA significance.
#' @param ppn_jaccard_cutoff Minimum Jaccard index to retain edges.
#' @param ppn_interaction_method Method for defining pathway interactions: "gsea_core", "database", "metabolite", or "rel_pathway".
#' @param ppn_compound_map Compound-pathway database for metabolite Jaccard: "kegg", "metacyc", or "custom". Only used when ppn_interaction_method = "metabolite".
#' @param ppn_compound_custom_map Path to custom compound-pathway CSV. Required when ppn_compound_map = "custom".
#' @param comparisons_list Optional list of pairwise group comparisons.
#' @return List with `gsea_paths` and `jaccard_paths` character vectors.
#' @keywords internal
con_ppn_int <- function(
  gene_abun_file, metadata_file, map_file, output_dir,
  ppn_da_method = c("maaslin2", "maaslin3", "deseq2", "edger", "aldex2", "ancombc", "wilcoxon"),
  ppn_map_database = c("kegg", "metacyc", "custom"),
  ppn_rank_by = c("signed_log_pvalue", "log2foldchange", "pvalue"),
  ppn_min_gs_size = 10, ppn_max_gs_size = 500, ppn_exponent = 1,
  ppn_eps = 1e-10, ppn_n_perm = 10000, ppn_seed = FALSE,
  ppn_padjust_method = "fdr", ppn_filter_by = c("none", "pvalue", "padjust"),
  ppn_pvalue_cutoff = 0.05, ppn_padjust_cutoff = 0.05,
  ppn_jaccard_cutoff = 0.2, ppn_interaction_method = c("gsea_core", "database", "metabolite", "rel_pathway"),
  ppn_compound_map = c("kegg", "metacyc", "custom"), ppn_compound_custom_map = NULL,
  comparisons_list = NULL
) {
  ppn_da_method <- match.arg(ppn_da_method)
  ppn_map_database <- match.arg(ppn_map_database)
  ppn_rank_by <- match.arg(ppn_rank_by)
  ppn_interaction_method <- match.arg(ppn_interaction_method)
  ppn_filter_by <- match.arg(ppn_filter_by)
  ppn_compound_map <- match.arg(ppn_compound_map)

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
    map_path <- system.file("extdata", "KEGG_pwy_to_KO.csv", package = "NUIMM")
  } else if (ppn_map_database == "metacyc") {
    map_path <- system.file("extdata", "metacyc_pwy_to_rxn_pro_cleaned.csv", package = "NUIMM")
  } else if (ppn_map_database == "custom") {
    map_path <- map_file
  }

  # Fallback to custom mapping if internal database is missing
  if (is.null(map_path) || map_path == "") {
    if (!is.null(map_file) && file.exists(map_file)) {
      message("Warning: Primary database unavailable. Using custom map_file.")
      map_path <- map_file
    } else {
      stop(paste("Mapping file for database", ppn_map_database, "not found. Please use 'custom' and upload a map file."))
    }
  }

  map_raw <- read_input_file(map_path, file_type = "csv", header = FALSE, stringsAsFactors = FALSE)

  if (ppn_map_database == "custom") {
    if (ncol(map_raw) <= 2) {
      stop("Custom mapping file must be in wide format (first column as Pathway ID, second column as Pathway Name, subsequent columns as Genes/KOs).")
    }
    long_df <- tidyr::pivot_longer(map_raw, cols = -c(1, 2), names_to = "drop", values_to = "gene")
    TERM2GENE <- data.frame(term = long_df[[1]], gene = long_df$gene, stringsAsFactors = FALSE)
    TERM2GENE <- TERM2GENE[!is.na(TERM2GENE$gene) & TERM2GENE$gene != "", ]
  } else {
    if (ncol(map_raw) > 2) {
      long_df <- tidyr::pivot_longer(map_raw, cols = -1, names_to = "drop", values_to = "gene")
      TERM2GENE <- data.frame(term = long_df[[1]], gene = long_df$gene, stringsAsFactors = FALSE)
      TERM2GENE <- TERM2GENE[!is.na(TERM2GENE$gene) & TERM2GENE$gene != "", ]
    } else {
      TERM2GENE <- data.frame(term = map_raw[, 1], gene = map_raw[, 2], stringsAsFactors = FALSE)
    }
  }

  term_groups <- NULL
  if (ppn_interaction_method == "database") term_groups <- split(TERM2GENE$gene, TERM2GENE$term)

  # Metabolite-based Jaccard: Load compound-pathway map
  compound_groups <- NULL
  if (ppn_interaction_method == "metabolite") {
    cpd_map_path <- NULL
    if (ppn_compound_map == "kegg") {
      cpd_map_path <- system.file("extdata", "kegg_compound_pathway_map.csv", package = "NUIMM")
      if (cpd_map_path == "") stop("Built-in KEGG compound-pathway mapping not found.")
    } else if (ppn_compound_map == "metacyc") {
      cpd_map_path <- system.file("extdata", "metacyc_compound_pathway_map.csv", package = "NUIMM")
      if (cpd_map_path == "") stop("Built-in MetaCyc compound-pathway mapping not found.")
    } else if (ppn_compound_map == "custom") {
      if (is.null(ppn_compound_custom_map) || !file.exists(ppn_compound_custom_map)) {
        stop("Custom compound mapping file must be provided and exist when ppn_compound_map = 'custom'.")
      }
      cpd_map_path <- ppn_compound_custom_map
    }
    cpd_df <- read.csv(cpd_map_path, stringsAsFactors = FALSE)
    # Build list: pathway -> vector of compound names
    compound_groups <- split(cpd_df$CompoundName, cpd_df$PathwayID)
  }
  
  # Knowledge-based connection: Load hierarchy map
  knowledge_map <- NULL
  if (ppn_interaction_method == "rel_pathway") {
    know_path <- system.file("extdata", "kegg_pathway_knowledge_map.csv", package = "NUIMM")
    if (know_path == "" || !file.exists(know_path)) stop("Built-in KEGG knowledge map not found. Did you generate it?")
    knowledge_map <- read.csv(know_path, stringsAsFactors = FALSE)
  }

  gsea_paths <- c()
  jaccard_paths <- c()

  for (comp in comparisons) {
    cond1 <- comp[1]
    cond2 <- comp[2]
    comp_name <- paste0(cond1, "_vs_", cond2)
    message(sprintf("  Analyzing Pathway-Pathway layer (Comparison: %s).", comp_name))

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
        input_data = t(sub_abun), input_metadata = sub_meta, output = maaslin_out, # Maaslin2 requires transposed abundance
        fixed_effects = "class", reference = c("class", cond1)
      )
      res_df <- fit_data$results
      res_df <- dplyr::rename(res_df, log2FoldChange = coef, pvalue = pval, gene = feature)
      
      # Restore original gene names mapped by Maaslin2
      orig_names <- rownames(sub_abun)
      names_map <- setNames(orig_names, make.names(orig_names))
      res_df$gene <- names_map[res_df$gene]
    } else if (ppn_da_method == "maaslin3") {
      if (!requireNamespace("maaslin3", quietly = TRUE)) stop("Install 'maaslin3' from Bioconductor.")
      maaslin3_out <- file.path(output_dir, paste0("maaslin3_results_", comp_name))
      fit_data <- maaslin3::maaslin3(
        input_data = t(sub_abun), input_metadata = sub_meta, output = maaslin3_out,
        formula = "~ class", reference = paste0("class,", cond1)
      )
      res_file <- file.path(maaslin3_out, "all_results.tsv")
      if (!file.exists(res_file)) stop("MaAsLin3 failed to generate results.")
      res_df <- read.delim(res_file, stringsAsFactors = FALSE)
      
      # Filter for abundance model results
      if ("model" %in% colnames(res_df)) {
        res_df <- res_df[res_df$model == "abundance", ]
      }
      res_df <- dplyr::rename(res_df, log2FoldChange = coef, pvalue = pval_individual, gene = feature)
      
      # Restore original gene names
      orig_names <- rownames(sub_abun)
      names_map <- setNames(orig_names, make.names(orig_names))
      res_df$gene <- names_map[res_df$gene]
    } else if (ppn_da_method == "aldex2") {
      if (!requireNamespace("ALDEx2", quietly = TRUE)) stop("Install 'ALDEx2' from Bioconductor: BiocManager::install('ALDEx2')")
      sub_abun_int <- round(sub_abun)
      conditions_vec <- as.character(sub_meta$class)
      aldex_res <- ALDEx2::aldex(sub_abun_int, conditions_vec, mc.samples = 128, test = "t", effect = TRUE, verbose = FALSE)
      res_df <- data.frame(
        gene = rownames(aldex_res),
        log2FoldChange = aldex_res$effect,
        pvalue = aldex_res$wi.ep,  # Raw p-value for GSEA correction
        stringsAsFactors = FALSE
      )
    } else if (ppn_da_method == "ancombc") {
      if (!requireNamespace("ANCOMBC", quietly = TRUE)) stop("Install 'ANCOMBC' from Bioconductor: BiocManager::install('ANCOMBC')")
      if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) stop("Install 'TreeSummarizedExperiment' from Bioconductor.")
      sub_abun_int <- round(sub_abun)
      sub_abun_int[sub_abun_int < 0] <- 0
      tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = list(counts = as.matrix(sub_abun_int)),
        colData = sub_meta
      )
      ancom_out <- ANCOMBC::ancombc2(
        data = tse, assay_name = "counts", fix_formula = "class",
        p_adj_method = "fdr", verbose = FALSE
      )
      ancom_res <- ancom_out$res
      # Extract results for the second level of the factor
      lfc_col <- grep("^lfc_class", colnames(ancom_res), value = TRUE)[1]
      p_col <- grep("^p_class", colnames(ancom_res), value = TRUE)[1]  # Raw p-value for GSEA correction
      res_df <- data.frame(
        gene = ancom_res$taxon,
        log2FoldChange = ancom_res[[lfc_col]],
        pvalue = ancom_res[[p_col]],
        stringsAsFactors = FALSE
      )
    } else if (ppn_da_method == "wilcoxon") {
      # Wilcoxon Mann-Whitney test
      grp1_vals <- sub_abun[, sub_meta$class == cond1, drop = FALSE]
      grp2_vals <- sub_abun[, sub_meta$class == cond2, drop = FALSE]
      pvals <- apply(sub_abun, 1, function(x) {
        g1 <- x[sub_meta$class == cond1]
        g2 <- x[sub_meta$class == cond2]
        tryCatch(wilcox.test(g1, g2, exact = FALSE)$p.value, error = function(e) NA)
      })
      # Using raw p-values for GSEA correction
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

    final_pvalue_cutoff <- 1.0
    final_p_adjust_method <- ppn_padjust_method

    if (ppn_filter_by == "pvalue") {
      final_pvalue_cutoff <- ppn_pvalue_cutoff
      final_p_adjust_method <- "none"
    } else if (ppn_filter_by == "padjust") {
      final_pvalue_cutoff <- ppn_padjust_cutoff
    }

    gsea_res <- tryCatch(
      {
        clusterProfiler::GSEA(
          geneList = gene_list, TERM2GENE = TERM2GENE, 
          minGSSize = ppn_min_gs_size, maxGSSize = ppn_max_gs_size,
          eps = ppn_eps, pvalueCutoff = final_pvalue_cutoff, 
          pAdjustMethod = final_p_adjust_method, 
          exponent = ppn_exponent, nPermSimple = ppn_n_perm,
          seed = ppn_seed, verbose = FALSE
        )
      },
      error = function(e) {
        message("Error during GSEA computation: ", e$message)
        NULL
      }
    )

    if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
      gsea_out <- as.data.frame(gsea_res)
      message(sprintf("    GSEA yielded %d significant pathways.", nrow(gsea_out)))
      fpath <- file.path(output_dir, paste0("gsea_results_", comp_name, ".csv"))
      write.csv(gsea_out, fpath, row.names = FALSE)
      gsea_paths <- c(gsea_paths, fpath)

      sig_paths <- gsea_out$ID
      jaccard_res <- data.frame()

      if (length(sig_paths) > 1) {
        combos <- combn(sig_paths, 2, simplify = FALSE)
        if (ppn_interaction_method == "gsea_core") {
          genes_list <- strsplit(gsea_out$core_enrichment, "/")
          names(genes_list) <- sig_paths
          for (cb in combos) {
            u <- length(union(genes_list[[cb[1]]], genes_list[[cb[2]]]))
            idx <- ifelse(u > 0, length(intersect(genes_list[[cb[1]]], genes_list[[cb[2]]])) / u, 0)
            if (idx >= ppn_jaccard_cutoff) jaccard_res <- rbind(jaccard_res, data.frame(FunctionID_1 = cb[1], FunctionID_2 = cb[2], jaccard_index = idx))
          }
        } else if (ppn_interaction_method == "database") {
          for (cb in combos) {
            g1 <- if (cb[1] %in% names(term_groups)) term_groups[[cb[1]]] else character(0)
            g2 <- if (cb[2] %in% names(term_groups)) term_groups[[cb[2]]] else character(0)
            u <- length(union(g1, g2))
            idx <- ifelse(u > 0, length(intersect(g1, g2)) / u, 0)
            if (idx >= ppn_jaccard_cutoff) jaccard_res <- rbind(jaccard_res, data.frame(FunctionID_1 = cb[1], FunctionID_2 = cb[2], jaccard_index = idx))
          }
        } else if (ppn_interaction_method == "metabolite") {
          for (cb in combos) {
            m1 <- if (cb[1] %in% names(compound_groups)) compound_groups[[cb[1]]] else character(0)
            m2 <- if (cb[2] %in% names(compound_groups)) compound_groups[[cb[2]]] else character(0)
            u <- length(union(m1, m2))
            idx <- ifelse(u > 0, length(intersect(m1, m2)) / u, 0)
            if (idx >= ppn_jaccard_cutoff) jaccard_res <- rbind(jaccard_res, data.frame(FunctionID_1 = cb[1], FunctionID_2 = cb[2], jaccard_index = idx))
          }
        } else if (ppn_interaction_method == "rel_pathway") {
          for (cb in combos) {
            # Check if this exact pair exists in the knowledge map
            pair_exists <- any(
              (knowledge_map$Pathway1 == cb[1] & knowledge_map$Pathway2 == cb[2]) |
              (knowledge_map$Pathway1 == cb[2] & knowledge_map$Pathway2 == cb[1])
            )
            if (pair_exists && 1 >= ppn_jaccard_cutoff) {
              jaccard_res <- rbind(jaccard_res, data.frame(FunctionID_1 = cb[1], FunctionID_2 = cb[2], jaccard_index = 1))
            }
          }
        }
      }
      message(sprintf("    Pathway-Pathway network generated: %d Jaccard edges.", nrow(jaccard_res)))
      jpath <- file.path(output_dir, paste0("pathway_jaccard_", comp_name, ".csv"))
      write.csv(jaccard_res, jpath, row.names = FALSE)
      jaccard_paths <- c(jaccard_paths, jpath)
    }
  }
  return(list(gsea_paths = gsea_paths, jaccard_paths = jaccard_paths))
}
