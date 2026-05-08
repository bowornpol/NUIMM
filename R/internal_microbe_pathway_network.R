#' Internal Microbe-Pathway Network Helper
#'
#' @details
#' Connects Microbes to Pathways utilizing data.table for rapid aggregation
#' and threshold filtering.
#'
#' @keywords internal
utils::globalVariables(c("relative_contribution", "FunctionID", "taxon_function_abun", "total_abun", "TaxonID"))

con_mpn_int <- function(
  path_con_file, metadata_file, taxonomy_file = NULL, output_dir,
  mpn_filtering = "top10%"
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Install 'data.table'.")

  contrib <- read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE)
  meta <- read_input_file(metadata_file, file_type = "csv", stringsAsFactors = FALSE)

  merged <- merge(contrib, meta, by = "SampleID")

  if (!is.null(taxonomy_file)) {
    taxonomy <- read_input_file(taxonomy_file, file_type = "csv", stringsAsFactors = FALSE)
    merged <- merge(merged, taxonomy, by = "FeatureID")
  }

  classes <- unique(merged$class)
  output_paths <- c()

  for (cls in classes) {
    sub_df <- merged[merged$class == cls, ]
    if (nrow(sub_df) == 0) next

    # Fast Aggregation using data.table
    dt <- data.table::as.data.table(sub_df)
    res <- dt[, .(taxon_function_abun = sum(taxon_function_abun)), by = .(FunctionID, TaxonID)]
    res[, total_abun := sum(taxon_function_abun), by = FunctionID]
    res[, relative_contribution := data.table::fifelse(total_abun == 0, 0, taxon_function_abun / total_abun)]

    if (mpn_filtering != "unfiltered") {
      if (mpn_filtering %in% c("mean", "median")) {
        FUN_used <- if (mpn_filtering == "mean") mean else median
        res[, threshold := FUN_used(relative_contribution), by = FunctionID]
        res <- res[relative_contribution >= threshold]
        res[, threshold := NULL]
      } else if (grepl("top", mpn_filtering)) {
        perc <- as.numeric(gsub("top|%", "", mpn_filtering)) / 100
        # Fast top N% slice per group
        res <- res[order(FunctionID, -relative_contribution)]
        res <- res[, .SD[1:max(1, round(.N * perc))], by = FunctionID]
      }
    }

    res_df <- as.data.frame(res)
    if (nrow(res_df) > 0) {
      fname <- file.path(output_dir, paste0("mpn_", cls, ".csv"))
      write.csv(res_df, fname, row.names = FALSE)
      output_paths <- c(output_paths, fname)
    }
  }
  return(output_paths)
}
