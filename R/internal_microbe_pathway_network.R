#' Internal Microbe-Pathway Network Helper
#'
#' Calculates relative contribution of taxa to functions.
#'
#' @param path_con_file Character path to processed contribution file.
#' @param metadata_file Character path to metadata file.
#' @param taxonomy_file Character path to taxonomy file.
#' @param output_dir Character path to output directory.
#' @param mpn_filtering Filtering options: "unfiltered", "mean", "median", "top10%", "top25%", "top50%", "top75%".
#' @return Vector of MPN file paths.
#' @keywords internal
con_mpn_int <- function(
    path_con_file,
    metadata_file,
    taxonomy_file,
    output_dir,
    mpn_filtering = "top10%"
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  contrib <- read.csv(path_con_file, stringsAsFactors = FALSE)
  meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
  taxonomy <- read.csv(taxonomy_file, stringsAsFactors = FALSE)

  merged <- merge(contrib, meta, by = "SampleID")
  merged <- merge(merged, taxonomy, by = "FeatureID")
  classes <- unique(merged$class)
  output_paths <- c()

  for (cls in classes) {
    sub_df <- merged[merged$class == cls, ]
    if(nrow(sub_df) == 0) next

    agg <- aggregate(taxon_function_abun ~ FunctionID + TaxonID, data = sub_df, sum)
    func_totals <- aggregate(taxon_function_abun ~ FunctionID, data = agg, sum)
    colnames(func_totals)[2] <- "total_abun"

    res <- merge(agg, func_totals, by = "FunctionID")
    res$relative_contribution <- ifelse(res$total_abun == 0, 0, res$taxon_function_abun / res$total_abun)

    if (mpn_filtering != "unfiltered") {
      if (mpn_filtering %in% c("mean", "median")) {
        FUN_used <- if (mpn_filtering == "mean") mean else median
        thresh <- aggregate(relative_contribution ~ FunctionID, res, FUN_used)
        colnames(thresh)[2] <- "threshold"
        res <- merge(res, thresh, by="FunctionID")
        res <- res[res$relative_contribution >= res$threshold, ]
      } else if (grepl("top", mpn_filtering)) {
        perc <- as.numeric(gsub("top|%", "", mpn_filtering)) / 100
        res <- res %>% dplyr::group_by(FunctionID) %>% dplyr::arrange(desc(relative_contribution)) %>% dplyr::slice_head(prop = perc) %>% dplyr::ungroup() %>% as.data.frame()
      }
    }

    fname <- file.path(output_dir, paste0("mpn_", cls, ".csv"))
    write.csv(res, fname, row.names = FALSE)
    output_paths <- c(output_paths, fname)
  }
  return(output_paths)
}
