utils::globalVariables(c("relative_contribution", "FunctionID", "taxon_function_abun", "total_abun", "TaxonID"))

#' Internal Microbe-Pathway Network Helper
#' @keywords internal

con_mpn_int <- function(
  path_con_file, metadata_file, taxonomy_file = NULL, output_dir,
  mpn_filtering = "top10%"
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Read Data
  contrib <- read_input_file(path_con_file, file_type = "csv", stringsAsFactors = FALSE)
  meta <- read_input_file(metadata_file, file_type = "csv", stringsAsFactors = FALSE)

  # Merge Metadata
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

    message(sprintf("  Processing Microbe-Pathway Layer for class: %s", cls))
    initial_pairs <- nrow(unique(sub_df[, c("FunctionID", "TaxonID")]))

    # BULLETPROOF MATH using dplyr
    res <- sub_df |>
      dplyr::group_by(FunctionID, TaxonID) |>
      dplyr::summarise(taxon_function_abun = sum(taxon_function_abun), .groups = "drop") |>
      dplyr::group_by(FunctionID) |>
      dplyr::mutate(total_abun = sum(taxon_function_abun)) |>
      dplyr::mutate(relative_contribution = ifelse(total_abun == 0, 0, taxon_function_abun / total_abun)) |>
      dplyr::ungroup() |>
      as.data.frame()

    # Filtering Logic
    if (mpn_filtering != "unfiltered") {
      if (mpn_filtering %in% c("mean", "median")) {
        FUN_used <- if (mpn_filtering == "mean") mean else median

        thresh <- aggregate(relative_contribution ~ FunctionID, res, FUN_used)
        colnames(thresh)[2] <- "threshold"
        res <- merge(res, thresh, by = "FunctionID")
        res <- res[res$relative_contribution >= res$threshold, ]
        res$threshold <- NULL
      } else if (grepl("top", mpn_filtering)) {
        perc <- as.numeric(gsub("top|%", "", mpn_filtering)) / 100
        res <- res |>
          dplyr::group_by(FunctionID) |>
          dplyr::arrange(dplyr::desc(relative_contribution)) |>
          dplyr::slice_head(prop = perc) |>
          dplyr::ungroup() |>
          as.data.frame()
      }
    }

    if (nrow(res) > 0) {
      message(sprintf("    Retained %d microbe-pathway associations involving %d unique microbial taxa.", nrow(res), length(unique(res$TaxonID))))
      fname <- file.path(output_dir, paste0("mpn_", cls, ".csv"))
      write.csv(res, fname, row.names = FALSE)
      output_paths <- c(output_paths, fname)
    }
  }
  return(output_paths)
}
