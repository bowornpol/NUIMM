# R/utils.R

#' Internal helper function to read CSV or TSV files
#'
#' This function provides a flexible way to read data from either
#' a comma-separated values (CSV) file or a tab-separated values (TSV) file.
#' It acts as a wrapper around `read.csv` and `read.delim` from base R.
#'
#' @param file_path A character string specifying the path to the input file.
#' @param file_type A character string indicating the file type. Must be either
#'   "csv" (for comma-separated values) or "tsv" (for tab-separated values).
#' @param ... Additional arguments to pass to `read.csv` or `read.delim`
#'   (e.g., `header`, `row.names`, `stringsAsFactors`).
#' @return A data frame containing the data from the specified file.
#' @keywords internal
read_input_file <- function(file_path, file_type = c("csv", "tsv"), ...) {
  file_type <- match.arg(file_type)

  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }

  if (file_type == "csv") {
    data <- read.csv(file_path, ...)
  } else if (file_type == "tsv") {
    data <- read.delim(file_path, ...)
  } else {
    stop("Invalid file_type. Must be 'csv' or 'tsv'.")
  }
  return(data)
}
