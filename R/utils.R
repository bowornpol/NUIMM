# R/utils.R

#' Internal helper function to read CSV or TSV files
#'
#' @details
#' Uses high-performance `data.table::fread` to read massive omics datasets
#' significantly faster than base R, returning a standard data.frame.
#'
#' @param file_path A character string specifying the path to the input file.
#' @param file_type A character string indicating the file type ("csv" or "tsv").
#' @param ... Additional arguments to pass to `data.table::fread`.
#' @return A data frame containing the data from the specified file.
#' @keywords internal
read_input_file <- function(file_path, file_type = c("csv", "tsv"), ...) {
  file_type <- match.arg(file_type)

  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The 'data.table' package is required for fast file I/O. Please install it.")
  }

  sep_char <- if (file_type == "csv") "," else "\t"

  # fread automatically handles large files efficiently
  data <- data.table::fread(file_path, sep = sep_char, data.table = FALSE, ...)

  return(data)
}
