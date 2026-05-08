# R/utils.R

#' Internal helper function to read CSV or TSV files
#'
#' @details
#' Uses high-performance `data.table::fread` to read massive omics datasets
#' significantly faster than base R, handling row.names safely.
#'
#' @param file_path A character string specifying the path to the input file.
#' @param file_type A character string indicating the file type ("csv" or "tsv").
#' @param ... Additional arguments to pass to the reader.
#' @return A data frame containing the data from the specified file.
#' @keywords internal
read_input_file <- function(file_path, file_type = c("csv", "tsv"), ...) {
  file_type <- match.arg(file_type)

  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }

  # Capture extra arguments (like row.names and check.names)
  args <- list(...)

  # Extract and remove row.names so fread doesn't crash
  row_col <- args$row.names
  args$row.names <- NULL

  if (requireNamespace("data.table", quietly = TRUE)) {
    sep_char <- if (file_type == "csv") "," else "\t"
    # Fast read using data.table
    data <- do.call(data.table::fread, c(list(file = file_path, sep = sep_char, data.table = FALSE), args))
  } else {
    # Fallback to base R
    sep_char <- if (file_type == "csv") "," else "\t"
    data <- read.delim(file_path, sep = sep_char, ...)
  }

  # Manually apply row names if the inner functions requested them
  if (!is.null(row_col) && row_col == 1) {
    rownames(data) <- as.character(data[[1]])
    data <- data[, -1, drop = FALSE]
  }

  return(data)
}
