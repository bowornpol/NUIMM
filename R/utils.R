# R/utils.R

#' Internal helper function to read CSV or TSV files
#'
#' @param file_path A character string specifying the path to the input file.
#' @param file_type A character string indicating the file type ("csv" or "tsv").
#' @param ... Additional arguments to pass to the reader.
#' @return A data frame containing the data from the specified file.
#' @keywords internal
read_input_file <- function(file_path, file_type = NULL, ...) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }

  # Capture extra arguments (like row.names and check.names)
  args <- list(...)

  # Extract and remove row.names so fread doesn't crash
  row_col <- args$row.names
  args$row.names <- NULL

  if (requireNamespace("data.table", quietly = TRUE)) {
    # Fast read using data.table (auto-detects separator)
    data <- do.call(data.table::fread, c(list(file = file_path, data.table = FALSE), args))
  } else {
    # Fallback to base R, guess separator by extension
    ext <- tolower(tools::file_ext(file_path))
    sep_char <- if (ext == "tsv" || (!is.null(file_type) && file_type == "tsv")) "\t" else ","
    data <- read.delim(file_path, sep = sep_char, ...)
  }

  # Manually apply row names if the inner functions requested them
  if (!is.null(row_col) && row_col == 1) {
    rownames(data) <- as.character(data[[1]])
    data <- data[, -1, drop = FALSE]
  }

  return(data)
}
