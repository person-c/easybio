#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param data dataframe.
#' @param matrix return matrix?
#' @param row.name keep rownames?
#' @param keep.factor keep factor?
#'
#' @return return a tidy dataframe or matrix
#' @examples
#' # ADD_EXAMPLES_HERE
tidy <- function(data, matrix = FALSE, row.name = TRUE,
  keep.factor = TRUE) {
  # return a numeric data frame with rownames
  if (is.character(data[[1]])) {
        # first column as rownames
        row_index <- which(duplicated(data[[1]]))
        data <- data[row_index, ]
        rownames(data) <- data[[1]]
        data <- data[, -1]

        # only perserve numeric column
        num_col <- purrr::map_lgl(data, is.numeric)
        return(data[, num_col])
      }

}
