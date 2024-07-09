#' set colnames and returns the renamed data.frame
#'
#' set colnames
#' @param  object matrix or data frame
#' @param nm names
#'
#' @return the renamed data.frame or matrix
#' @export
setcolnames <- function(object = nm, nm) {
  colnames(object) <- nm
  object
}

#' set rownames and returns the renamed data.frame
#'
#' set colnames
#' @param  object matrix or data frame
#' @param nm names
#'
#' @return the renamed data.frame or matrix
#' @export
setrownames <- function(object = nm, nm) {
  rownames(object) <- nm
  object
}

#' Convert a named list with vector values in each element to a long data.table
#'
#' @param  x a list.
#'
#' @return data.table
#' @importFrom data.table data.table
#' @export
list2dt <- function(x) {
  data.table(name = rep(names(x), sapply(x, length)), value = unlist(x))
}


#' Split big matrix to multiple smaller matrices by column.
#'
#' @param matrix raw matrix.
#' @param chunk_size number of column of each smaller matrix.
#'
#' @return a list containg multiple smaller matrix
#' @export
split_matrix <- function(matrix, chunk_size) {
  chunk_number <- ifelse(ncol(matrix) %% chunk_size == 0,
    ncol(matrix) / chunk_size - 1,
    floor(ncol(matrix) / chunk_size)
  )
  message(sprintf("matrix was divided to %f chunks", chunk_number + 1))
  start_end <- lapply(0:chunk_number, function(x) {
    c(1, chunk_size) + (chunk_size * x)
  })
  start_end[[chunk_number + 1]][[2]] <- ncol(matrix)
  start_end
  matrix_divided <- lapply(start_end, function(x) matrix[, x[[1]]:x[[2]]])

  matrix_divided
}

#' Get attributes.
#'
#' @param x R object with attributes.
#' @param attr_name attributes name.
#'
#' @return attributes name.
#' @export
get_attr <- function(x, attr_name) {
  attributes(x)[[attr_name]]
}
