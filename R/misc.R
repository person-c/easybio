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


#' Get graph from a named list according to the overlap between each list element.
#'
#' @param nodes named list.
#'
#' @import data.table
#' @return attributes name.
#' @export
list2graph <- function(nodes) {
  node_x <- c()
  node_y <- c()
  weight <- numeric()
  for (i in seq_along(nodes)) {
    j <- i + 1
    while (j <= length(nodes)) {
      node_x <- append(node_x, names(nodes[i]))
      node_y <- append(node_y, names(nodes[j]))
      weight <- append(weight, intersect(nodes[[i]], nodes[[i]]) |> length())

      j <- j + 1
    }
  }

  net <- data.table(node_x = node_x, node_y = node_y, weight = weight)
  net
}

#' Summary by group according to the regex.
#'
#' @param f function.
#' @param x data
#' @param xname names of x.
#' @param patterns regex to group.
#'
#' @return summary data by group
#' @export
groupStat <- function(f, x, xname = names(x), patterns) {
  i <- lapply(patterns, \(.x) which(xname %like% .x))
  sapply(i, \(.x) force(f)(x[.x]), simplify = FALSE)
}

#' Summary by group according to the index.
#'
#' @param f function.
#' @param x data
#' @param i group index.
#'
#' @return summary data by group
#' @export
groupStatI <- function(f, x, i) {
  sapply(i, \(.x) force(f)(x[.x]), simplify = FALSE)
}
