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

#' Insert character value in specified index  for character in a weird way.
#'
#' @param x expression
#' @param na the default value for the NA.
#'
#' @return characte
#' @export
finsert <- function(x = expression(c(1, 3) == "a", c(2, 4) == "b"), na = "Unknown") {
  x <- eval(substitute(x))
  x <- lapply(x, as.list)
  x <- rapply(x, eval, classes = "call", how = "replace")
  x <- unlist(x, recursive = FALSE)
  itor <- 1
  v <- character()
  while (itor < length(x)) {
    v[x[[itor + 1]]] <- x[[itor + 2]]
    itor <- itor + 3
  }
  v[is.na(v)] <- na
  v
}
