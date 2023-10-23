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
  colnames(object) <- nm
  object
}
