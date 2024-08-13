#' Rename Column Names of a Data Frame or Matrix
#'
#' This function renames the column names of a data frame or matrix to the
#' specified names.
#'
#' @param object A data frame or matrix whose column names will be renamed.
#' @param nm A character vector containing the new names for the columns.
#'
#' @return A data frame or matrix with the new column names.
#' @export
setcolnames <- function(object = nm, nm) {
  colnames(object) <- nm
  object
}

#' Rename Row Names of a Data Frame or Matrix
#'
#' This function renames the row names of a data frame or matrix to the
#' specified names.
#'
#' @param object A data frame or matrix whose row names will be renamed.
#' @param nm A character vector containing the new names for the rows.
#'
#' @return A data frame or matrix with the new row names.
#' @export
setrownames <- function(object = nm, nm) {
  rownames(object) <- nm
  object
}

#' Convert a List with Vector Values to a Long Data.table
#'
#' This function converts a named list with vector values in each element to a
#' long data.table. The list is first flattened into a single vector, and then
#' the data.table is created with two columns: one for the name of the original
#' list element and another for the value.
#'
#' @param x A named list where each element contains a vector of values.
#'
#' @return A long data.table with two columns: 'name' and 'value'.
#' @importFrom data.table data.table
#' @export
list2dt <- function(x) {
  data.table(name = rep(names(x), sapply(x, length)), value = unlist(x))
}


#' Split a Matrix into Smaller Submatrices by Column
#'
#' This function splits a matrix into multiple smaller matrices by column.
#' It is useful for processing large matrices in chunks, such as when performing
#' analysis on a single computer with limited memory.
#'
#' @param matrix A numeric or logical matrix to be split.
#' @param chunk_size The number of columns to include in each smaller matrix.
#'
#' @return A list of smaller matrices, each with `chunk_size` columns.
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

#' Retrieve Attributes from an R Object
#'
#' This function extracts a specified attribute from an R object.
#'
#' @param x An R object that has attributes.
#' @param attr_name The name of the attribute to retrieve.
#'
#' @return The value of the attribute with the given name.
#' @export
get_attr <- function(x, attr_name) {
  attributes(x)[[attr_name]]
}


#' Convert a Named List into a Graph Based on Overlap
#'
#' This function creates a graph from a named list, where the edges are determined
#' by the overlap between the elements of the list. Each node in the graph represents
#' an element of the list, and the weight of the edge between two nodes is the number
#' of overlapping elements between the two corresponding lists.
#'
#' @param nodes A named list where each element is a vector.
#'
#' @return A data.table representing the graph, with columns for the node names
#'   (`node_x` and `node_y`) and the weight of the edge (`weight`).
#' @import data.table
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


#' Perform Summary Analysis by Group Using an Index
#'
#' This function applies a specified function to each group defined by an index,
#' and returns a summary of the results. It is useful for summarizing data by
#' group when the groups are defined by an index rather than a named column.
#'
#' @param f A function that takes a single argument and returns a summary of the data.
#' @param x A data frame or matrix containing the data to be summarized.
#' @param idx A vector of indices or group names that define the groups.
#'
#' @return A data frame or matrix containing the summary statistics for each group.
#' @export
groupStatI <- function(f, x, idx) {
  sapply(idx, \(.x) force(f)(x[.x]), simplify = FALSE)
}

#' Perform Summary Analysis by Group Using Regular Expressions
#'
#' This function applies a specified function to each group defined by a regular expression
#' pattern applied to the names of a data object. It is useful for summarizing data when
#' groups are defined by a pattern in the names rather than a specific column or index.
#'
#' @param f A function that takes a single argument and returns a summary of the data.
#' @param x A data frame or matrix containing the data to be summarized.
#' @param xname A character vector containing the names of the variables in `x`.
#' @param patterns A character vector of regular expressions that define the groups.
#'
#' @return A data frame or matrix containing the summary statistics for each group.
#' @export
groupStat <- function(f, x, xname = names(x), patterns) {
  idx <- lapply(patterns, \(.x) which(xname %like% .x))
  groupStatI(f, x, idx)
}


#' Set a Directory for Saving Files
#'
#' This function sets a directory path for saving files, creating the directory if it
#' does not already exist. The directory path is created with the given arguments, which
#' are passed directly to `file.path()`.
#'
#' @param ... Arguments to be passed to `file.path()` to construct the directory path.
#'
#' @return The path to the newly created or existing directory.
#' @export
setSavedir <- function(...) {
  savedir <- file.path(...)
  if (!dir.exists(savedir)) dir.create(savedir, recursive = TRUE)

  return(savedir)
}

#' Perform Operations in a Specified Directory and Return to the Original Directory
#'
#' This function allows you to perform operations in a specified directory and then
#' return to the original directory. It is useful when you need to work with files or
#' directories that are located in a specific location, but you want to return to the
#' original working directory after the operation is complete.
#'
#' @param dir The directory path in which to operate. If the directory does not exist,
#'   it will be created recursively.
#' @param expr An R expression to be evaluated within the specified directory.
#'
#' @return The result of evaluating the expression within the specified directory.
#' @export
workIn <- function(dir, expr) {
  .tmp <- getwd()
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  setwd(dir)

  res <- eval(substitute(expr))
  setwd(.tmp)

  res
}
