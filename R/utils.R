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
setcolnames <- function(object, nm) {
  if (length(nm) != ncol(object)) {
    stop("Length of 'nm' must equal the number of columns of 'object'")
  }
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
setrownames <- function(object, nm) {
  if (length(nm) != nrow(object)) {
    stop("Length of 'nm' must equal the number of rows of 'object'")
  }
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
#' @param col_names The colnames of the returned result.
#'
#' @return A long data.table with two columns: 'name' and 'value'.
#' @export
#' @examples
#' library(easybio)
#' list2dt(list(a = c(1, 1), b = c(2, 2)))
list2dt <- function(x, col_names = c("name", "value")) {
  res <- data.table(name = rep(names(x), sapply(x, length)), value = unlist(x))
  setnames(res, new = col_names)
  res
}

#' Split a Matrix into Smaller Sub-matrices by Column or Row
#'
#' This function splits a matrix into multiple smaller matrices by column or row.
#' It is useful for processing large matrices in chunks, such as when performing
#' analysis on a single computer with limited memory.
#'
#' @param matrix A numeric or logical matrix to be split.
#' @param chunk_size The number of columns or rows to include in each smaller matrix.
#' @param column  Divided by column(default is `TRUE`)
#'
#' @return A list of smaller matrices, each with `chunk_size` columns or rows.
#' @export
#' @examples
#' library(easybio)
#' split_matrix(mtcars, chunk_size = 2)
#' split_matrix(mtcars, chunk_size = 5, column = FALSE)
split_matrix <- function(matrix, chunk_size, column = TRUE) {
  n <- if (column) ncol(matrix) else nrow(matrix)
  starts <- seq(1, n, by = chunk_size)
  ends <- pmin(starts + chunk_size - 1, n)
  num_chunks <- length(starts)
  message(sprintf("Matrix was divided into %d chunks", num_chunks))

  lapply(seq_len(num_chunks), function(i) {
    s <- starts[i]
    e <- ends[i]
    if (column) matrix[, s:e, drop = FALSE] else matrix[s:e, , drop = FALSE]
  })
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
#'   (`node_1` and `node_2`) and the weight of the edge (`interWeight`).
#' @export
list2graph <- function(nodes) {
  comb2 <- combn(names(nodes), m = 2, simplify = FALSE)
  inter <- lapply(comb2, \(x) length(intersect(nodes[[x[[1]]]], nodes[[x[[2]]]])))

  data.table(
    node1 = sapply(comb2, \(x) x[[1]]),
    node2 = sapply(comb2, \(x) x[[2]]),
    interWeight = as.integer(inter)
  )
}


#' Perform Summary Analysis by Group Using an column Index
#'
#' This function applies a specified function to each group defined by an column index,
#' and returns a summary of the results. It is useful for summarizing data by
#' group when the groups are defined by an  column index.
#'
#' @param f A function that takes a single argument and returns a summary of the data.
#' @param x A data frame or matrix containing the data to be summarized.
#' @param idx A list of indices or group names that define the column groups.
#'
#' @return A list containing the summary statistics for each group.
#' @export
#' @examples
#' library(easybio)
#' groupStatI(f = \(x) x + 1, x = mtcars, idx = list(c(1, 10), 2))
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
#' @param patterns A list of regular expressions that define the groups.
#'
#' @return A list containing the summary statistics for each group.
#' @export
#' @examples
#' library(easybio)
#' groupStat(f = \(x) x + 1, x = mtcars, patterns = list("mp", "t"))
groupStat <- function(f, x, xname = colnames(x), patterns) {
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
  oldwd <- getwd()
  on.exit(setwd(oldwd))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  setwd(dir)

  res <- eval(substitute(expr))
  res
}

#' Extract Unique Elements from a Column with Optional Filtering
#'
#' Retrieves the unique, non-missing values from a specified column of a data frame.
#' An optional expression can be provided to filter the rows of the data frame
#' before extracting the values.
#'
#' @param data A data frame from which to extract values.
#' @param col_name A single string specifying the name of the target column.
#' @param subset An optional logical expression used to subset the data frame.
#'   This expression is evaluated in the context of the `data`, so columns can be
#'   referred to by their names directly (e.g., `Sepal.Length > 5`).
#'
#' @return A vector containing the unique, non-NA values from the specified
#'   column after the optional filtering has been applied.
#'
#' @export
#'
#' @examples
#' # Example 1: Get all unique species from the iris dataset
#' available_ele(iris, "Species")
#'
#' # Example 2: Get unique species for flowers with Sepal.Length > 7
#' available_ele(iris, "Species", subset = Sepal.Length > 7)
#'
#' # Example 3: Get unique carb values for cars with 6 cylinders
#' available_ele(mtcars, "carb", subset = cyl == 6)
available_ele <- function(data, col_name, subset) {
  assert_data_frame(data)
  assert_string(col_name)
  assert_subset(col_name, choices = names(data))

  if (!missing(subset)) {
    subset_expr <- substitute(subset)
    row_idx <- eval(subset_expr, envir = data, enclos = parent.frame())

    data <- data[row_idx & !is.na(row_idx), , drop = FALSE]
  }

  values <- data[[col_name]]
  unique(na.omit(values))
}


#' Suggest Best Matches for a String from a Vector of Choices
#'
#' This function provides intelligent suggestions for a user's input string by
#' finding the best matches from a given vector of choices. It follows a
#' multi-layered approach:
#' 1.  Performs normalization (case-insensitivity, trimming whitespace).
#' 2.  Checks for an exact match first for maximum performance and accuracy.
#' 3.  If no exact match, it uses a combination of fuzzy string matching
#'     (Levenshtein distance via `adist`) to catch typos and partial/substring
#'     matching (`grep`) to handle incomplete input.
#' 4.  Ranks the potential matches and returns the best suggestion(s).
#'
#' @param x A single character string; the user input to find matches for.
#' @param choices A character vector of available, valid options.
#' @param n An integer specifying the maximum number of suggestions to return.
#'   Defaults to 1.
#' @param threshold An integer; the maximum Levenshtein distance to consider a
#'   choice a "close" match. A lower value is stricter. Defaults to 2.
#' @param ignore.case A logical value. If `TRUE`, matching is case-insensitive.
#'   Defaults to `TRUE`.
#' @param return_distance A logical value. If `TRUE`, the output is a data.frame
#'   containing the suggestions and their calculated distance/score. Defaults to
#'   `FALSE`.
#'
#' @return
#' By default (`return_distance = FALSE`), returns a character vector of the
#' best `n` suggestions. If no suitable match is found, returns `NA`.
#' If `return_distance = TRUE`, returns a `data.frame` with columns
#' `suggestion` and `distance`, or `NULL` if no match is found.
#'
#' @export
#'
#' @examples
#' # --- Setup ---
#' cell_types <- c(
#'   "B cell", "T cell", "Macrophage", "Monocyte", "Neutrophil",
#'   "Natural Killer T-cell", "Dendritic cell"
#' )
#'
#' # --- Usage ---
#' # 1. Exact match (after normalization)
#' suggest_best_match("t cell", cell_types)
#' #> [1] "T cell"
#'
#' # 2. Typo correction (fuzzy match)
#' suggest_best_match("Macrophaeg", cell_types)
#' #> [1] "Macrophage"
#'
#' # 3. Partial input (substring match)
#' suggest_best_match("Mono", cell_types)
#' #> [1] "Monocyte"
#'
#' # 4. Requesting multiple suggestions
#' suggest_best_match("t", cell_types, n = 3)
#' #> [1] "T cell" "Neutrophil" "Natural Killer T-cell"
#'
#' # 5. No good match found
#' suggest_best_match("Erythrocyte", cell_types)
#' #> [1] NA
#'
#' # 6. Returning suggestions with their distance score
#' suggest_best_match("t ce", cell_types, n = 3, return_distance = TRUE)
#' #>              suggestion distance
#' #> 1                T cell        1
#' #> 2        Dendritic cell        2
#' #> 3 Natural Killer T-cell        2
suggest_best_match <- function(x,
                               choices,
                               n = 1,
                               threshold = 2,
                               ignore.case = TRUE,
                               return_distance = FALSE) {
  # --- 1. Input Validation and Normalization ---
  stopifnot(
    is.character(x), length(x) == 1,
    is.character(choices)
  )

  if (length(choices) == 0) {
    return(if (return_distance) NULL else NA_character_)
  }

  # Normalize input and choices
  input_norm <- if (ignore.case) tolower(trimws(x)) else trimws(x)
  choices_norm <- if (ignore.case) tolower(trimws(choices)) else trimws(choices)

  # --- 2. Exact Match ---
  exact_match_idx <- which(choices_norm == input_norm)
  if (length(exact_match_idx) > 0) {
    if (return_distance) {
      return(data.frame(suggestion = choices[exact_match_idx[1]], distance = 0))
    } else {
      return(choices[exact_match_idx[1]])
    }
  }

  # --- 3. Gather Candidates from Fuzzy and Partial Matching ---
  # Fuzzy matching (Levenshtein distance) for typos
  distances <- adist(input_norm, choices_norm, ignore.case = FALSE)
  fuzzy_idx <- which(distances <= threshold)

  # Partial matching (grep) for substrings
  partial_idx <- grep(input_norm, choices_norm, ignore.case = FALSE)

  # Combine candidates into a data.frame with their scores
  # We give partial matches a low, fixed score (e.g., 0.5) to rank them highly.
  candidates <- rbind(
    if (length(fuzzy_idx) > 0) data.frame(idx = fuzzy_idx, score = distances[fuzzy_idx]),
    if (length(partial_idx) > 0) data.frame(idx = partial_idx, score = 0.5)
  )

  if (is.null(candidates) || nrow(candidates) == 0) {
    return(if (return_distance) NULL else NA_character_)
  }

  # --- 4. Rank and Select Best Matches ---
  # Order by score (lower is better), then remove duplicates, keeping the best score
  candidates <- candidates[order(candidates$score), ]
  best_candidates <- candidates[!duplicated(candidates$idx), ]

  # Get the top N results
  top_n <- head(best_candidates, n)

  if (nrow(top_n) == 0) {
    return(if (return_distance) NULL else NA_character_)
  }

  # --- 5. Format Output ---
  if (return_distance) {
    data.frame(
      suggestion = choices[top_n$idx],
      distance = top_n$score,
      row.names = NULL
    )
  } else {
    choices[top_n$idx]
  }
}
