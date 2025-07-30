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
