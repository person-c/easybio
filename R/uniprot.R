.isJobReady <- function(jobId) {
  url <- paste0("https://rest.uniprot.org/idmapping/status/", jobId)
  status <- httr2::request(url) |>
    httr2::req_retry(max_tries = 10, max_seconds = 60) |>
    httr2::req_perform() |>
    httr2::resp_body_json()

  if (any(!is.null(status[["results"]]), !is.null(status[["failedIds"]]))) {
    message("Job completed!")
    return(TRUE)
  }

  if (!is.null(status[["messages"]])) {
    message(status[["messages"]])
    return(FALSE)
  }
}

.getResultsURL <- function(redirectURL) {
  url <- fifelse(
    redirectURL %flike% "/idmapping/results/",
    gsub("/idmapping/results/", "/idmapping/stream/", redirectURL),
    gsub("/results/", "/results/stream/", redirectURL)
  )

  url
}


#' Map UniProt IDs to Other Identifiers
#'
#' This function maps UniProt IDs to other identifiers using UniProt's ID mapping service.
#' It sends a request to the UniProt API to perform the mapping and retrieves the results in a tabular format.
#'
#' @param ... Parameters to be passed in the request body.
#'
#' @return A `data.table` containing the mapped identifiers.
#' @import data.table
#' @export
#' @examples
#' uniprot_id_map(
#'   ids = "P21802,P12345",
#'   from = "UniProtKB_AC-ID",
#'   to = "UniRef90"
#' )
#'
uniprot_id_map <- function(...) {
  submission <- httr2::request("https://rest.uniprot.org/idmapping/run") |>
    httr2::req_body_form(...) |>
    httr2::req_perform() |>
    httr2::resp_body_json()

  if (.isJobReady(submission[["jobId"]])) {
    url <- paste0("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]])
    details <- httr2::req_perform(httr2::request(url)) |>
      httr2::resp_body_json()
    url <- .getResultsURL(details[["redirectURL"]])
    url <- paste0(url, "?format=tsv")
    text <- httr2::req_perform(httr2::request(url)) |> httr2::resp_body_string()
    resultsTable <- fread(text)

    return(resultsTable)
  }

  stop("Max tries elapsed, error!")
}
