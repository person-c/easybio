.is_job_ready <- function(job_id) {
  url <- paste0("https://rest.uniprot.org/idmapping/status/", job_id)
  status <- tryCatch(
    {
      httr2::request(url) |>
        httr2::req_retry(max_tries = 10, max_seconds = 60) |>
        httr2::req_perform() |>
        httr2::resp_body_json()
    },
    error = \(e) list(messages = "error")
  )

  if (any(!is.null(status[["results"]]), !is.null(status[["failedIds"]]))) {
    message("Job completed!")
    return(TRUE)
  }

  if (!is.null(status[["messages"]])) {
    message(status[["messages"]])
    FALSE
  }
}

.get_results_url <- function(redirect_url) {
  url <- fifelse(
    redirect_url %flike% "/idmapping/results/",
    gsub("/idmapping/results/", "/idmapping/stream/", redirect_url),
    gsub("/results/", "/results/stream/", redirect_url)
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
#' @export
#' @examples
#' \dontrun{
#' uniprot_id_map(
#'   ids = "P21802,P12345",
#'   from = "UniProtKB_AC-ID",
#'   to = "UniRef90"
#' )
#' }
uniprot_id_map <- function(...) {
  submission <- httr2::request("https://rest.uniprot.org/idmapping/run") |>
    httr2::req_body_form(...) |>
    httr2::req_perform() |>
    httr2::resp_body_json()

  if (.is_job_ready(submission[["job_id"]])) {
    url <- paste0("https://rest.uniprot.org/idmapping/details/", submission[["job_id"]])
    details <- httr2::req_perform(httr2::request(url)) |>
      httr2::resp_body_json()
    url <- .get_results_url(details[["redirect_url"]])
    url <- paste0(url, "?format=tsv")
    text <- httr2::req_perform(httr2::request(url)) |> httr2::resp_body_string()
    results_table <- fread(text)

    return(results_table)
  }

  warning("Maximum number of tries has been reached!")
  NULL
}
