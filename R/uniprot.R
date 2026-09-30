.is_job_ready <- function(job_id, timeout = 60, interval = 2) {
  url <- paste0("https://rest.uniprot.org/idmapping/status/", job_id)
  deadline <- Sys.time() + timeout

  repeat {
    status <- httr2::request(url) |>
      httr2::req_retry(max_tries = 5, max_seconds = 30) |>
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

    # no results yet, so the job is queued or running ("jobStatus" is "NEW" or
    # "RUNNING"): keep asking until the deadline
    if (Sys.time() >= deadline) {
      message("Job is still running after ", timeout, " seconds; giving up")
      return(FALSE)
    }
    Sys.sleep(interval)
  }
}

#' Map UniProt IDs to Other Identifiers
#'
#' This function maps UniProt IDs to other identifiers using UniProt's ID mapping service.
#' It sends a request to the UniProt API to perform the mapping and retrieves the results in a tabular format.
#'
#' @param ... Parameters to be passed in the request body.
#' @param timeout Numeric, how many seconds to wait for the mapping job to
#'   finish. Defaults to 60; large submissions can take considerably longer.
#' @param interval Numeric, how many seconds to wait between status requests.
#'   Defaults to 2.
#'
#' @return A `data.table` containing the mapped identifiers. `NULL`, with a
#'   warning, if the job did not finish within `timeout`.
#' @export
#' @examples
#' \dontrun{
#' uniprot_id_map(
#'   ids = "P21802,P12345",
#'   from = "UniProtKB_AC-ID",
#'   to = "UniRef90"
#' )
#' }
uniprot_id_map <- function(..., timeout = 60, interval = 2) {
  submission <- httr2::request("https://rest.uniprot.org/idmapping/run") |>
    httr2::req_body_form(...) |>
    httr2::req_perform() |>
    httr2::resp_body_json()

  # the API answers with "jobId"; reading a snake_case field silently produced
  # NULL, which sent every follow-up request to an empty job id and made the
  # function report a timeout for perfectly good submissions
  job_id <- submission[["jobId"]]
  if (is.null(job_id)) {
    stop(
      "UniProt did not return a job id; the response had: ",
      paste(names(submission), collapse = ", "),
      call. = FALSE
    )
  }

  if (!.is_job_ready(job_id, timeout = timeout, interval = interval)) {
    warning(
      "The ID mapping job did not finish within ", timeout, " seconds",
      call. = FALSE
    )
    return(NULL)
  }

  # the streaming endpoint takes the job id directly, so the details call (and
  # its "redirectURL" field, also camelCase) is not needed
  url <- paste0("https://rest.uniprot.org/idmapping/stream/", job_id, "?format=tsv")
  text <- httr2::req_perform(httr2::request(url)) |> httr2::resp_body_string()

  fread(text)
}
