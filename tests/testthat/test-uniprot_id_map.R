# the UniProt API is mocked throughout: these tests pin the response shapes the
# service actually returns (camelCase jobId, results under "results")

json_response <- function(body) {
  httr2::response(
    200L,
    headers = "content-type: application/json",
    body = charToRaw(body)
  )
}

text_response <- function(body) {
  httr2::response(
    200L,
    headers = "content-type: text/plain",
    body = charToRaw(body)
  )
}

test_that("uniprot_id_map() reads the job id the API returns", {
  seen <- character()

  httr2::local_mocked_responses(function(req) {
    seen <<- c(seen, req$url)
    if (grepl("/idmapping/run$", req$url)) {
      return(json_response('{"jobId":"JOB123"}'))
    }
    if (grepl("/idmapping/status/", req$url, fixed = TRUE)) {
      return(json_response('{"results":[{"from":"P21802"}]}'))
    }
    if (grepl("/idmapping/stream/", req$url, fixed = TRUE)) {
      return(text_response("From\tTo\nP21802\tUniRef90_P21802\n"))
    }
    httr2::response(404L)
  })

  res <- uniprot_id_map(ids = "P21802", from = "UniProtKB_AC-ID", to = "UniRef90")

  expect_s3_class(res, "data.table")
  expect_equal(colnames(res), c("From", "To"))
  expect_equal(res[["To"]], "UniRef90_P21802")
  # the id from the run response has to reach the follow-up requests
  expect_true(any(grepl("JOB123", seen)))
})

test_that("uniprot_id_map() keeps asking while the job is still running", {
  status_calls <- 0L

  httr2::local_mocked_responses(function(req) {
    if (grepl("/idmapping/run$", req$url)) {
      return(json_response('{"jobId":"JOB9"}'))
    }
    if (grepl("/idmapping/status/", req$url, fixed = TRUE)) {
      status_calls <<- status_calls + 1L
      if (status_calls < 3L) {
        return(json_response('{"jobStatus":"RUNNING"}'))
      }
      return(json_response('{"results":[{"from":"P21802"}]}'))
    }
    text_response("From\tTo\nP21802\tX\n")
  })

  res <- uniprot_id_map(
    ids = "P21802", from = "UniProtKB_AC-ID", to = "UniRef90",
    interval = 0.01
  )

  expect_equal(status_calls, 3L)
  expect_s3_class(res, "data.table")
})

test_that("uniprot_id_map() gives up when the job never finishes", {
  httr2::local_mocked_responses(function(req) {
    if (grepl("/idmapping/run$", req$url)) {
      return(json_response('{"jobId":"JOB9"}'))
    }
    json_response('{"jobStatus":"RUNNING"}')
  })

  expect_warning(
    res <- uniprot_id_map(
      ids = "P21802", from = "UniProtKB_AC-ID", to = "UniRef90",
      timeout = 0.2, interval = 0.01
    ),
    "did not finish"
  )
  expect_null(res)
})

test_that("uniprot_id_map() reports a response without a job id", {
  httr2::local_mocked_responses(function(req) json_response('{"error":"nope"}'))

  expect_error(
    uniprot_id_map(ids = "P21802", from = "UniProtKB_AC-ID", to = "UniRef90"),
    "job id"
  )
})
