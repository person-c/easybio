# library(httr2)
# uniprot_id_mapping <- function(ids, from, to, ...) {
#   baseurl <- "https://rest.uniprot.org/idmapping/run"
#   resp1 <- request(baseurl) |>
#     req_body_form(
#       from = from,
#       to = to,
#       ids = ids,
#       ...
#     ) |>
#     req_perform() |>
#     resp_body_json()


#   status_url <- paste0(baseurl, "status", resp1$jobId, collapse = "/")
#   resp2 <- request(status_url) |>
#     req_perform() |>
#     resp_body_json()

#   resp2
# }
# request("https://rest.uniprot.org/idmapping/run") |>
#   req_body_form(
#     from = "UniProtKB_AC-ID",
#     to = "Gene_Name",
#     ids = "P21802,P12345"
#   ) |>
#   req_perform() |>
#   resp_body_json()

# resp <- request("https://rest.uniprot.org/idmapping/status/fd23a245dd8d39454fcecec721629d2dc439eda5") |>
#   req_perform() |>
#   resp_body_json()
