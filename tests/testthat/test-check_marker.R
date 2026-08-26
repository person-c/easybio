data("pbmc.markers", package = "easybio")

matched <- match_ref(pbmc.markers, n = 50, spc = "Human")

test_that("check_marker trans mode returns canonical markers", {
  res <- check_marker(matched, cl = 0, top_cell_n = 1, cis = FALSE)
  expect_type(res, "list")
  expect_true(length(res) > 0)
  expect_true(all(vapply(res, is.character, logical(1))))
})

test_that("check_marker cis mode returns local markers from input", {
  res <- check_marker(matched, cl = 0, top_cell_n = 1, cis = TRUE)
  expect_type(res, "list")
  expect_true(length(res) > 0)
  expect_true(all(vapply(res, is.character, logical(1))))
})

test_that("check_marker returns results for multiple clusters", {
  res <- check_marker(matched, cl = c(0, 1), top_cell_n = 2)
  expect_type(res, "list")
  expect_true(length(res) >= 2)
})

test_that("check_marker returns results for top_cell_n > 1", {
  res <- check_marker(matched, cl = 0, top_cell_n = 3)
  expect_type(res, "list")
  expect_true(length(res) >= 2)
})

test_that("check_marker cis mode shows matched genes that led to annotation", {
  res_cis <- check_marker(matched, cl = 0, top_cell_n = 1, cis = TRUE)
  expect_type(res_cis, "list")
  expect_true(length(res_cis) >= 1)

  marker_genes <- unique(pbmc.markers$gene)
  for (cell_genes in res_cis) {
    expect_true(all(cell_genes %in% marker_genes))
  }
})

test_that("check_marker errors when spc is missing from match_ref result", {
  matched_no_spc <- match_ref(pbmc.markers, n = 10, ref = data.frame(
    cell_name = "A", marker = "RPS12", stringsAsFactors = FALSE
  ))
  expect_error(
    check_marker(matched_no_spc, cl = 0, cis = FALSE),
    "species"
  )
})

test_that("check_marker rejects non-cellmarker_match input", {
  plain_dt <- data.table::data.table(
    cluster = factor(0), cell_name = "A", uniqueN = 1,
    N = 1, ordered_symbol = list("X"), orderN = list(1)
  )
  expect_error(
    check_marker(plain_dt, cl = 0),
    "match_ref"
  )
})
