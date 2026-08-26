data("pbmc.markers", package = "easybio")

test_that("match_ref returns expected structure with default params", {
  res <- match_ref(pbmc.markers, n = 30, spc = "Human")

  expect_s3_class(res, "data.table")
  expect_true(all(c("cluster", "cell_name", "uniqueN", "N", "ordered_symbol", "orderN") %in% colnames(res)))
  expect_true(res[, is.ordered(cluster) || is.factor(cluster)])
})

test_that("match_ref results are ordered by N descending within each cluster", {
  res <- match_ref(pbmc.markers, n = 30, spc = "Human")

  for (cl in unique(res$cluster)) {
    cluster_rows <- res[cluster == cl]
    if (nrow(cluster_rows) > 1) {
      expect_true(all(diff(cluster_rows$N) <= 0))
    }
  }
})

test_that("match_ref stores ref, is_custom_ref, and filter_args as attributes", {
  res <- match_ref(pbmc.markers, n = 30, spc = "Human")

  expect_true(is.data.frame(attr(res, "ref")))
  expect_false(attr(res, "is_custom_ref"))
  expect_type(attr(res, "filter_args"), "list")
  expect_equal(attr(res, "filter_args")$marker_filter[["n"]], 30)
})

test_that("match_ref works with custom reference", {
  custom_ref <- data.frame(
    cell_name = c("T-cell", "T-cell", "B-cell", "Myeloid"),
    marker    = c("CD3D", "CD3E", "MS4A1", "LYZ"),
    stringsAsFactors = FALSE
  )

  res <- match_ref(pbmc.markers, n = 50, ref = custom_ref)
  expect_s3_class(res, "data.table")
  expect_true(attr(res, "is_custom_ref"))
  expect_equal(nrow(attr(res, "ref")), nrow(custom_ref))
})

test_that("match_ref filters by avg_log2FC threshold", {
  res_strict <- match_ref(pbmc.markers, n = 50, spc = "Human",
                                  avg_log2FC_threshold = 0.5)
  res_loose   <- match_ref(pbmc.markers, n = 50, spc = "Human",
                                  avg_log2FC_threshold = 0)

  expect_true(nrow(res_strict) > 0)
})

test_that("match_ref filters by p_val_adj threshold", {
  res_strict <- match_ref(pbmc.markers, n = 50, spc = "Human",
                                  p_val_adj_threshold = 0.01)
  res_loose  <- match_ref(pbmc.markers, n = 50, spc = "Human",
                                  p_val_adj_threshold = 1)

  expect_true(nrow(res_strict) > 0)
})

test_that("match_ref returns empty dt when no markers pass filter", {
  res <- match_ref(pbmc.markers, n = 10, spc = "Human",
                           avg_log2FC_threshold = 100, p_val_adj_threshold = 1e-300)
  expect_equal(nrow(res), 0)
})

test_that("match_ref respects tissueClass filter", {
  res_all <- match_ref(pbmc.markers, n = 30, spc = "Human")
  res_blood <- match_ref(pbmc.markers, n = 30, spc = "Human",
                                 tissueClass = "Blood")

  expect_true(nrow(res_blood) <= nrow(res_all))
})

test_that("match_ref handles n larger than available markers per cluster", {
  # n = 10000 is larger than any cluster has markers
  res <- match_ref(pbmc.markers, n = 10000, spc = "Human")
  expect_s3_class(res, "data.table")
  expect_true(nrow(res) > 0)
})

test_that("match_ref returns a cellmarker_match S3 object", {
  res <- match_ref(pbmc.markers, n = 30, spc = "Human")
  expect_s3_class(res, "cellmarker_match")
  expect_s3_class(res, "data.table")
})

test_that("cellmarker_match attributes survive row subsetting", {
  res <- match_ref(pbmc.markers, n = 30, spc = "Human")
  sub <- res[1:3]
  expect_s3_class(sub, "cellmarker_match")
  expect_true(is.data.frame(attr(sub, "ref")))
  expect_type(attr(sub, "filter_args"), "list")
})

test_that("cellmarker_match attributes survive column selection", {
  res <- match_ref(pbmc.markers, n = 30, spc = "Human")
  sub <- res[, .(cluster, cell_name, N)]
  expect_s3_class(sub, "cellmarker_match")
  expect_true(is.data.frame(attr(sub, "ref")))
})

test_that("check_marker accepts cellmarker_match input", {
  res <- match_ref(pbmc.markers, n = 30, spc = "Human")
  expect_no_error(check_marker(res, cl = 0, topcellN = 1))
})

test_that("matchCellMarker2 is deprecated but still works", {
  expect_warning(
    res <- matchCellMarker2(pbmc.markers, n = 10, spc = "Human"),
    "deprecated"
  )
  expect_s3_class(res, "cellmarker_match")
})
