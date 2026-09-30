data("pbmc.markers", package = "easybio")

test_that("match_ref returns expected structure with default params", {
  res <- match_ref(pbmc.markers, n = 30, spc = "Human")

  expect_s3_class(res, "data.table")
  expect_true(all(c("cluster", "cell_name", "uniqueN", "N", "ordered_symbol", "orderN") %in% colnames(res)))
  expect_true(res[, is.ordered(cluster) || is.factor(cluster)])
})

test_that("match_ref ranks candidates by uniqueN descending, then N, within each cluster", {
  res <- match_ref(pbmc.markers, n = 30, spc = "Human")

  for (cl in unique(res$cluster)) {
    cluster_rows <- res[cluster == cl]
    if (nrow(cluster_rows) > 1) {
      expect_true(all(diff(cluster_rows$uniqueN) <= 0))
      # ties on uniqueN are broken by N
      expect_true(all(cluster_rows[, all(diff(N) <= 0), by = uniqueN][["V1"]]))
    }
  }
})

test_that("match_ref reports pct.1 aligned with ordered_symbol", {
  res <- match_ref(pbmc.markers, n = 30, spc = "Human")

  row <- res[cluster == 0][1]
  expected_pct <- pbmc.markers$pct.1[match(
    row$ordered_symbol[[1]],
    pbmc.markers$gene[pbmc.markers$cluster == 0]
  )]

  expect_length(row$pct_with[[1]], length(row$ordered_symbol[[1]]))
  expect_equal(row$pct_with[[1]], expected_pct)
})

test_that("match_ref reports NA pct_with when the input has no pct.1", {
  no_pct <- pbmc.markers[, c("cluster", "gene", "avg_log2FC", "p_val_adj")]
  res <- match_ref(no_pct, n = 30, spc = "Human")

  expect_true(all(vapply(res$pct_with, function(x) all(is.na(x)), logical(1))))
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
    marker = c("CD3D", "CD3E", "MS4A1", "LYZ"),
    stringsAsFactors = FALSE
  )

  res <- match_ref(pbmc.markers, n = 50, ref = custom_ref)
  expect_s3_class(res, "data.table")
  expect_true(attr(res, "is_custom_ref"))
  expect_equal(nrow(attr(res, "ref")), nrow(custom_ref))
})

test_that("match_ref filters by avg_log2FC threshold", {
  res_strict <- match_ref(pbmc.markers,
    n = 50, spc = "Human",
    avg_log2fc_threshold = 0.5
  )
  res_loose <- match_ref(pbmc.markers,
    n = 50, spc = "Human",
    avg_log2fc_threshold = 0
  )

  expect_true(nrow(res_strict) > 0)
})

test_that("match_ref filters by p_val_adj threshold", {
  res_strict <- match_ref(pbmc.markers,
    n = 50, spc = "Human",
    p_val_adj_threshold = 0.01
  )
  res_loose <- match_ref(pbmc.markers,
    n = 50, spc = "Human",
    p_val_adj_threshold = 1
  )

  expect_true(nrow(res_strict) > 0)
})

test_that("match_ref drops markers below min_pct", {
  loose <- match_ref(pbmc.markers, n = 50, spc = "Human")
  strict <- match_ref(pbmc.markers, n = 50, spc = "Human", min_pct = 0.5)

  # every marker kept as evidence is detected in at least half of the cluster
  expect_true(all(vapply(strict$pct_with, \(x) all(is.na(x) | x >= 0.5), logical(1))))

  # n is applied after the gate, so the markers used are not a subset: genes
  # ranked below the cut move into the top n once the weak ones are gone
  expect_false(identical(
    sort(unique(unlist(loose$ordered_symbol))),
    sort(unique(unlist(strict$ordered_symbol)))
  ))

  # min_pct = 0 keeps everything
  expect_equal(match_ref(pbmc.markers, n = 50, spc = "Human", min_pct = 0)$uniqueN, loose$uniqueN)
})

test_that("match_ref records min_pct for provenance and validates it", {
  res <- match_ref(pbmc.markers, n = 10, spc = "Human", min_pct = 0.25)
  expect_equal(attr(res, "filter_args")$marker_filter[["min_pct"]], 0.25)

  # not applied -> nothing recorded (c() drops the NULL, so the entry is absent)
  plain <- match_ref(pbmc.markers, n = 10, spc = "Human")
  expect_false("min_pct" %in% names(attr(plain, "filter_args")$marker_filter))

  expect_error(match_ref(pbmc.markers, n = 10, spc = "Human", min_pct = 2), "min_pct")
  expect_error(match_ref(pbmc.markers, n = 10, spc = "Human", min_pct = "0.5"), "min_pct")
})

test_that("match_ref ignores min_pct with a message when pct.1 is missing", {
  no_pct <- pbmc.markers[, c("cluster", "gene", "avg_log2FC", "p_val_adj")]

  expect_message(res <- match_ref(no_pct, n = 10, spc = "Human", min_pct = 0.5), "min_pct")
  expect_equal(nrow(res), nrow(match_ref(no_pct, n = 10, spc = "Human")))
})

test_that("match_ref returns empty dt when no markers pass filter", {
  res <- match_ref(pbmc.markers,
    n = 10, spc = "Human",
    avg_log2fc_threshold = 100, p_val_adj_threshold = 1e-300
  )
  expect_equal(nrow(res), 0)
})

test_that("match_ref respects tissue_class filter", {
  res_all <- match_ref(pbmc.markers, n = 30, spc = "Human")
  res_blood <- match_ref(pbmc.markers,
    n = 30, spc = "Human",
    tissue_class = "Blood"
  )

  expect_true(nrow(res_blood) <= nrow(res_all))
})

test_that("match_ref warns when the tissue filters select no reference entry", {
  # "Bone marrow" is a tissue_class, not a tissue_type recorded under Blood,
  # and the two filters are ANDed, so this plausible pair matches nothing
  expect_warning(
    res <- match_ref(pbmc.markers,
      n = 10, spc = "Human",
      tissue_class = "Blood", tissue_type = "Bone marrow"
    ),
    "No reference entry has both"
  )
  expect_equal(nrow(res), 0)
})

test_that("match_ref summarises a long tissue filter instead of printing it", {
  # tissue_type defaults to every type of the species, which must not be
  # spelled out in the message
  expect_warning(
    match_ref(pbmc.markers, n = 10, spc = "Human", tissue_class = "NoSuchTissue"),
    "[0-9]+ values"
  )
})

test_that("match_ref does not warn about the tissue filters when they select entries", {
  expect_no_warning(
    match_ref(pbmc.markers, n = 10, spc = "Human", tissue_class = c("Blood", "Bone marrow"))
  )
  # a custom reference does not go through the tissue filters at all
  expect_no_warning(
    match_ref(pbmc.markers, n = 10, ref = data.frame(cell_name = "T-cell", marker = "CD3D"))
  )
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
  expect_no_error(check_marker(res, cl = 0, top_cell_n = 1))
})

test_that("matchCellMarker2 is deprecated but still works", {
  expect_warning(
    res <- matchCellMarker2(pbmc.markers, n = 10, spc = "Human"),
    "deprecated"
  )
  expect_s3_class(res, "cellmarker_match")
})
