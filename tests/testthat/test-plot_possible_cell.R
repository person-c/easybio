data("pbmc.markers", package = "easybio")

matched <- match_ref(pbmc.markers, n = 30, spc = "Human")

# a hand-made match_ref() result, so the checks below do not depend on the
# database: A matches 2 markers (one detected, one not), B matches 3 (one)
fake_matched <- data.table::data.table(
  cluster = factor(c(0, 0)),
  cell_name = c("A", "B"),
  uniqueN = c(2L, 3L),
  N = c(2L, 3L),
  ordered_symbol = list(c("X", "Y"), c("X", "Y", "Z")),
  orderN = list(c(1L, 1L), c(1L, 1L, 1L)),
  pct_with = list(c(0.9, 0.1), c(0.9, NA, 0.05))
)

test_that("plot_possible_cell keeps its point plot by default", {
  p <- plot_possible_cell(matched, min_unique_n = 2)

  expect_s3_class(p, "ggplot")
  expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
  expect_true(all(c("N", "uniqueN", "cell_name", "cluster") %in% colnames(p$data)))
})

test_that("plot_possible_cell can size the points by uniqueN", {
  p <- plot_possible_cell(matched, min_unique_n = 2, value = "uniqueN")

  expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
  expect_equal(sort(unique(p$data$uniqueN)), sort(unique(matched[uniqueN >= 2, uniqueN])))
})

test_that("min_unique_n includes candidates matching exactly that many markers", {
  expect_equal(sort(plot_possible_cell(fake_matched, min_unique_n = 2)$data$cell_name), c("A", "B"))
  expect_equal(plot_possible_cell(fake_matched, min_unique_n = 3)$data$cell_name, "B")
})

test_that("value = 'pct' shows the share of markers detected at min_pct", {
  p <- plot_possible_cell(fake_matched, min_unique_n = 2, value = "pct", min_pct = 0.5)

  expect_s3_class(p$layers[[1]]$geom, "GeomTile")
  # A: X is detected in 90% of the cells, Y in 10% -> 1 of 2 markers
  expect_equal(p$data[cell_name == "A", pct_supported], 0.5)
  # B: X detected, Y has no recorded rate (NA counts as not detected), Z not
  expect_equal(p$data[cell_name == "B", pct_supported], 1 / 3)

  # the threshold is inclusive
  p_low <- plot_possible_cell(fake_matched, min_unique_n = 2, value = "pct", min_pct = 0.9)
  expect_equal(p_low$data[cell_name == "A", pct_supported], 0.5)
  p_high <- plot_possible_cell(fake_matched, min_unique_n = 2, value = "pct", min_pct = 0.95)
  expect_equal(p_high$data[cell_name == "A", pct_supported], 0)
})

test_that("value = 'pct' errors when the input carries no detection rate", {
  no_rate <- data.table::data.table(
    cluster = factor(0), cell_name = "A", uniqueN = 2L, N = 2L,
    ordered_symbol = list(c("X", "Y")), orderN = list(c(1L, 1L)),
    pct_with = list(c(NA_real_, NA_real_))
  )
  expect_error(plot_possible_cell(no_rate, value = "pct"), "no detection rate")

  # objects from before pct_with existed have no column at all
  stripped <- matched[, !"pct_with"]
  expect_error(plot_possible_cell(stripped, value = "pct"), "pct_with")
})

test_that("plot_possible_cell validates its arguments", {
  expect_error(plot_possible_cell(matched, min_pct = 2), "min_pct")
  expect_error(plot_possible_cell(matched, min_pct = "0.5"), "min_pct")
  expect_error(plot_possible_cell(matched, value = "nope"), "arg")
})
