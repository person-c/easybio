data("pbmc.markers", package = "easybio")

test_that("deprecated utils aliases warn and work", {
  expect_warning(setcolnames(mtcars, paste0("c", seq_len(ncol(mtcars)))), "deprecated")
  expect_warning(setrownames(mtcars, rownames(mtcars)), "deprecated")
  expect_warning(list2dt(list(a = 1, b = 2)), "deprecated")
  expect_warning(
    list2graph(list(a = c("x", "y"), b = c("y", "z"))),
    "deprecated"
  )
  expect_warning(groupStatI(f = \(x) x, x = mtcars, idx = list(1, 2)), "deprecated")
  expect_warning(groupStat(f = \(x) x, x = mtcars, patterns = list("mp")), "deprecated")
  expect_warning(setSavedir(tempfile("easybio-dir-")), "deprecated")
  expect_warning(workIn(tempfile("easybio-dir-"), 1 + 1), "deprecated")
})

test_that("deprecated plotting aliases warn and work", {
  expect_warning(plotRank(setNames(rnorm(5), letters[1:5])), "deprecated")

  ora_data <- data.frame(
    gene = c("A", "B"),
    p = c(0.01, 0.02),
    count = c(5, 3),
    group = c("x", "x")
  )
  expect_warning(
    plotORA(ora_data, x = count, y = gene, size = p, fill = group),
    "deprecated"
  )

  expect_warning(plotMarkerDistribution("CD14"), "deprecated")

  matched <- match_ref(pbmc.markers, n = 30, spc = "Human")
  expect_warning(plotPossibleCell(matched), "deprecated")
})
