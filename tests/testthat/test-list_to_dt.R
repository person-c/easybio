test_that("list_to_dt converts a named list to a data.table", {
  x <- list(a = c(1, 2), b = c(3, 4, 5))
  res <- list_to_dt(x)
  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), 5)
  expect_equal(colnames(res), c("name", "value"))
  expect_equal(res$name, c("a", "a", "b", "b", "b"))
  expect_equal(res$value, c(1, 2, 3, 4, 5))
})

test_that("list_to_dt accepts custom column names", {
  x <- list(a = "foo", b = c("bar", "baz"))
  res <- list_to_dt(x, col_names = c("group", "item"))
  expect_equal(colnames(res), c("group", "item"))
})

test_that("list_to_dt handles single-element list", {
  x <- list(only = c(1, 2, 3))
  res <- list_to_dt(x)
  expect_equal(nrow(res), 3)
  expect_equal(unique(res$name), "only")
})

test_that("list_to_dt handles empty elements in list", {
  x <- list(a = c(1, 2), b = integer(0), c = 3)
  res <- list_to_dt(x)
  expect_equal(nrow(res), 3)
  expect_false("b" %in% res$name)
})

test_that("list_to_dt handles character vectors", {
  x <- list(g1 = c("a", "b"), g2 = "c")
  res <- list_to_dt(x)
  expect_equal(res$value, c("a", "b", "c"))
})
