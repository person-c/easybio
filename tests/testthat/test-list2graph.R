test_that("list2graph weights edges by the overlap between elements", {
  nodes <- list(
    a = c("x", "y", "z"),
    b = c("y", "z"),
    c = "w"
  )
  res <- list2graph(nodes)

  expect_s3_class(res, "data.table")
  expect_equal(colnames(res), c("node1", "node2", "interWeight"))
  expect_equal(nrow(res), 3)
  expect_equal(res[node1 == "a" & node2 == "b", interWeight], 2L)
  expect_equal(res[node1 == "a" & node2 == "c", interWeight], 0L)
  expect_equal(res[node1 == "b" & node2 == "c", interWeight], 0L)
})

test_that("list2graph keeps the order of the input names", {
  res <- list2graph(list(first = "x", second = c("x", "y")))
  expect_equal(res[["node1"]], "first")
  expect_equal(res[["node2"]], "second")
  expect_equal(res[["interWeight"]], 1L)
})
