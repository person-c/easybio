expectValue <- c(
  "Uncertain",
  "Uncertain",
  "Monocyte",
  "Uncertain",
  "Macrophage",
  "Monocyte",
  "Uncertain",
  "Uncertain",
  "Macrophage",
  "Uncertain"
)
expectValue <- setNames(expectValue, as.character(0:9))

test_that("finsert works with expression format (legacy)", {
  expect_equal(
    finsert(
      expression(
        c(2, 5) == "Monocyte",
        c(4, 8) == "Macrophage"
      ),
      len = 10,
      na = "Uncertain"
    ),
    expectValue
  )
})

test_that("finsert works with formula list format", {
  res <- finsert(
    list(
      c(2, 5) ~ "Monocyte",
      c(4, 8) ~ "Macrophage"
    ),
    len = 10,
    na = "Uncertain"
  )
  expect_equal(res, expectValue)
})

test_that("finsert auto-extends when len not specified", {
  res <- finsert(
    list(c(0, 3) ~ "A", 5 ~ "B"),
    na = "NA"
  )
  expect_length(res, 6)
  expect_equal(res[["0"]], "A")
  expect_equal(res[["3"]], "A")
  expect_equal(res[["5"]], "B")
  expect_equal(res[["1"]], "NA")
})

test_that("finsert supports setname = FALSE", {
  res <- finsert(
    list(c(0, 1) ~ "Label"),
    len = 3,
    setname = FALSE
  )
  expect_null(names(res))
  expect_equal(res, c("Label", "Label", "Unknown"))
})

test_that("finsert handles overlapping index ranges (last wins)", {
  res <- finsert(
    list(c(0, 1) ~ "First", c(1, 2) ~ "Second"),
    len = 3,
    na = "Other"
  )
  expect_equal(res[["0"]], "First")
  expect_equal(res[["1"]], "Second")
  expect_equal(res[["2"]], "Second")
})

test_that("finsert default parameters work", {
  res <- finsert(list(c(0, 1) ~ "X", 2 ~ "Y"))
  expect_equal(names(res), c("0", "1", "2"))
  expect_equal(res[[1]], "X")
  expect_equal(res[[3]], "Y")
})
