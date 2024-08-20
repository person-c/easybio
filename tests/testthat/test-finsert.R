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
test_that("finsert works", {
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
