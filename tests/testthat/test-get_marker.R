test_that("get_marker returns markers for valid human cell types", {
  res <- get_marker(spc = "Human", cell = c("Macrophage", "Monocyte"),
                     number = 5, min_count = 1)
  expect_type(res, "list")
  expect_equal(names(res), c("Macrophage", "Monocyte"))
  expect_true(all(vapply(res, is.character, logical(1))))
  expect_true(all(vapply(res, length, integer(1)) <= 5))
})

test_that("get_marker respects number parameter", {
  res <- get_marker(spc = "Human", cell = "B cell", number = 3)
  expect_true(length(res[["B cell"]]) <= 3)
})

test_that("get_marker suggests corrections for typos", {
  expect_message(
    get_marker(spc = "Human", cell = c("Macrophae", "Monocyte")),
    "did you mean"
  )
})

test_that("get_marker returns only valid cells when mix of valid and invalid", {
  res <- get_marker(spc = "Human", cell = c("Monocyte", "NonExistentCellXYZ"))
  expect_type(res, "list")
  expect_true("Monocyte" %in% names(res))
  expect_false("NonExistentCellXYZ" %in% names(res))
})

test_that("get_marker returns NULL when no valid cells provided", {
  expect_null(suppressMessages(
    get_marker(spc = "Human", cell = c("NotACellTypeXYZ"))
  ))
})

test_that("get_marker works for Mouse species", {
  res <- get_marker(spc = "Mouse", cell = "B cell", number = 3)
  expect_type(res, "list")
  expect_true(length(res[["B cell"]]) > 0)
})

test_that("get_marker returns NULL gracefully for invalid species", {
  res <- suppressMessages(
    get_marker(spc = "Zebrafish", cell = "B cell")
  )
  expect_null(res)
})

test_that("get_marker filters by tissue_class", {
  res <- get_marker(spc = "Human", cell = "B cell",
                     tissue_class = "Blood", number = 5)
  expect_type(res, "list")
})

test_that("get_marker excludes markers below min_count", {
  res_loose <- get_marker(spc = "Human", cell = "B cell",
                           number = 50, min_count = 1)
  res_strict <- get_marker(spc = "Human", cell = "B cell",
                            number = 50, min_count = 1000)

  expect_true(length(res_strict[["B cell"]]) <= length(res_loose[["B cell"]]))
})
