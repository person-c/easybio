test_that("exact match returns single best result", {
  choices <- c("B cell", "T cell", "Macrophage", "Monocyte", "Neutrophil")
  expect_equal(suggest_best_match("T cell", choices), "T cell")
  expect_equal(suggest_best_match("B cell", choices), "B cell")
})

test_that("case-insensitive matching works", {
  choices <- c("B cell", "T cell", "Macrophage")
  expect_equal(suggest_best_match("t cell", choices), "T cell")
  expect_equal(suggest_best_match("MACROPHAGE", choices), "Macrophage")
})

test_that("whitespace is trimmed during matching", {
  choices <- c("B cell", "T cell")
  expect_equal(suggest_best_match("  T cell  ", choices), "T cell")
})

test_that("fuzzy matching catches typos within threshold", {
  choices <- c("Macrophage", "Monocyte", "Neutrophil")
  expect_equal(suggest_best_match("Macrophaeg", choices), "Macrophage")
  expect_equal(suggest_best_match("Macrophagee", choices), "Macrophage")
})

test_that("fuzzy matching with ignore.case = FALSE", {
  choices <- c("Monocyte", "Macrophage")
  expect_equal(
    suggest_best_match("monocyte", choices, ignore.case = FALSE),
    "Monocyte"
  )
})

test_that("substring matching finds partial inputs", {
  choices <- c("B cell", "T cell", "Macrophage", "Monocyte")
  expect_equal(suggest_best_match("Mono", choices), "Monocyte")
  expect_equal(suggest_best_match("Macro", choices), "Macrophage")
})

test_that("returns multiple suggestions when n > 1", {
  choices <- c("T cell", "Neutrophil", "Natural Killer T-cell", "Dendritic cell")
  res <- suggest_best_match("t", choices, n = 3)
  expect_length(res, 3)
  expect_true("T cell" %in% res)
})

test_that("returns NA when no match found", {
  choices <- c("B cell", "T cell")
  expect_true(is.na(suggest_best_match("Erythrocyte", choices)))
})

test_that("return_distance = TRUE returns data.frame", {
  choices <- c("T cell", "B cell", "Macrophage")
  res <- suggest_best_match("t cell", choices, return_distance = TRUE)
  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("suggestion", "distance"))
  expect_equal(res$suggestion[1], "T cell")
  expect_equal(res$distance[1], 0)
})

test_that("return_distance = TRUE returns NULL when no match", {
  choices <- c("Xylophone", "Querulous")
  res <- suggest_best_match("Zz", choices, return_distance = TRUE)
  expect_null(res)
})

test_that("empty choices returns NA or NULL", {
  expect_true(is.na(suggest_best_match("x", character(0))))
  expect_null(suggest_best_match("x", character(0), return_distance = TRUE))
})

test_that("invalid input types are rejected", {
  expect_error(suggest_best_match(1, c("a", "b")))
  expect_error(suggest_best_match("x", 1:5))
})

test_that("exact match is preferred over fuzzy even when fuzzy exists", {
  choices <- c("Monocyte", "Monocytes")
  expect_equal(suggest_best_match("Monocyte", choices), "Monocyte")
})
