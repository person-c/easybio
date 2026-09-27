skip_if_not_installed("limma")
skip_if_not_installed("edgeR")

fit_with_groups <- function(labels, per_group = 3L) {
  set.seed(1)
  n <- length(labels) * per_group
  counts <- matrix(
    rpois(200 * n, lambda = 50),
    nrow = 200, ncol = n,
    dimnames = list(paste0("gene", 1:200), paste0("s", 1:n))
  )
  sample_info <- data.frame(group = rep(labels, each = per_group), row.names = colnames(counts))
  dge <- dge_list(counts, sample_info, data.frame(row.names = rownames(counts)))

  # limma_fit() draws the voom and mean-variance plots, send them nowhere
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  suppressMessages(limma_fit(dge, "group"))
}

test_that("limma_fit works with syntactically valid group labels", {
  fit <- fit_with_groups(c("Control", "Treated"))

  expect_s4_class(fit, "MArrayLM")
  expect_equal(colnames(fit), "ControlvsTreated")
  expect_equal(rownames(attr(fit, "contrast")), c("Control", "Treated"))
  expect_equal(unname(attr(fit, "contrast")[, 1]), c(1, -1))
  expect_equal(nrow(fit$coefficients), 200)
})

test_that("limma_fit builds contrasts from labels that are not valid R names", {
  # the TCGA sample types, as handed over by prepare_tcga()
  fit <- fit_with_groups(c("Primary Tumor", "Solid Tissue Normal"))

  expect_s4_class(fit, "MArrayLM")
  expect_equal(colnames(fit), "Primary TumorvsSolid Tissue Normal")
  expect_equal(
    rownames(attr(fit, "contrast")),
    c("Primary.Tumor", "Solid.Tissue.Normal")
  )
  expect_equal(unname(attr(fit, "contrast")[, 1]), c(1, -1))
  # contrasts.fit() matches the two by name
  expect_equal(colnames(attr(fit, "design")), rownames(attr(fit, "contrast")))

  # "Non-tumor - Tumor" used to be parsed as a subtraction of three symbols
  fit_hyphen <- fit_with_groups(c("Non-tumor", "Tumor"))
  expect_equal(colnames(fit_hyphen), "Non-tumorvsTumor")
  expect_equal(unname(attr(fit_hyphen, "contrast")[, 1]), c(1, -1))
})

test_that("limma_fit reports every pairwise contrast", {
  fit <- fit_with_groups(c("A", "B", "C"), per_group = 2L)

  expect_setequal(colnames(fit), c("AvsB", "AvsC", "BvsC"))

  contrast_matrix <- attr(fit, "contrast")
  expect_equal(ncol(contrast_matrix), 3)
  expect_equal(unname(contrast_matrix[c("A", "B"), "AvsB"]), c(1, -1))
  expect_equal(unname(contrast_matrix[c("A", "C"), "AvsC"]), c(1, -1))
  expect_equal(unname(contrast_matrix[c("B", "C"), "BvsC"]), c(1, -1))
})
