skip_if_not_installed("Seurat")
skip_if_not_installed("patchwork")

# Seurat emits its own ggplot2 deprecation warnings under current versions;
# they are muffled here so a test failure means something
set.seed(1)
genes <- paste0("gene", 1:20)
counts <- matrix(
  abs(rnorm(length(genes) * 30, mean = 2)),
  nrow = length(genes),
  dimnames = list(genes, paste0("cell", 1:30))
)
srt <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
srt$cluster <- rep(c("0", "1"), each = 15)
suppressWarnings(Seurat::Idents(srt) <- "cluster")

features <- list(A = c("gene1", "gene2"), B = c("gene3", "gene4"))

test_that("plot_seurat_dot draws one dot plot for all features", {
  p <- suppressWarnings(plot_seurat_dot(features, srt))

  expect_s3_class(p, "ggplot")
})

test_that("plot_seurat_dot draws one plot per feature when split", {
  p <- suppressWarnings(plot_seurat_dot(features, srt, split = TRUE))

  expect_s3_class(p, "plot_layout")
})

test_that("plot_seurat_dot drops duplicated markers from a single plot", {
  duplicated_features <- list(A = c("gene1", "gene2"), B = c("gene1", "gene3"))

  # Seurat::DotPlot() cannot handle the duplicates itself, it fails with
  # "factor level [3] is duplicated"
  suppressWarnings(expect_warning(
    p <- plot_seurat_dot(duplicated_features, srt),
    "Duplicated markers"
  ))
  expect_s3_class(p, "ggplot")
})
