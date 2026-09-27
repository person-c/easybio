skip_if_not_installed("edgeR")
skip_if_not_installed("limma")

run_process <- function() {
  suppressMessages(process_dge_list(dge, "group"))
}

set.seed(1)
genes <- paste0("gene", 1:500)
samples <- paste0("s", 1:12) # more than 10, so the sample subset is exercised
counts <- matrix(
  rpois(length(genes) * length(samples), lambda = 80),
  nrow = length(genes),
  dimnames = list(genes, samples)
)
# half of the genes are barely expressed in every sample, so filterByExpr()
# has something to drop
counts[1:250, ] <- rpois(250 * length(samples), lambda = 0.5)
sample_info <- data.frame(
  group = rep(c("Control", "Treated"), each = 6),
  row.names = samples
)
dge <- dge_list(counts, sample_info, data.frame(row.names = genes))

test_that("process_dge_list returns a normalized DGEList", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  res <- run_process()

  # DGEList is an S4 class since edgeR 4
  expect_s4_class(res, "DGEList")
  expect_false(is.null(res$samples$norm.factors))
  expect_lt(nrow(res), length(genes)) # low-expressed genes were filtered out
})

test_that("process_dge_list does not draw from the random number stream", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  set.seed(1)
  before <- runif(1)
  set.seed(1)
  invisible(run_process())
  after <- runif(1)

  # the sample subset and the line colours used to be sampled, which made the
  # plots differ between runs
  expect_equal(after, before)
})
