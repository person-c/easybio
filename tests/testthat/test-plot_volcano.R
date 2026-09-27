degs <- data.frame(
  gene = c("A", "B", "C", "D"),
  log2FC = c(-2.5, -0.2, 1.8, 3.1),
  p_val_adj = c(0.01, 0.7, 0.001, 1e-05),
  group = c("down", "ns", "up", "up")
)
to_label <- degs[degs$group == "up", ]

test_that("plot_volcano returns the points on their own", {
  p <- plot_volcano(degs, x = log2FC, y = p_val_adj, color = group)

  expect_s3_class(p, "ggplot")
  expect_length(p$layers, 1)
  expect_s3_class(p$layers[[1]]$geom, "GeomPoint")
  expect_equal(nrow(p$data), nrow(degs))
})

test_that("plot_volcano adds the labels of data_text to the same plot", {
  skip_if_not_installed("ggrepel")

  p <- plot_volcano(
    degs,
    data_text = to_label,
    x = log2FC, y = p_val_adj, color = group, label = gene
  )

  expect_length(p$layers, 2)
  expect_s3_class(p$layers[[2]]$geom, "GeomTextRepel")
  # the labels are drawn from data_text, the points from the input
  expect_equal(sort(p$layers[[2]]$data$gene), c("C", "D"))
  expect_equal(nrow(p$data), nrow(degs))
})
