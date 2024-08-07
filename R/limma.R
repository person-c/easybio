#' construct DGEList
#'
#' @param count count matrix with rownames(features in row, sample in column).
#' @param sample.info sample information.
#' @param feature.info features information.
#'
#' @return DGEList
#' @import data.table
#' @export
dgeList <- function(count, sample.info, feature.info) {
  stopifnot(nrow(count) == nrow(feature.info))
  stopifnot(ncol(count) == nrow(sample.info))
  sample.info <- sample.info[order(sample.info[[1]], decreasing = TRUE), , drop = FALSE]
  count <- count[rownames(feature.info), rownames(sample.info)]
  x <- edgeR::DGEList(count)
  x$genes <- feature.info
  x$samples <- cbind(x$samples, sample.info)

  x
}


#' filter low expressed genes and normalize data
#'
#' @param x DGEList.
#' @param min.count see \code{edgeR::filterByExpr}
#' @param group.column group column name.
#'
#' @return DGEList
#' @import data.table
#' @import grDevices
#' @import graphics
#' @import stats
#' @import utils
#' @export
dprocess.DGEList <- function(x, group.column, min.count = 10) {
  lcpm <- edgeR::cpm(x, log = TRUE, prior.count = 2)
  # filter low expressed genes
  keep_exprs <- edgeR::filterByExpr(x, group = x$samples[[group.column]], min.count = min.count)
  x <- x[keep_exprs, , keep.lib.sizes = FALSE]

  nsamples <- ncol(x)
  if (nsamples > 10) nsamples <- sample(ncol(x), 10)

  par(mfrow = c(1, 2))
  plot(density(lcpm[, nsamples[[1]]]), lwd = 2, ylim = c(0, 1), las = 2, main = "", xlab = "")
  title(main = "A. Raw data", xlab = "Log-cpm")
  for (i in nsamples) {
    den <- density(lcpm[, i])
    lines(den$x, den$y, lwd = 2, col = sample(colors(), 1))
  }
  lcpm <- edgeR::cpm(x, log = TRUE)
  plot(density(lcpm[, nsamples[[1]]]), lwd = 2, ylim = c(0, 1), las = 2, main = "", xlab = "")
  title(main = "B. Filtered data", xlab = "Log-cpm")
  for (i in nsamples) {
    den <- density(lcpm[, i])
    lines(den$x, den$y, lwd = 2, col = sample(colors(), 1))
  }


  # Normalize the data
  x <- edgeR::calcNormFactors(x)
  lcpm <- edgeR::cpm(x, log = TRUE)
  boxplot(lcpm, las = 2)
  title("Normalized data")
  limma::plotMDS(lcpm,
    label = x$samples[[group.column]],
    col = rleid(x$samples[[group.column]]), dim = c(1, 2)
  )
  title("MDS")

  x
}


#' limma fit
#'
#' @param x processed DGEList.
#' @param group.column group column name.
#' @import data.table
#' @importFrom limma makeContrasts
#' @export
limmaFit <- function(x, group.column) {
  design <- model.matrix(~ 0 + x$samples[[group.column]])
  colnames(design) <- gsub(".*\\]\\]", "", colnames(design))

  all_vs <- utils::combn(unique(x$samples[[group.column]]), 2, simplify = TRUE)
  all_vs2 <- str2expression(paste0(all_vs[1, ], "-", all_vs[2, ]))
  all_vs2 <- setNames(as.list(all_vs2), paste0(all_vs[1, ], "vs", all_vs[2, ]))
  contr_matrix <- do.call("makeContrasts", c(all_vs2, levels = list(colnames(design))))

  v <- limma::voom(x, design, plot = TRUE)
  vfit <- limma::lmFit(v, design)
  vfit <- limma::contrasts.fit(vfit, contrasts = contr_matrix)
  efit <- limma::eBayes(vfit)
  limma::plotSA(efit, main = "Final model: Mean-variance trend")

  setattr(efit, name = "design", design)
  setattr(efit, name = "contrast", contr_matrix)

  efit
}


#' Plot volcano
#'
#' @param data DEGs result.
#' @param data.text labeled data.
#' @param x axis x.
#' @param y axis y.
#' @param color color aesthetic.
#' @param label text label aesthetic.
#' @import ggplot2
#' @export
plotVolcano <- function(data, data.text, x, y, color, label) {
  p <- ggplot(data, aes(x = {{ x }}, y = {{ y }})) +
    geom_point(aes(color = {{ color }})) +
    theme_classic() +
    scale_color_manual(
      values = c("#414788FF", "darkgray", "#22A884FF"),
      guide = guide_legend(override.aes = list(size = 4))
    )

  if (!missing(data.text)) {
    p + ggrepel::geom_text_repel(
      data = data.text, aes(label = {{ label }}),
      arrow = arrow(length = unit(0.1, "cm")),
      max.overlaps = 20,
      size = 4,
      box.padding = unit(0.5, "lines"),
      min.segment.length = 0,
      point.padding = unit(0.8, "lines"),
      segment.color = "black",
      show.legend = FALSE
    )
  }
}
