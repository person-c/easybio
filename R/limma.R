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

#' Volcano plot
#'
#' @param x processed DGEList.
#' @param fd_name Fold change symbol.
#' @param fd_hold Folding change threshold.
#' @param ap_name significant name.
#' @param ap_hold significant threshold.
#' @param top number of term to show.
#' @param label aesthetics \code{ggrepel::geom_label_repel}
#' @param legendName new column for color aesthetics.
#' @import data.table
#' @import ggplot2
#' @export
view.volcano <- function(x, fd_name, fd_hold, ap_name, ap_hold, top = FALSE, label, legendName) {
  dt <- setDT(copy(x), keep.rownames = TRUE)
  expr <- substitute(dt[, legendName := fcase(
    ap_name < ap_hold & fd_name > fd_hold, "Up",
    ap_name < ap_hold & fd_name < -fd_hold, "Down",
    default = "Unchanged"
  )])
  eval(expr)

  p <- ggplot(
    data = dt,
    aes(x = {{ fd_name }}, y = -log10({{ ap_name }}), color = {{ legendName }})
  ) +
    geom_point(size = 1.2, alpha = 1, na.rm = TRUE) +
    scale_color_manual(
      values = c(Up = "#414788FF", Unchanged = "darkgray", Down = "#22A884FF"),
      guide = guide_legend(override.aes = list(size = 4))
    ) +
    scale_x_continuous(name = expression(LogFC)) +
    scale_y_continuous(name = substitute(-log10(ap_name))) +
    geom_vline(xintercept = c(-fd_hold, fd_hold), lty = 4, col = "darkgray", lwd = 0.6) +
    geom_hline(yintercept = -log10(fd_hold), lty = 4, col = "darkgray", lwd = 0.6) +
    theme_classic()

  if (top) {
    dttp <- eval(substitute(dt[.(c("Up", "Down")), on = .(legendName)][
      , utils::head(.SD[order(-abs(fd_name))], 5),
      by = .(legendName)
    ]))

    p + ggrepel::geom_label_repel(
      data = dttp, aes(label = {{ label }}),
      label.size = NA, arrow = arrow(length = unit(0.1, "cm")), max.overlaps = 20,
      size = 4, box.padding = unit(0.5, "lines"), min.segment.length = 0,
      point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE
    )
  }
}
