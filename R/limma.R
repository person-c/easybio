#' Construct a DGEList Object
#'
#' This function creates a `DGEList` object from a count matrix, sample
#' information, and feature information. It is designed to facilitate the
#' analysis of differential gene expression using the edgeR package.
#'
#' @param count A numeric matrix where rows represent features (e.g., genes) and
#'   columns represent samples. Row names should correspond to feature identifiers,
#'   and column names should correspond to sample identifiers.
#' @param sample.info A data frame containing information about the samples. The
#'   number of rows should match the number of columns in the `count` matrix.
#' @param feature.info A data frame containing information about the features. The
#'   number of rows should match the number of rows in the `count` matrix.
#'
#' @return A `DGEList` object as defined by the edgeR package, which includes the
#'   count data, sample information, and feature information.
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


#' Filter Low-Expressed Genes and Normalize DGEList Data
#'
#' This function filters out low-expressed genes from a `DGEList` object and
#' normalizes the count data. It also provides diagnostic plots for raw and
#' filtered data.
#'
#' @param x A `DGEList` object containing raw count data.
#' @param group.column The name of the column in `x$samples` that contains the
#'   grouping information for the samples.
#' @param min.count The minimum number of counts required for a gene to be
#'   considered expressed. Genes with counts below this threshold in any group
#'   will be filtered out. Defaults to 10.
#'
#' @return The function returns a `DGEList` object with low-expressed genes
#'   filtered out and normalization factors calculated.
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


#' Fit a Linear Model for RNA-seq data using limma
#'
#' This function fits a linear model to processed `DGEList` data using the
#' `limma` package. It defines contrasts between groups and performs
#' differential expression analysis.
#'
#' @param x A processed `DGEList` object containing normalized count data.
#' @param group.column The name of the column in `x$samples` that contains
#'   the grouping information for the samples.
#'
#' @return An `eBayes` object containing the fitted linear model and
#'   results of the differential expression analysis.
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


#' Plot Volcano Plot for Differentially Expressed Genes
#'
#' This function generates a volcano plot for differentially expressed genes
#' (DEGs) using `ggplot2`. It allows for customization of the plot with
#' different aesthetic parameters.
#'
#' @param data A data frame containing the DEGs result.
#' @param data.text A data frame containing labeled data for text annotation.
#' @param x variable representing the aesthetic for the x-axis.
#' @param y variable representing the aesthetic for the y-axis.
#' @param color variable representing the column name for the color aesthetic.
#' @param label variable representing the column name for the text label aesthetic.
#'
#' @return A `ggplot` object representing the volcano plot.
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
