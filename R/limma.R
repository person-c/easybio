#' Construct a DGEList Object
#'
#' This function creates a `DGEList` object from a count matrix, sample
#' information, and feature information. It is designed to facilitate the
#' analysis of differential gene expression using the edgeR package.
#'
#' @param count A numeric matrix where rows represent features (e.g., genes) and
#'   columns represent samples. Row names should correspond to feature identifiers,
#'   and column names should correspond to sample identifiers.
#' @param sample_info A data frame containing information about the samples. The
#'   number of rows should match the number of columns in the `count` matrix.
#' @param feature_info A data frame containing information about the features. The
#'   number of rows should match the number of rows in the `count` matrix.
#'
#' @return A `DGEList` object as defined by the edgeR package, which includes the
#'   count data, sample information, and feature information.
#' @export
dge_list <- function(count, sample_info, feature_info) {
  assert_names(rownames(count), identical.to = rownames(feature_info))
  assert_names(colnames(count), identical.to = rownames(sample_info))

  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop(
      "To construct DGEList object, dge_list() requires 'edgeR' package which ",
      "cannot be found. Please install 'edgeR' using 'BiocManager::install('edgeR')'",
      call. = FALSE
    )
  }
  x <- edgeR::DGEList(count, samples = sample_info, genes = feature_info)

  x
}

#' Filter Low-Expressed Genes and Normalize DGEList Data
#'
#' This function filters out low-expressed genes from a `DGEList` object and
#' normalizes the count data. It also provides diagnostic plots for raw and
#' filtered data.
#'
#' @details
#' At most the first 10 samples are shown in the density plots, and each of them
#' keeps its colour, so repeated runs give the same plots.
#'
#' @param x A `DGEList` object containing raw count data.
#' @param group_column The name of the column in `x$samples` that contains the
#'   grouping information for the samples.
#' @param min_count The minimum number of counts required for a gene to be
#'   considered expressed. Genes with counts below this threshold in any group
#'   will be filtered out. Defaults to 10.
#'
#' @return The function returns a `DGEList` object with low-expressed genes
#'   filtered out and normalization factors calculated.
#' @export
process_dge_list <- function(x, group_column, min_count = 10) {
  lcpm <- edgeR::cpm(x, log = TRUE, prior.count = 2)
  # filter low expressed genes
  keep_exprs <- edgeR::filterByExpr(x, group = x$samples[[group_column]], min.count = min_count)
  x <- x[keep_exprs, , keep.lib.sizes = FALSE]

  # at most 10 samples, taken in order: both the subset and the line colours
  # used to be drawn at random, so the diagnostic plots differed between runs
  samples <- seq_len(min(ncol(x), 10L))
  line_cols <- rainbow(length(samples))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(1, 2))
  plot(density(lcpm[, samples[[1]]]), lwd = 2, ylim = c(0, 1), las = 2, main = "", xlab = "")
  title(main = "A. Raw data", xlab = "Log-cpm")
  for (i in seq_along(samples)) {
    den <- density(lcpm[, samples[[i]]])
    lines(den$x, den$y, lwd = 2, col = line_cols[[i]])
  }
  lcpm <- edgeR::cpm(x, log = TRUE)
  plot(density(lcpm[, samples[[1]]]), lwd = 2, ylim = c(0, 1), las = 2, main = "", xlab = "")
  title(main = "B. Filtered data", xlab = "Log-cpm")
  for (i in seq_along(samples)) {
    den <- density(lcpm[, samples[[i]]])
    lines(den$x, den$y, lwd = 2, col = line_cols[[i]])
  }

  # Normalize the data
  x <- edgeR::calcNormFactors(x)
  lcpm <- edgeR::cpm(x, log = TRUE)
  boxplot(lcpm, las = 2)
  title("Normalized data")

  if (!requireNamespace("limma", quietly = TRUE)) {
    stop(
      "To plot MDS plot, 'plotMDS' requires 'limma' package which cannot be ",
      "found. Please install 'limma' using 'BiocManager::install('limma')'",
      call. = FALSE
    )
  }
  limma::plotMDS(lcpm,
    label = x$samples[[group_column]],
    col = as.integer((x$samples[[group_column]])), dim = c(1, 2)
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
#' @details
#' `limma::makeContrasts()` only accepts syntactically valid names, so the group
#' labels are normalised with [make.names()] (the design's coefficient names
#' with them). Labels that are not valid R names, such as the TCGA sample types
#' `"Primary Tumor"` and `"Solid Tissue Normal"`, would otherwise fail: building
#' the contrasts as expressions from those labels either does not parse or
#' silently turns `"Non-tumor" - "Tumor"` into a subtraction of three symbols.
#' The contrasts keep the original labels in their names, e.g.
#' `"Primary TumorvsSolid Tissue Normal"`.
#'
#' @param x A processed `DGEList` object containing normalized count data.
#' @param group_column The name of the column in `x$samples` that contains
#'   the grouping information for the samples.
#'
#' @return An `eBayes` object containing the fitted linear model and
#'   results of the differential expression analysis. The design and contrast
#'   matrices are attached as the `design` and `contrast` attributes.
#' @export
limma_fit <- function(x, group_column) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (!requireNamespace("limma", quietly = TRUE)) {
    stop(
      "To fit linear model, 'limma_fit' requires 'limma' package which cannot ",
      "be found. Please install 'limma' using 'BiocManager::install('limma')'",
      call. = FALSE
    )
  }

  design <- model.matrix(~ 0 + x$samples[[group_column]])
  colnames(design) <- make.names(gsub(".*\\]\\]", "", colnames(design)))

  group_labels <- as.character(unique(x$samples[[group_column]]))
  group_levels <- make.names(group_labels)
  pairs <- combn(seq_along(group_levels), 2)
  contr_strings <- paste0(group_levels[pairs[1, ]], "-", group_levels[pairs[2, ]])

  contr_matrix <- limma::makeContrasts(contrasts = contr_strings, levels = design)
  # makeContrasts() names the columns after the contrast strings, so restore the
  # "AvsB" names built from the original labels
  colnames(contr_matrix) <- paste0(group_labels[pairs[1, ]], "vs", group_labels[pairs[2, ]])

  par(mfrow = c(1, 2))
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
#' @param data_text An optional data frame of genes to label. When given, the
#'   labels are added to the same plot; `x`, `y` and `label` are mapped for it
#'   as well, so the columns must be named like the ones in `data`.
#' @param x variable representing the aesthetic for the x-axis.
#' @param y variable representing the aesthetic for the y-axis.
#' @param color variable representing the column name for the color aesthetic.
#' @param label variable representing the column name for the text label aesthetic.
#'
#' @return A `ggplot` object representing the volcano plot.
#' @export
plot_volcano <- function(data, data_text, x, y, color, label) {
  p <- ggplot(data, aes(x = {{ x }}, y = {{ y }})) +
    geom_point(aes(color = {{ color }})) +
    theme_classic() +
    scale_color_manual(
      values = c("#414788FF", "darkgray", "#22A884FF"),
      guide = guide_legend(override.aes = list(size = 4))
    )

  if (!missing(data_text)) {
    # the labels have to be assigned back, otherwise the plot is built and
    # thrown away, and the caller gets the unlabelled plot
    p <- p + ggrepel::geom_text_repel(
      data = data_text, aes(label = {{ label }}),
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

  p
}

# --- Deprecated aliases ---

#' Construct a DGEList Object (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `dgeList()` was renamed to [dge_list()] to follow the snake_case naming
#' style. It will be removed in version 1.4.0.
#'
#' @param ... Arguments passed on to [dge_list()].
#' @return See [dge_list()].
#' @export
dgeList <- function(...) { # nolint: object_name_linter.
  lifecycle::deprecate_warn("1.3.0", "dgeList()", "dge_list()")
  dge_list(...)
}

#' Filter and Normalize DGEList Data (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `dprocess_dgeList()` was renamed to [process_dge_list()] to follow the
#' snake_case naming style. It will be removed in version 1.4.0.
#'
#' @param ... Arguments passed on to [process_dge_list()].
#' @return See [process_dge_list()].
#' @export
dprocess_dgeList <- function(...) { # nolint: object_name_linter.
  lifecycle::deprecate_warn("1.3.0", "dprocess_dgeList()", "process_dge_list()")
  process_dge_list(...)
}

#' Fit a Linear Model for RNA-seq Data (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `limmaFit()` was renamed to [limma_fit()] to follow the snake_case naming
#' style. It will be removed in version 1.4.0.
#'
#' @param ... Arguments passed on to [limma_fit()].
#' @return See [limma_fit()].
#' @export
limmaFit <- function(...) { # nolint: object_name_linter.
  lifecycle::deprecate_warn("1.3.0", "limmaFit()", "limma_fit()")
  limma_fit(...)
}

#' Plot a Volcano Plot (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `plotVolcano()` was renamed to [plot_volcano()] to follow the snake_case
#' naming style. It will be removed in version 1.4.0.
#'
#' @param ... Arguments passed on to [plot_volcano()].
#' @return See [plot_volcano()].
#' @export
plotVolcano <- function(...) { # nolint: object_name_linter.
  lifecycle::deprecate_warn("1.3.0", "plotVolcano()", "plot_volcano()")
  plot_volcano(...)
}
