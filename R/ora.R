#' Plot Enrichment for a Specific Pathway in fgsea
#'
#' This function creates a plot of enrichment scores for a specified pathway.
#' It provides a visual representation of the enrichment
#' score (ES) along with the ranks and ticks indicating the GSEA walk length.
#'
#' @param pathways A list of pathways.
#' @param pwayname The name of the pathway for which to plot enrichment.
#' @param stats A rank vector obtained from the 'fgsea' package.
#' @param gsea_param The GSEA walk length parameter. Default is 1.
#' @param ticks_size The size of the tick marks. Default is 0.2.
#'
#' @return A ggplot object representing the enrichment plot.
#' @import ggplot2
#' @export
plot_enrichment <- function(pathways, pwayname, stats, gsea_param = 1, ticks_size = 0.2) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("To get plot data, plot_enrichment() requires 'fgsea' package which cannot be found. Please install 'fgsea' using 'BiocManager::install('fgsea')'.")
  }
  pd <- fgsea::plotEnrichmentData(
    pathway = pathways[[pwayname]], stats = stats,
    gsea_param = gsea_param
  )
  with(pd, ggplot(data = curve) +
    geom_line(aes(x = rank, y = ES),
      color = "blue"
    ) +
    geom_segment(data = ticks, mapping = aes(
      x = rank,
      y = -spreadES / 16, xend = rank, yend = spreadES / 16
    ), linewidth = ticks_size) +
    geom_hline(
      yintercept = posES, colour = "gray90",
      linetype = "dashed", linewidth = 0.2
    ) +
    geom_hline(
      yintercept = negES, colour = "gray90",
      linetype = "dashed", linewidth = 0.2
    ) +
    geom_hline(yintercept = 0, colour = "black") +
    scale_x_discrete(expand = expansion(0, 0)) +
    labs(y = "Enrichment Score")) +
    theme(
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      panel.grid.major = element_line(color = "grey92"),
      text = element_text(size = 6)
    )
}

#' @title Visualization of GSEA Rank Statistics
#'
#' @description
#' The `plot_rank` function visualizes the ranked statistics of a GSEA (Gene Set Enrichment Analysis) analysis.
#' The function creates a plot where the x-axis represents the rank of each gene, and the y-axis shows
#' the corresponding ranked list metric.
#'
#' @param stats A numeric vector containing the ranked statistics from a GSEA analysis.
#'
#' @import ggplot2
#' @return ggplot2 object
#' @export
plot_rank <- function(stats) {
  ranks <- y <- NULL

  rank_data <- data.table(ranks = 1:length(stats), y = fsort(stats, TRUE))
  ggplot(data = rank_data) +
    scale_x_discrete(expand = expansion(0, 0)) +
    geom_segment(aes(x = ranks, y = 0, xend = ranks, yend = y), color = "gray60") +
    theme_classic(base_size = 6) +
    labs(x = "Rank", y = "Ranked List Metric")
}


#' @title Visualization of GSEA Result from [fgsea::fgsea()]
#'
#' @description
#' The `plot_gsea` function visualizes the results of a GSEA (Gene Set Enrichment Analysis) using data from
#' the `fgsea` package. It generates a composite plot that includes an enrichment plot and a ranked metric plot.
#'
#' @param fgsea_res A data table containing the GSEA results from the `fgsea` package.
#' @param pathways A list of all pathways used in the GSEA analysis.
#' @param pwayname The name of the pathway to visualize.
#' @param stats A numeric vector representing the ranked statistics.
#' @param save A logical value indicating whether to save the plot as a PDF file. Default is `FALSE`.
#'
#' @import ggplot2
#' @return ggplot2 object.
#' @export
plot_gsea <- function(fgsea_res, pathways, pwayname, stats, save = FALSE) {
  . <- NULL
  NES <- padj <- pathway <- NULL

  anno_text <- fgsea_res[.(pwayname), c(NES, padj), on = .(pathway)]
  p1 <- plot_enrichment(pathways, pwayname, stats)
  p2 <- plot_rank(stats)
  p3 <- patchwork::wrap_plots(p1, p2, heights = c(0.8, 0.3)) +
    patchwork::plot_annotation(
      title = pwayname,
      subtitle = sprintf("NES: %f, padj: %f", anno_text[[1]], anno_text[[2]]),
      theme = theme(
        plot.title = element_text(size = 8),
        plot.subtitle = element_text(size = 5, face = "italic")
      )
    )

  if (save) ggsave(paste0(pwayname, ".pdf"), width = 4, height = 2.7)

  p3
}

#' @title Visualization of ORA Test Results
#'
#' @description
#' The `plot_ora` function visualizes the results of an ORA (Over-Representation Analysis) test.
#' It generates a plot with customizable aesthetics for x, y, point size, and fill, with an option to flip the axes.
#'
#' @param data A data frame containing the ORA results to be visualized.
#' @param x The column in `data` to map to the x-axis.
#' @param y The column in `data` to map to the y-axis.
#' @param size The column in `data` to map to the size of the points.
#' @param fill The column in `data` to map to the fill color of the bars or points.
#'        Use a constant value for a single category.
#' @param flip A logical value indicating whether to flip the axes of the plot. Default is `FALSE`.
#'
#' @import ggplot2
#' @return ggplot2 object.
#' @export
plot_ora <- function(data, x, y, size, fill, flip = FALSE) {
  is_character <- try(is.character(fill), silent = TRUE)
  if (inherits(is_character, "try-error")) is_character <- FALSE

  p <- ggplot(
    data,
    aes(
      x = {{ x }},
      y = {{ y }}
    )
  ) +
    geom_col(aes(fill = {{ fill }}), width = 0.5, show.legend = !is_character) +
    geom_point(aes(size = {{ size }})) +
    scale_x_continuous(expand = expansion(add = c(0, 1))) +
    scale_fill_brewer(palette = "Paired") +
    theme_classic()

  if (flip) {
    p <- p + scale_y_discrete(
      expand = expansion(0, 0),
      guide = guide_axis(angle = 80)
    ) + coord_flip()
  }

  p
}

# Deprecated aliases ----------------------------------------------------------

#' Plot Enrichment for a Pathway (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `plotEnrichment2()` was renamed to [plot_enrichment()] to follow the
#' snake_case naming style. It will be removed in the next version.
#'
#' @param ... Arguments passed on to [plot_enrichment()].
#' @return See [plot_enrichment()].
#' @export
plotEnrichment2 <- function(...) {
  lifecycle::deprecate_warn("1.2.4", "plotEnrichment2()", "plot_enrichment()")
  plot_enrichment(...)
}

#' Visualize GSEA Rank Statistics (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `plotRank()` was renamed to [plot_rank()] to follow the snake_case naming
#' style. It will be removed in the next version.
#'
#' @param ... Arguments passed on to [plot_rank()].
#' @return See [plot_rank()].
#' @export
plotRank <- function(...) {
  lifecycle::deprecate_warn("1.2.4", "plotRank()", "plot_rank()")
  plot_rank(...)
}

#' Visualize GSEA Results (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `plotGSEA()` was renamed to [plot_gsea()] to follow the snake_case naming
#' style. It will be removed in the next version.
#'
#' @param ... Arguments passed on to [plot_gsea()].
#' @return See [plot_gsea()].
#' @export
plotGSEA <- function(...) {
  lifecycle::deprecate_warn("1.2.4", "plotGSEA()", "plot_gsea()")
  plot_gsea(...)
}

#' Visualize ORA Test Results (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `plotORA()` was renamed to [plot_ora()] to follow the snake_case naming
#' style. It will be removed in the next version.
#'
#' @param ... Arguments passed on to [plot_ora()].
#' @return See [plot_ora()].
#' @export
plotORA <- function(...) {
  lifecycle::deprecate_warn("1.2.4", "plotORA()", "plot_ora()")
  plot_ora(...)
}
