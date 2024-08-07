#' Modified plotEnrichment function in R package 'fgsea'
#'
#' @param pathways all pathways.
#' @param pwayname pathway name you want to show.
#' @param stats rank vector.
#' @param gseaParam walk length.
#' @param ticksSize default=0.2.
#' @import data.table
#' @import ggplot2
#' @export
plotEnrichment2 <- function(pathways, pwayname, stats, gseaParam = 1, ticksSize = 0.2) {
  pd <- fgsea::plotEnrichmentData(
    pathway = pathways[[pwayname]], stats = stats,
    gseaParam = gseaParam
  )
  with(pd, ggplot(data = curve) +
    geom_line(aes(x = rank, y = ES),
      color = "blue"
    ) +
    geom_segment(data = ticks, mapping = aes(
      x = rank,
      y = -spreadES / 16, xend = rank, yend = spreadES / 16
    ), linewidth = ticksSize) +
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

#' Visualization of GSEA rank stats.
#'
#' @param stats rank data.
#'
#' @import data.table
#' @import ggplot2
#' @export
plotRank <- function(stats) {
  ranks <- y <- NULL

  rankData <- data.table(ranks = 1:length(stats), y = fsort(stats, TRUE))
  ggplot(data = rankData) +
    scale_x_discrete(expand = expansion(0, 0)) +
    geom_segment(aes(x = ranks, y = 0, xend = ranks, yend = y), color = "gray60") +
    theme_classic(base_size = 6) +
    labs(x = "Rank", y = "Ranked List Metric")
}


#' Visualization of GSEA result from fgsea
#'
#' @param fgseaRes GSEA result from fgsea.
#' @param pathways all pathways.
#' @param pwayname pathway names you want to show.
#' @param stats rank list.
#' @param save Wheather to save plot.
#'
#'
#' @import data.table
#' @import ggplot2
#' @export
plotGSEA <- function(fgseaRes, pathways, pwayname, stats, save = FALSE) {
  . <- NULL
  NES <- padj <- pathway <- NULL

  anno_text <- fgseaRes[.(pwayname), c(NES, padj), on = .(pathway)]
  p1 <- plotEnrichment2(pathways, pwayname, stats)
  p2 <- plotRank(stats)
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

#' Visualization of ORA test.
#'
#' @param data ORA result.
#' @param x x-axis.
#' @param y y-axis.
#' @param size point size.
#' @param fill aesthetic fill. Use constant value for one categroy.
#' @param flip wheather to flip plot
#'
#' @import ggplot2
#' @export
plotORA <- function(data, x, y, size, fill, flip = FALSE) {
  p <- ggplot(
    data,
    aes(
      x = {{ x }},
      y = {{ y }}
    )
  ) +
    geom_col(aes(fill = {{ fill }}), width = 0.5, show.legend = !is.character(fill)) +
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
