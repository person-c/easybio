# #' GO ENRICHMENT
# #'
# #' @param degs processed DGEList.
# #' @param .subset subset expression.
# #' @param .column gene name column.
# #' @param org species.
# #' @param pathway defined gene sets.
# #' @import data.table
# #' @importFrom clusterProfiler enricher
# #' @importFrom clusterProfiler enrichGO
# richGO <- function(degs, .subset, .column, org, pathway = NULL) {
#   sig <- eval(substitute(subset(degs, subset = .subset, .column, drop = TRUE)))
#   message(sprintf("%d genes used to do GO enrichment", length(sig)))
#   if (!is.null(pathway)) {
#     return(clusterProfiler::enricher(sig, TERM2GENE = pathway, pvalueCutoff = 1, qvalueCutoff = 1)@result)
#   }
#   res <- clusterProfiler::enrichGO(sig, org, keyType = "SYMBOL", ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1)
#   res@result
# }

#' GO ENRICHMENT VISIULIZATION
#'
#' @param data GO result data.frame.
#' @param y significant symbol.
#' @param n number of terms to show.
#' @import data.table
#' @export
view.go <- function(data, y, n = 8) {
  # x <- x@result
  x <- data.table::setDT(copy(data))

  if ("ONTOLOGY" %in% colnames(x)) {
    x <- eval(substitute(x[, head(.SD[order(y)], n), keyby = ONTOLOGY]))
    x[, Description := factor(Description, levels = Description)]
    p <- ggplot2::ggplot(
      x,
      ggplot2::aes(
        y = -log10({{ y }}),
        x = Description,
        fill = ONTOLOGY
      )
    ) +
      ggplot2::facet_grid(. ~ ONTOLOGY, scales = "free_x")
  } else {
    x <- eval(substitute(x[, head(.SD[order(y)], n)]))
    x[, Description := factor(Description, levels = Description)]
    p <- ggplot2::ggplot(
      x,
      ggplot2::aes(
        y = -log10({{ y }}),
        x = Description
      )
    )
  }

  p <- p +
    ggplot2::geom_col(width = 1) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(0.01))
  p + scale_x_discrete(guide = guide_axis(angle = 75))
}


# #' KEGG ENRICHMENT
# #'
# #' @param degs degs.
# #' @param .subset subset expression.
# #' @param .column gene name column.
# #' @param  org species
# #' @param pathway gene sets.
# #' @import data.table
# richKEGG <- function(degs, .subset, .column, org, pathway = NULL) {
#   sig <- eval(substitute(subset(degs, subset = .subset, .column, drop = TRUE)))
#   message(sprintf("%d genes used to do KEGG enrichment", length(sig)))
#   if (!is.null(pathway)) {
#     return(clusterProfiler::enricher(sig, TERM2GENE = pathway, pvalueCutoff = 1, qvalueCutoff = 1)@result)
#   }
#   sig <- clusterProfiler::bitr(sig, "SYMBOL", "ENTREZID", org)
#   sig <- sig[[2]]
#   res <- clusterProfiler::enrichKEGG(sig, pvalueCutoff = 1, qvalueCutoff = 1)

#   res@result
# }

#' KEGG ENRICHMENT VISUALIZATION
#'
#' @param data KEGG result.
#' @param color color mapping.
#' @param n number of terms to show.
#' @import data.table
#' @import ggplot2
#' @export
view.kegg <- function(data, color, n = 8) {
  dt <- data.table::setDT(copy(data))
  if (nrow(dt) < n) n <- nrow(dt)
  dt <- eval(substitute(dt[, head(.SD, n)][order(color)]))

  ggplot(dt, aes(x = GeneRatio, y = Description, color = {{ color }})) +
    geom_point(aes(size = Count)) +
    scale_color_distiller(palette = "YlOrRd") +
    theme_bw()
}

#' Visualization of enrichment result
#'
#' @param data enrichment dataframe.
#' @param y the y value's aesthetic.
#' @param n number of terms to show.
#' @param fill the filled value for ggplot2.
#' @param size size aes.
#' @import data.table
#' @import ggplot2
#' @export
view.rich <- function(data, y, n = 8, fill, size) {
  x <- data.table::setDT(data)
  if ("ONTOLOGY" %in% colnames(x)) {
    x <- eval(substitute(x[, head(.SD[order(y)], n), keyby = ONTOLOGY]))
    x[, `:=`(Description, factor(Description, levels = rev(Description)))]
    p <- ggplot2::ggplot(x, ggplot2::aes(y = -log10({{ y }}), x = Description, fill = ONTOLOGY)) +
      coord_flip() +
      ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free_y") +
      geom_col(aes(fill = ONTOLOGY), width = 0.5) +
      scale_fill_brewer(palette = "Paired")
  } else {
    x <- eval(substitute(x[, head(.SD[order(y)], n)]))
    x[, `:=`(Description, factor(Description, levels = rev(Description)))]
    p <- ggplot2::ggplot(x, ggplot2::aes(y = -log10({{ y }}), x = Description)) +
      ggplot2::geom_col(fill = fill, width = 0.5) +
      coord_flip()
  }
  p <- p +
    ggplot2::geom_point(aes(size = log10({{ size }}))) +
    ggplot2::scale_size_continuous(range = c(0.5, 3)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(0, 0)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(0, 1.01)))
  p <- p + scale_x_discrete() +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      text = element_text(family = "sans"),
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.spacing = unit(0, "cm")
    )

  p
}


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
  pd <- plotEnrichmentData(
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
