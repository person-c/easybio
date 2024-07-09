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
#' @param .fill the filled value for ggplot2.
#' @import data.table
#' @import ggplot2
#' @export
view.rich <- function(data, y, n = 8, .fill) {
  x <- data.table::setDT(copy(data))
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
      ggplot2::geom_col(fill = .fill, width = 0.5) +
      coord_flip()
  }
  p <- p +
    ggplot2::geom_point(aes(size = log10(Count))) +
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
