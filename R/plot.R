mytheme <-
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 13, color = "black"),
    axis.title = ggplot2::element_text(size = 13),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_line(linewidth = 1),
    axis.ticks.length = ggplot2::unit(3, "mm"),
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(fill = "transparent"),
    legend.direction = "horizontal",
    legend.position = "top",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 10),
    plot.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 10, unit = "pt")
    )


#' plot valcano
#'
#' used to plot valcano
#'
#' @param x differential results.
#' @param fd_name fold change name.
#' @param fd_hold fold change therehold.
#' @param ap_name adjusted p value.
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @export
#' @examples
#' # plot_valcano(data, logFC, 1, adj.P.Val)
#'
plot.limma <- function(x, fd_name, fd_hold, ap_name, ...) {
  data <- x[["diff"]]

  data <- with(data, {
    exp <- ifelse(logFC > 1 & adj.P.Val < 0.05,
      "Up", ifelse(logFC < -1 & adj.P.Val < 0.05, "Down", "Unchanged"))
    data.frame(data, exp)
    }
  )

  ggplot2::ggplot(data, ggplot2::aes(logFC, -log(adj.P.Val, 10))) +
    ggplot2::geom_point(ggplot2::aes(color = exp), size = 1.2) + # nolint
    ggplot2::xlab(expression("log"[2] * "FC")) +
    ggplot2::ylab(expression("-log"[10] * "adj.P.Val")) +
    ggplot2::scale_color_manual(values = c("green", "grey", "red")) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(size = 1.5))
      ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.position = "top",
      legend.text = ggplot2::element_text(size = 10, color = "black"),
      axis.text = ggplot2::element_text(size = 12, color = "black"),
      axis.title = ggplot2::element_text(size = 12, color = "black"),
      axis.ticks.length = ggplot2::unit(3, "mm"),
      axis.line.x.bottom = ggplot2::element_line(),
      axis.line.y.left = ggplot2::element_line(),
    )
}


#' plot survival curve
#'
#' used to plot survival curve
#'
#' @param x surv result
#' @param time month or year.
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @export
#' @examples
#' # plot_surv(fit, data, 'y')
plot.surv <- function(x, time = 'y', ...) {

  fit <- x$fit
  data <- x$data

  attr <- strsplit(x = attributes(fit$strata)$names, split = "=")
  legend_labs <- c(attr[[1]][[2]], attr[[2]][[2]])

  raw_curve <- survminer::ggsurvplot(
    fit,
    data,
    size = 1,
    palette = c("#440154", "#0d0887"),
    conf.int = TRUE,
    pval = TRUE,
    xlab = paste0('Time(', time, ')'),
    risk.table = TRUE,
    risk.table.col = "strata",
    legend.labs = legend_labs,
    risk.table.height = 0.25,

    ggtheme = ggplot2::theme_bw() +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
      ),

    tables.theme = ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect()
    )
  )

  raw_curve
}



#' plot gsea plot
#'
#' used to gsea
#'
#' @param x gsea results
#' @param name pathway name
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @export
#' @examples
#' # plot_gsea()
plot.gsea <- function(x, name, ...) {
  fgsea::plotEnrichment(x$pathways[[name]], x$ranks) +
  ggplot2::labs(title = name)
}

#' plot go_rich plot
#'
#' used to plot GO rich colplot
#'
#' @param x go rich result
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @export
#' @examples
#' # no example
plot.go <- function(x, ...) {

  x <- x@result
  x <- by(x, INDICES = x$ONTOLOGY, FUN = function(x) head(x, n = 8))
  x <- rlang::exec('rbind', !!!x)
  x <- transform(x, Description = factor(Description, levels = Description))

  ggplot2::ggplot(x,
    ggplot2::aes(
      x = Description,
      y = -log10(p.adjust),
      fill = ONTOLOGY)
  ) +
  ggplot2::facet_wrap(~ ONTOLOGY, scales = "free_x", nrow = 1) +
  ggplot2::geom_col(width = 1, ggplot2::aes(colour = "black")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(0)) +
  ggplot2::scale_color_manual(values = "black") +
  ggplot2::labs(y = "-log10(p.adjust)") +
  mytheme +
  ggplot2::theme(
    aspect.ratio = 15 / 10,
    axis.title.y = ggplot2::element_text(size = 8),
    axis.text.x = ggplot2::element_text(
      size = 8, angle = 90, hjust = 1, vjust = .5
      ),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.length = ggplot2::unit(2, "mm"),
    axis.ticks.x = ggplot2::element_line(linewidth = .5),
    legend.position = "none"
    )
}



#' plot kegg rich plot
#'
#' used to plot kegg rich colplot
#'
#' @param x kegg result
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @export
#' @examples
#' # plot_gorich()
plot.kegg <- function(x, ...) {

  x <- x@result
  x <- head(subset(x, p.adjust < 0.05), n = 8)
  x <- transform(x, Description = factor(Description, levels = Description))

  ggplot2::ggplot(x,
    ggplot2::aes(
      x = GeneRatio,
      y = Description,
      color = p.adjust)
    ) +
  ggplot2::geom_point(ggplot2::aes(size = Count)) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      size = 8, angle = 90, hjust = 1, vjust = .5)
    )
}
