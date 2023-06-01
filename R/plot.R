#' Plot valcano
#'
#' Used to plot valcano for limma analysis.
#'
#' @param x differential results.
#' @param fd_name fold change name.
#' @param fd_hold fold change therehold.
#' @param ap_name adjusted p value.
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @importFrom data.table `:=`
#' @export
#' @examples
#' data(expr)
#' y <- analyze(expr, "limma", "cc", "array")
#' plot(y)
plot.limma <- function(x, fd_name, fd_hold, ap_name, ...) {
  data <- x$diff

  lfun <- function(x, y) {
    ifelse(x > 1 & y < 0.05, "Up", ifelse(x < -1 & y < 0.05, "Down", "Unchanged"))
  }

  data <- data.table::setDT(data)
  data[, exp := list(val = lfun(logFC, adj.P.Val))]

  ggplot2::ggplot(data, ggplot2::aes(logFC, -log(adj.P.Val, 10))) +
    ggplot2::geom_point(ggplot2::aes(color = exp), size = 1.2) +
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


#' Plot survival curve
#'
#' Used to plot survival curve for surv task.
#'
#' @param x surv result
#' @param time month or year.
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @export
#' @examples
#' library(survival)
#' y <- analyze(lung, "surv", Surv(time, status) ~ sex)
#' plot(y, time = "y")
plot.surv <- function(x, time = "y", ...) {
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
    xlab = paste0("Time(", time, ")"),
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
#' library(fgsea)
#' data(examplePathways)
#' data(exampleRanks)
#' set.seed(42)
#' y <- analyze(exampleRanks, "gsea", examplePathways)
#' plot(y, name = "5991130_Programmed_Cell_Death")
plot.gsea <- function(x, name, ...) {
  fgsea::plotEnrichment(x$pathways[[name]], x$ranks) +
    ggplot2::labs(title = name)
}

#' plot go_rich plot
#'
#' used to plot GO rich colplot
#'
#' @param x go rich result
#' @param n number of pathways you want to show
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @export
#' @examples
#' data(gene_vector)
#' library(org.Hs.eg.db)
#' y <- analyze(gene_vector, "go",
#'   db = org.Hs.eg.db,
#'   from = "SYMBOL",
#'   to = "ENSEMBL",
#'   ont = "ALL"
#' )
#' plot(y)
plot.go <- function(x, n = 8, ...) {
  x <- x@result
  x <- data.table::setDT(x)

  if ("ONTOLOGY" %in% colnames(x)) {
    x <- x[, head(.SD, n), keyby = ONTOLOGY]
    p <- ggplot2::ggplot(
      x,
      ggplot2::aes(
        x = -log10(p.adjust),
        y = Description,
        fill = ONTOLOGY
      )
    ) +
      ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free_y")
  } else {
    x <- x[, head(.SD, n)]
    p <- ggplot2::ggplot(
      x,
      ggplot2::aes(
        x = -log10(p.adjust),
        y = Description,
        fill = "constant"
      )
    )
  }

  x[, Description := list(val = factor(Description, levels = rev(Description)))]

  p <- p +
    ggplot2::geom_col(width = 1, ggplot2::aes(colour = "black")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(0)) +
    ggplot2::scale_color_manual(values = "black") +
    ggplot2::labs(x = "-log10(p.adjust)") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.ticks.x = ggplot2::element_line(linewidth = .5),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks.length = ggplot2::unit(2, "mm"),
      panel.background = ggplot2::element_blank(),
      legend.position = "none"
    )

  return(p)
}



#' Plot kegg rich plot
#'
#' used to plot kegg rich colplot for kegg task.
#'
#' @param x kegg result
#' @param n number of pathways you want to show
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @export
#' @examples
#' library(org.Hs.eg.db)
#' data(kegg)
#' y <- analyze(
#'   object = kegg, task = "kegg", db = org.Hs.eg.db,
#'   from = "SYMBOL",
#'   to = "ENTREZID"
#' )
plot.kegg <- function(x, n = 8, ...) {
  x <- x@result
  x <- data.table::setDT(x)
  x <- x[, head(.SD, n)][order(GeneRatio)]
  x[, Description := list(val = factor(Description, levels = Description))]

  ggplot2::ggplot(
    x,
    ggplot2::aes(
      x = GeneRatio,
      y = Description,
      color = p.adjust
    )
  ) +
    ggplot2::geom_point(ggplot2::aes(size = Count)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        size = 8, angle = 90, hjust = 1, vjust = .5
      )
    )
}


#' plot forest plot
#'
#' used to plot forest plot
#'
#' @param x coxph result
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @export
#' @examples
#' library(survival)
#' y <- analyze(lung, "cox", Surv(time, status) ~ sex)
#' plot(y, time = "y")
plot.cox <- function(x, ...) {
  survminer::ggforest(model = x$fit, data = x$data)
}
