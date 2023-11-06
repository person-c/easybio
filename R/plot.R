#' Plot valcano
#'
#' Used to plot valcano for limma analysis.
#'
#' @param x differential results.
#' @param fd_name fold change name.
#' @param fd_hold fold change therehold.
#' @param ap_name adjusted p value.
#' @param ap_hold adjusted p value therehold.
#' @param top top number genes to show lables.
#' @param ... additional arguments
#' @return ggplot2 object
#' @importFrom data.table `:=`
#' @export
#' @examples
#' data(expr)
#' y <- analyze(expr, "limma", "cc", "array")
#' plot(y, logFC, 1, adj.P.Val, 0.05, 10)
plot.limma <- function(x, fd_name, fd_hold, ap_name, ap_hold, top = FALSE, ...) {
  lfun <- function(x, y) {
    ifelse(x > fd_hold & y < ap_hold, "Up",
      ifelse(x < -fd_hold & y < ap_hold, "Down", "Unchanged")
    )
  }

  dt <- data.table::setDT(x$diff, keep.rownames = TRUE)
  eval(substitute(dt[, exp := list(val = lfun(fd_name, ap_name))]))


  p <- ggplot2::ggplot(
    data = dt,
    ggplot2::aes(x = {{ fd_name }}, y = -log10({{ ap_name }}), color = exp)
  ) +
    ggplot2::geom_point(size = 1.2, alpha = 0.4, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c(Down = "seagreen", Unchanged = "darkgray", Up = "firebrick2"),
      guide = guide_legend(override.aes = list(size = 4))
    ) +
    ggplot2::scale_x_continuous(
      name = expression("log"[2] * "FC")
    ) +
    ggplot2::scale_y_continuous(name = expression("-log"[10] * "adj.P.Val")) +
    ggplot2::geom_vline(
      xintercept = c(-fd_hold, fd_hold), lty = 4, col = "darkgray", lwd = 0.6
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(0.05), lty = 4, col = "darkgray", lwd = 0.6
    ) +
    ggplot2::theme_bw(base_size = 12, base_family = "Times") +
    ggplot2::theme(
      legend.position = "top",
      panel.grid = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      text = element_text(
        face = "bold", color = "black", size = 18
      )
    )

  if (top) {
    dttp <- eval(substitute(dt[c(
      order(logFC)[1:top],
      order(-logFC)[1:top]
    )]))
    p + ggrepel::geom_label_repel(
      data = dttp, aes(label = rn),
      max.overlaps = 20, size = 4,
      box.padding = unit(0.5, "lines"), min.segment.length = 0,
      point.padding = unit(0.8, "lines"),
      segment.color = "black", show.legend = FALSE
    )
  }
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
#' @param name pathway you want to plot
#' @param top top p.adjust number to show
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @export
plot.gsea <- function(x, name = FALSE, top = 5, ...) {
  if (name) {
    fgsea::plotEnrichment(x$pathways[[name]], x$ranks) +
      ggplot2::labs(title = name)
  } else {
    sortedgsea <- data.table::as.data.table(x[[1]])[order(p.adjust)]
    enrichplot::gseaplot2(x[[1]], sortedgsea[["ID"]][1:top],
      base_size = 10,
      color = c(
        "#7B68EE", "#CD3333", "#20B2AA",
        "#FF8C00", "#FF6666", "#8CC785"
      ),
      rel_heights = c(1.5, 0.3, 0.5)
    )
  }
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


#' Visualize PPI analysis
#'
#' used to plot PPI network
#'
#' @param x PPI analyzze result
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @export
plot.ppi <- function(x, ...) {
  x <- tidygraph::as_tbl_graph(x)
  igraph::V(x)$deg <- igraph::degree(x)
  igraph::V(x)$size <- igraph::degree(x) / 5
  igraph::E(x)$width <- igraph::E(x)$weight / 10
  ggraph::ggraph(x, layout = "linear", circular = TRUE) +
    ggraph::geom_edge_fan(ggplot2::aes(edge_width = width),
      color = "lightblue", show.legend = FALSE
    ) +
    ggraph::geom_node_point(aes(size = size), color = "orange", alpha = 0.7) +
    ggraph::geom_node_text(aes(filter = deg > 0, label = name),
      size = 5, repel = TRUE
    ) +
    ggraph::scale_edge_width(range = c(0.2, 1)) +
    ggplot2::scale_size_continuous(range = c(1, 10)) +
    ggplot2::guides(size = FALSE) +
    ggraph::theme_graph()
}
