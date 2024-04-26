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
