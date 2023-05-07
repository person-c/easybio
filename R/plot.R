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
#' @param data differential results.
#' @param fd_name fold change name.
#' @param fd_hold fold change therehold.
#' @param ap_name adjusted p value.
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @examples
#' # plot_valcano(data, logFC, 1, adj.P.Val)
#'
plot.limma <- function(data, fd_name, fd_hold, ap_name, ...) {

  data <- data[["diff"]]

  data <- with(data, {
    exp <- ifelse(logFC > 1 & adj.P.Val < 0.05,
      "Up", ifelse(logFC < -1 & adj.P.Val < 0.05, "Down", "Unchanged"))
    data.frame(data, exp)
    }
  )

  ggplot2::ggplot(data, ggplot2::aes(logFC, -log(adj.P.Val, 10))) +
    ggplot2::geom_point(ggplot2::aes(color = exp), size = 1.2) + # nolint
    ggplot2::xlab(expression("log"[2] * "FC")) +
    ggplot2::ylab(expression("-log"[10] * "P.Val")) +
    ggplot2::scale_color_manual(values = c("red", "blue", "green", "grey")) +
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
#' @param fit surv fit result.
#' @param data fit data.
#' @param time month or year.
#' @param ... additional arguments
#'
#' @return ggplot2 object
#' @examples
#' # plot_surv(fit, data, 'y')
plot.surv <- function(fit, data, time, ...) {

  force(fit)
  force(data)
  force(time)

  attr <- strsplit(x = attributes(fit$strata)$names, split = "=")
  legend_labs <- c(attr[[1]][[2]], attr[[2]][[2]])

  raw_curve <- survminer::ggsurvplot(
    fit,
    data,
    size = 1,
    palette = c("#440154", "#0d0887"),
    conf.int = TRUE,
    pval = TRUE,
    xlab = glue::glue("Time({time})"),
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
#' @param pathways gsea results
#' @param rank gene vector
#' @param ... additional arguments
#'
#' @return figure
#' @examples
#' # plot_gsea()
plot.gsea <- function(pathways, rank, ...) {
  fgsea::plotEnrichment(pathways, rank) +
  ggplot2::labs(title = "Programmed Cell Death")
}

#' plot go_rich plot
#'
#' used to plot GO rich colplot
#'
#' @param go_result go rich result
#' @param ... additional arguments
#' @importFrom magrittr `%>%`
#'
#' @return figure
#' @export
#' @examples
#' # plot_gorich()
plot.go <- function(go_result, ...) {

  go_result <- data.table::setDT(go_result@result)
  go_result[order(p.adjust), head(.SD, 8), by = ONTOLOGY] %>%  #nolint

  dplyr::group_by(ONTOLOGY) %>%
  dplyr::mutate(Description = forcats::fct_relevel(Description, Description)) %>% # nolint
  dplyr::ungroup() %>%
    ggplot2::ggplot(
      ggplot2::aes(x = Description, y = -log10(p.adjust), fill = ONTOLOGY)) + # nolint
    ggplot2::facet_wrap(~ ONTOLOGY, scales = "free_x", nrow = 1) +
    ggplot2::geom_col(width = 1, ggplot2::aes(colour = "black")) +
    ggplot2::labs(y = "-log10(p.adjust)") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(0)) +
    viridis::scale_fill_viridis(discrete = TRUE) +
    ggplot2::scale_color_manual(values = "black") +
    mytheme + # nolint
    ggplot2::theme(
      aspect.ratio = 15 / 10,
      axis.title.y = ggplot2::element_text(size = 8),
      axis.text.x = ggplot2::element_text(
        size = 8, angle = 90, hjust = 1, vjust = .5),
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
#' @param kegg_result kegg result
#' @param ... additional arguments
#' @importFrom magrittr `%>%`
#'
#' @return figure
#' @examples
#' # plot_gorich()
plot.kegg <- function(kegg_result, ...) {



  kegg_result <- kegg_result@result
  kegg_result <- data.table::setDT(kegg_result)

  kegg_result[p.adjust < 0.05][order(GeneRatio), head(.SD, 8)] |>
  tibble::as_tibble() %>%
  dplyr::mutate(Description = forcats::fct_relevel(Description, Description)) %>%

    ggplot2::ggplot(ggplot2::aes(x = GeneRatio, y = Description,
      color = p.adjust)) +
    ggplot2::geom_point(ggplot2::aes(size = Count)) + # nolint
    viridis::scale_color_viridis(discrete = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        size = 8, angle = 90, hjust = 1, vjust = .5)
      )
}


#' plot heatmap
#'
#' used to plot heatmap
#'
#' @param tidy_data long data.
#' @param .x x axis.
#' @param .y  y axis.
#' @param .fill fill column
#'
#' @return ggplot2 object
#' @export
#' @examples
#' # plot_heatmap(data, x, y, group)
plot_heatmap <- function(tidy_data, .x, .y, .fill) {
  ggplot2::ggplot(tidy_data,
    ggplot2::aes(x = {{.x}}, y = {{.y}}, fill = {{.fill}})) +
  ggplot2::geom_raster() +
  ggplot2::scale_x_discrete(expand = ggplot2::expansion(0)) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(0)) +
  ggplot2::scale_fill_gradient2(
    midpoint = 0,
    limit = c(-1, 1)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 90)
  )
}

#' plot venn
#'
#' used to plot venn plot
#'
#' @param venn_set a list containing set
#' @param file_name file saved name
#'
#' @return figure
#' @export
#' @examples
#' # plot_venn(list, 'venn')
plot_venn <- function(venn_set, file_name) {
  VennDiagram::venn.diagram( # nolint
    venn_set,
    filename = file_name,
    disable.logging = TRUE,
    height = 800,
    width = 1000,
    resolution = 500,
    imagetype = "png",
    output = TRUE,
    compression = "lzw",
    lwd = .2,
    col = c("black", "black"),
    fill = c("#440154ff", "#006600"),
    cex = 0.3,
    fontfamily = "sans",
    cat.cex = 0.5,
    cat.default.pos = "text",
    cat.fontfamily = "sans",
    cat.pos = c(-35, 27),
    cat.dist = c(0.05, 0.05)
  )
}