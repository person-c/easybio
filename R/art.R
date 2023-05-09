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
art.heatmap <- function(tidy_data, .x, .y, .fill) {
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
art.venn <- function(venn_set, file_name) {
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