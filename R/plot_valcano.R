#' plot valcano
#'
#' used to plot valcano 
#'
#' @param data differential results.
#' @param fd_name fold change name.
#' @param fd_hold fold change therehold.
#' @param ap_name adjusted p value.
#'
#' @return ggplot2 object
#' @examples
#' # plot_valcano(data, logFC, 1, adj.P.Val)
#'
plot_valcano <- function(data, fd_name, fd_hold, ap_name) {
  data |>
  dplyr::mutate(
    Expression = dplyr::case_when(
      {{fd_name}} >= fd_hold & {{ap_name}} <= 0.05 ~ "Up-regulated",
      {{fd_name}} <= - fd_hold & {{ap_name}} <= 0.05 ~ "Down-regulated",
      TRUE ~ "Unchanged"),
    Significance = dplyr::case_when(
      abs({{fd_name}}) >= fd_hold & {{ap_name}} <= 0.05 &
        {{ap_name}} > 0.01 ~ "P.val 0.05",
      abs({{fd_name}}) >= fd_hold & {{ap_name}} <= 0.01 &
        {{ap_name}} > 0.001 ~ "P.val 0.01",
      abs({{fd_name}}) >= fd_hold & {{ap_name}} <= 0.001 ~ "P.val 0.001",
      TRUE ~ "Unchanged")
  ) |>

  ggplot2::ggplot(ggplot2::aes({{fd_name}}, -log({{ap_name}}, 10))) +
    ggplot2::geom_point(ggplot2::aes(color = Significance), size = 1.2) +
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
