#' Custom ggplot2 Theme for Academic Publications
#'
#' `theme_publication` creates a custom ggplot2 theme designed for academic publications, ensuring clarity, readability, and a professional appearance.
#' It is based on `theme_classic()` and includes additional refinements to axis lines, text, and other plot elements to meet the standards of high-quality academic figures.
#'
#' @param base_size numeric, the base font size. Default is 12.
#' @param base_family character, the base font family. Default is "sans".
#'
#' @return A ggplot2 theme object that can be applied to ggplot2 plots.
#' @import ggplot2
#' @export
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   theme_publication()
#' print(p)
theme_publication <- function(base_size = 12, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text = element_text(color = "black", size = base_size * 0.8),
      axis.title = element_text(color = "black", size = base_size * 0.9, face = "bold"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = base_size * 0.8),
      plot.title = element_text(size = base_size * 1.2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size * 0.9, hjust = 0.5),
      plot.caption = element_text(size = base_size * 0.8, hjust = 1),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      complete = TRUE
    )
}
