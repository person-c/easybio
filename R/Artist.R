#' @title Reductor Class
#'
#' @description
#' The base class to plot.
#' @import ggplot2
#' @export
Artist <- R6::R6Class("Artist",
  public = list(
    #' @description
    #' dumbbbell plot
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param col map to col
    #' @param ... additional aesthetics properties mapping
    dumbbbell = function(data, x, y, col, ...) {
      ggplot(data, aes(x = {{ x }}, y = {{ y }}), ...) +
        geom_line() +
        geom_point(aes(col = {{ col }}), size = 3)
    },
    #' @description
    #' bubble plot
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param size map to size
    #' @param col map to col
    #' @param ... additional aesthetics properties mapping
    bubble = function(data, x, y, size, col, ...) {
      ggplot(
        data,
        aes(
          x = {{ x }}, y = {{ y }},
          size = {{ size }}, col = {{ col }},
          ...
        )
      ) +
        geom_point() +
        scale_size(name = "Size", range = c(1, 10))
    },
    #' @description
    #' divergence barchart
    #' @param data data
    #' @param group map to group
    #' @param y map to y-axis
    #' @param fill map to fill
    #' @param ... additional aesthetics properties mapping
    barchart_divergence = function(data, group, y, fill, ...) {
      ggplot(
        df,
        aes(
          x = reorder({{ group }}, {{ y }}),
          y = {{ y }}, ...
        )
      ) +
        geom_bar(
          stat = "identity",
          show.legend = FALSE,
          fill = ifelse(y >= 0, "lightblue", "lightpink"),
          col = "white"
        ) +
        geom_hline(yintercept = 0, col = 1, lwd = 0.2) +
        geom_text(aes(
          label = group, # Text with groups
          hjust = ifelse(value < 0, 1.5, -1),
          vjust = 0.5
        ), size = 2.5) +
        xlab("Group") +
        ylab("Value") +
        scale_y_continuous(
          breaks = seq(-2, 2, by = 1),
          limits = c(-2.5, 2.5)
        ) +
        coord_flip() +
        theme_minimal() +
        theme(
          axis.texty = element_blank(), # Remove Y-axis texts
          axis.ticksy = element_blank(), # Remove Y-axis ticks
          panel.grid.majory = element_blank()
        ) # Remove horizontal grid
    },
    #' @description
    #' lollipot plot
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param ... additional aesthetics properties mapping
    lollipop = function(data, x, y, ...) {
      ggplot(data, aes(x = {{ x }}, y = {{ y }}, ...)) +
        geom_segment(aes(x = {{ x }}, xend = {{ x }}, y = 0, yend = {{ y }}),
          col = "gray", lwd = 1
        ) +
        geom_point(size = 7.5, pch = 21, bg = 4, col = 1) +
        geom_text(aes(label = y), col = "white", size = 3) +
        scale_x_discrete(labels = paste0("G_", 1:10)) +
        coord_flip() +
        theme_minimal()
    },
    #' @description
    #' contour
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param ... additional aesthetics properties mapping
    contour = function(data, x, y, ...) {
      ggplot(data, aes(x = {{ x }}, y = {{ y }}, ...)) +
        geom_point() +
        geom_density_2d_filled(alpha = 0.4) +
        geom_density_2d(colour = "black")
    },
    #' @description
    #' scatter plot with ellipses
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param col map to col
    #' @param ... additional aesthetics properties mapping
    scatter_ellipses = function(data, x, y, col, ...) {
      ggplot(data, aes(
        x = {{ x }},
        y = {{ y }}, col = {{ col }}, ...
      )) +
        geom_point() +
        stat_ellipse(
          geom = "polygon",
          aes(fill = group),
          alpha = 0.25
        )
    },
    #' @description
    #' donut plot
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param fill map to fill
    #' @param ... additional aesthetics properties mapping
    donut = function(data, x, y, fill, ...) {
      hsize <- 3
      ggplot(data, aes(
        x = {{ x }}, y = {{ y }},
        fill = {{ fill }}, ...
      )) +
        geom_col(col = "black") +
        geom_text(aes(label = value),
          position = position_stack(vjust = 0.5)
        ) +
        coord_polar(theta = "y") +
        scale_fill_brewer(palette = "GnBu") +
        xlim(c(0.2, hsize + 0.5)) +
        theme(
          panel.background = element_rect(fill = "white"),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()
        )
    },
    #' @description
    #' pie plot
    #' @param data data
    #' @param y map to y-axis
    #' @param fill map to fill
    #' @param ... additional aesthetics properties mapping
    pie = function(data, y, fill) {
      ggplot(
        data,
        aes(
          x = "", y = {{ y }},
          fill = fct_inorder({{ fill }})
        )
      ) +
        geom_col(width = 1, col = 1) +
        geom_text(aes(label = value),
          position = position_stack(vjust = 0.5)
        ) +
        coord_polar(theta = "y") +
        guides(fill = guide_legend(title = "Group")) +
        scale_y_continuous(breaks = df2$pos, labels = df$group) +
        theme(
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 15),
          legend.position = "none", # Removes the legend
          panel.background = element_rect(fill = "white")
        )
    }
  )
)
