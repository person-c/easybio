#' @title Visualization Artist for Custom Plots
#'
#' @description
#' The `Artist` class provides a collection of methods to create various types of plots using `ggplot2`.
#' It includes methods for generating dumbbell plots, bubble plots, divergence bar charts, lollipop plots,
#' contour plots, scatter plots with ellipses, donut plots, and pie charts. Each method is designed
#' to map data to specific aesthetics and apply additional customizations.
#'
#' @import ggplot2
#' @export
Artist <- R6::R6Class("Artist",
  public = list(
    #' @description
    #' Create a dumbbell plot
    #'
    #' This method generates a dumbbell plot using the provided data, mapping the specified columns
    #' to the x-axis, y-axis, and color aesthetic.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param x The column in `data` to map to the x-axis.
    #' @param y The column in `data` to map to the y-axis.
    #' @param col The column in `data` to map to the color aesthetic.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return A ggplot object representing the dumbbell plot.
    dumbbbell = function(data, x, y, col, ...) {
      ggplot(data, aes(x = {{ x }}, y = {{ y }}), ...) +
        geom_line() +
        geom_point(aes(col = {{ col }}), size = 3)
    },

    #' @description
    #' Create a bubble plot
    #'
    #' This method generates a bubble plot where points are mapped to the x and y axes, with their
    #' size and color representing additional variables.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param x The column in `data` to map to the x-axis.
    #' @param y The column in `data` to map to the y-axis.
    #' @param size The column in `data` to map to the size of the points.
    #' @param col The column in `data` to map to the color of the points.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return A ggplot object representing the bubble plot.
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
    #' Create a divergence bar chart
    #'
    #' This method generates a divergence bar chart where bars are colored based on their
    #' positive or negative value.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param group The column in `data` representing the grouping variable.
    #' @param y The column in `data` to map to the y-axis.
    #' @param fill The column in `data` to map to the fill color of the bars.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return A ggplot object representing the divergence bar chart.
    barchart_divergence = function(data, group, y, fill, ...) {
      ggplot(
        data,
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
          label = {{ group }},
          hjust = ifelse(y < 0, 1.5, -1),
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
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank()
        )
    },

    #' @description
    #' Create a lollipop plot
    #'
    #' This method generates a lollipop plot, where points are connected to a baseline by vertical
    #' segments, with customizable colors and labels.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param x The column in `data` to map to the x-axis.
    #' @param y The column in `data` to map to the y-axis.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return A ggplot object representing the lollipop plot.
    lollipop = function(data, x, y, ...) {
      ggplot(data, aes(x = {{ x }}, y = {{ y }}, ...)) +
        geom_segment(aes(x = {{ x }}, xend = {{ x }}, y = 0, yend = {{ y }}),
          col = "gray", lwd = 1
        ) +
        geom_point(size = 7.5, pch = 21, bg = 4, col = 1) +
        geom_text(aes(label = {{ y }}), col = "white", size = 3) +
        scale_x_discrete(labels = paste0("G_", 1:10)) +
        coord_flip() +
        theme_minimal()
    },

    #' @description
    #' Create a contour plot
    #'
    #' This method generates a contour plot that includes filled and outlined density contours,
    #' with data points overlaid.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param x The column in `data` to map to the x-axis.
    #' @param y The column in `data` to map to the y-axis.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return A ggplot object representing the contour plot.
    contour = function(data, x, y, ...) {
      ggplot(data, aes(x = {{ x }}, y = {{ y }}, ...)) +
        geom_point() +
        geom_density_2d_filled(alpha = 0.4) +
        geom_density_2d(colour = "black")
    },

    #' @description
    #' Create a scatter plot with ellipses
    #'
    #' This method generates a scatter plot where data points are colored by group, with ellipses
    #' representing the confidence intervals for each group.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param x The column in `data` to map to the x-axis.
    #' @param y The column in `data` to map to the y-axis.
    #' @param col The column in `data` to map to the color aesthetic.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return A ggplot object representing the scatter plot with ellipses.
    scatter_ellipses = function(data, x, y, col, ...) {
      ggplot(data, aes(
        x = {{ x }},
        y = {{ y }}, col = {{ col }}, ...
      )) +
        geom_point() +
        stat_ellipse(
          geom = "polygon",
          aes(fill = {{ col }}),
          alpha = 0.25
        )
    },

    #' @description
    #' Create a donut plot
    #'
    #' This method generates a donut plot, which is a variation of a pie chart with a hole in the center.
    #' The sections of the donut represent the proportion of categories in the data.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param x The column in `data` to map to the x-axis.
    #' @param y The column in `data` to map to the y-axis.
    #' @param fill The column in `data` to map to the fill color of the sections.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return A ggplot object representing the donut plot.
    donut = function(data, x, y, fill, ...) {
      hsize <- 3
      ggplot(data, aes(
        x = {{ x }}, y = {{ y }},
        fill = {{ fill }}, ...
      )) +
        geom_col(col = "black") +
        geom_text(aes(label = {{ y }}),
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
    #' Create a pie chart
    #'
    #' This method generates a pie chart where sections represent the proportion of categories in the data.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param y The column in `data` to map to the y-axis.
    #' @param fill The column in `data` to map to the fill color of the sections.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return A ggplot object representing the pie chart.
    pie = function(data, y, fill, ...) {
      ggplot(
        data,
        aes(
          x = "", y = {{ y }},
          fill = fct_inorder({{ fill }})
        )
      ) +
        geom_col(width = 1, col = 1) +
        geom_text(aes(label = {{ y }}),
          position = position_stack(vjust = 0.5)
        ) +
        coord_polar(theta = "y") +
        guides(fill = guide_legend(title = "Group")) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
        theme(
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 15),
          legend.position = "none",
          panel.background = element_rect(fill = "white")
        )
    }
  )
)
