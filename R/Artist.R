#' @title Visualization Artist for Custom Plots
#'
#' @description
#' The `Artist` class offers a suite of methods designed to create a variety of plots using `ggplot2` for
#' data exploration. All methods log their calls and results, allowing you to review all outcomes later
#' via the `get_all_results()` method.
#'
#' Each `plot_*` method displays the generating command as the plot title. When a `plot_*` method
#' immediately follows a `test_*` method, the test result (p-value) is automatically added as the subtitle.
#'
#' All methods return `invisible(self)`, enabling fluent method chaining.
#'
#' @import ggplot2 R6 data.table
#' @return The `R6` class [Artist].
#' @export
#'
#' @examples
#' library(data.table)
#' air <- subset(airquality, Month %in% c(5, 6))
#' setDT(air)
#' cying <- Artist$new(data = air)
#' cying$plot_scatter(x = Wind, y = Temp)
#' cying$test_wilcox(formula = Ozone ~ Month)
#' cying$plot_scatter(x = Wind, y = Temp)
#'
Artist <- R6::R6Class("Artist",
  public = list(
    #' @field data Stores the dataset used for plotting.
    data = NULL,
    #' @field command recode the command.
    command = list(),
    #' @field result record the plot.
    result = list(),
    #' @description
    #' Initializes the Artist class with an optional dataset.
    #'
    #' @param data A data frame  containing the dataset to be used for plotting. Default is `NULL`.
    #' @return An instance of the Artist class.
    initialize = function(data = NULL) {
      message("Welcome to the amazing Artist's world, enjoy exploring your data in a new way!")
      self$data <- data
    },
    #' @description
    #' Get all history result
    #'
    #' @return a data.table object
    get_all_result = function() {
      data.table(command = self$command, result = self$result)
    },
    #' @description
    #' Conduct wilcox.test
    #'
    #' @param formula [wilcox.test()] formula arguments
    #' @param data A data frame containing the data to be plotted. Default is `self$data`.
    #' @param ... Additional aesthetic mappings passed to [wilcox.test()].
    #' @return The Artist object invisibly.
    test_wilcox = function(formula, data = self$data, ...) {
      mc <- match.call()
      htestRes <- wilcox.test(formula = formula, data = data, ...)
      private$finalize_test(htestRes, mc)
    },
    #' @description
    #' Conduct t.test
    #'
    #' @param formula [t.test()] formula arguments
    #' @param data A data frame containing the data to be plotted. Default is `self$data`.
    #' @param ... Additional aesthetic mappings passed to [t.test()].
    #' @return The Artist object invisibly.
    test_t = function(formula, data = self$data, ...) {
      mc <- match.call()
      htestRes <- t.test(formula = formula, data = data, ...)
      private$finalize_test(htestRes, mc)
    },
    #' @description
    #' Creates a scatter plot.
    #'
    #' @param data A data frame containing the data to be plotted. Default is `self$data`.
    #' @param x The column name for the x-axis.
    #' @param y The column name for the y-axis.
    #' @param add whether to add the test result as subtitle.
    #' @param fun function to process the `self$data`.
    #' @param ... Additional aesthetic mappings passed to `aes()`.
    #' @return The Artist object invisibly.
    plot_scatter = function(data = self$data, fun = \(x) x, x, y, ..., add = private$is_htest()) {
      mc <- match.call()
      data <- force(fun)(data)

      p <- ggplot(data, aes(x = {{ x }}, y = {{ y }}, ...)) +
        geom_point()

      private$finalize_plot(p, mc, add)
    },
    #' @description
    #' Creates a box plot.
    #'
    #' @param data A data frame or tibble containing the data to be plotted. Default is `self$data`.
    #' @param x The column name for the x-axis.
    #' @param add whether to add the test result as subtitle.
    #' @param fun function to process the `self$data`.
    #' @param ... Additional aesthetic mappings passed to `aes()`.
    #' @return The Artist object invisibly.
    plot_box = function(data = self$data, fun = \(x) x, x, ..., add = private$is_htest()) {
      mc <- match.call()
      data <- force(fun)(data)

      p <- ggplot(data, aes(x = {{ x }}, ...)) +
        geom_boxplot()

      private$finalize_plot(p, mc, add)
    },
    #' @description
    #' Creates a dumbbell plot.
    #'
    #' This method generates a dumbbell plot using the provided data, mapping the specified columns
    #' to the x-axis, y-axis, and color aesthetic.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param x The column in `data` to map to the x-axis.
    #' @param y The column in `data` to map to the y-axis.
    #' @param col The column in `data` to map to the color aesthetic.
    #' @param add whether to add the test result as subtitle.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return The Artist object invisibly.
    plot_dumbbell = function(data = self$data, x, y, col, add = private$is_htest(), ...) {
      mc <- match.call()
      p <- ggplot(data, aes(x = {{ x }}, y = {{ y }}), ...) +
        geom_line() +
        geom_point(aes(col = {{ col }}), size = 3)

      private$finalize_plot(p, mc, add)
    },

    #' @description
    #' Creates a bubble plot.
    #'
    #' This method generates a bubble plot where points are mapped to the x and y axes, with their
    #' size and color representing additional variables.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param x The column in `data` to map to the x-axis.
    #' @param y The column in `data` to map to the y-axis.
    #' @param size The column in `data` to map to the size of the points.
    #' @param col The column in `data` to map to the color of the points.
    #' @param add whether to add the test result as subtitle.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return The Artist object invisibly.
    plot_bubble = function(data = self$data, x, y, size, col, add = private$is_htest(), ...) {
      mc <- match.call()
      p <- ggplot(
        data,
        aes(
          x = {{ x }}, y = {{ y }},
          size = {{ size }}, col = {{ col }},
          ...
        )
      ) +
        geom_point() +
        scale_size(name = "Size", range = c(1, 10))

      private$finalize_plot(p, mc, add)
    },

    #' @description
    #' Creates a divergence bar chart.
    #'
    #' This method generates a divergence bar chart where bars are colored based on their
    #' positive or negative value.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param group The column in `data` representing the grouping variable.
    #' @param y The column in `data` to map to the y-axis.
    #' @param add whether to add the test result as subtitle.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return The Artist object invisibly.
    plot_barchart_divergence = function(data = self$data, group, y, add = private$is_htest(), ...) {
      mc <- match.call()
      y_vec <- data[[deparse(substitute(y))]]
      p <- ggplot(
        data,
        aes(
          x = reorder({{ group }}, {{ y }}),
          y = {{ y }}, ...
        )
      ) +
        geom_bar(
          stat = "identity",
          show.legend = FALSE,
          fill = ifelse(y_vec >= 0, "lightblue", "lightpink"),
          col = "white"
        ) +
        geom_hline(yintercept = 0, col = 1, lwd = 0.2) +
        geom_text(aes(
          label = {{ group }},
          hjust = ifelse({{ y }} < 0, 1.5, -1),
          vjust = 0.5
        ), size = 2.5) +
        xlab("Group") +
        ylab("Value") +
        coord_flip() +
        theme_minimal() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank()
        )

      private$finalize_plot(p, mc, add)
    },

    #' @description
    #' Creates a lollipop plot.
    #'
    #' This method generates a lollipop plot, where points are connected to a baseline by vertical
    #' segments.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param x The column in `data` to map to the x-axis.
    #' @param y The column in `data` to map to the y-axis.
    #' @param add whether to add the test result as subtitle.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return The Artist object invisibly.
    plot_lollipop = function(data = self$data, x, y, add = private$is_htest(), ...) {
      mc <- match.call()
      p <- ggplot(data, aes(x = {{ x }}, y = {{ y }}, ...)) +
        geom_segment(aes(x = {{ x }}, xend = {{ x }}, y = 0, yend = {{ y }}),
          col = "gray", lwd = 1
        ) +
        geom_point(size = 7.5, pch = 21, bg = 4, col = 1) +
        geom_text(aes(label = {{ y }}), col = "white", size = 3) +
        coord_flip() +
        theme_minimal()

      private$finalize_plot(p, mc, add)
    },

    #' @description
    #' Creates a contour plot.
    #'
    #' This method generates a contour plot that includes filled and outlined density contours,
    #' with data points overlaid.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param x The column in `data` to map to the x-axis.
    #' @param y The column in `data` to map to the y-axis.
    #' @param add whether to add the test result as subtitle.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return The Artist object invisibly.
    plot_contour = function(data = self$data, x, y, add = private$is_htest(), ...) {
      mc <- match.call()
      p <- ggplot(data, aes(x = {{ x }}, y = {{ y }}, ...)) +
        geom_point() +
        geom_density_2d_filled(alpha = 0.4) +
        geom_density_2d(colour = "black")

      private$finalize_plot(p, mc, add)
    },

    #' @description
    #' Creates a scatter plot with ellipses.
    #'
    #' This method generates a scatter plot where data points are colored by group, with ellipses
    #' representing the confidence intervals for each group.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param x The column in `data` to map to the x-axis.
    #' @param y The column in `data` to map to the y-axis.
    #' @param col The column in `data` to map to the color aesthetic.
    #' @param add whether to add the test result as subtitle.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return The Artist object invisibly.
    plot_scatter_ellipses = function(data = self$data, x, y, col, add = private$is_htest(), ...) {
      mc <- match.call()
      p <- ggplot(data, aes(
        x = {{ x }},
        y = {{ y }}, col = {{ col }}, ...
      )) +
        geom_point() +
        stat_ellipse(
          geom = "polygon",
          aes(fill = {{ col }}),
          alpha = 0.25
        )

      private$finalize_plot(p, mc, add)
    },

    #' @description
    #' Creates a donut plot.
    #'
    #' This method generates a donut plot, which is a variation of a pie chart with a hole in the center.
    #' The sections of the donut represent the proportion of categories in the data.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param x The column in `data` to map to the x-axis.
    #' @param y The column in `data` to map to the y-axis.
    #' @param fill The column in `data` to map to the fill color of the sections.
    #' @param add whether to add the test result as subtitle.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return The Artist object invisibly.
    plot_donut = function(data = self$data, x, y, fill, add = private$is_htest(), ...) {
      mc <- match.call()
      hsize <- 3
      p <- ggplot(data, aes(
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

      private$finalize_plot(p, mc, add)
    },

    #' @description
    #' Creates a pie chart.
    #'
    #' This method generates a pie chart where sections represent the proportion of categories in the data.
    #'
    #' @param data A data frame containing the data to be plotted.
    #' @param y The column in `data` to map to the y-axis.
    #' @param fill The column in `data` to map to the fill color of the sections.
    #' @param add whether to add the test result as subtitle.
    #' @param ... Additional aesthetic mappings or other arguments passed to `ggplot`.
    #'
    #' @return The Artist object invisibly.
    plot_pie = function(data = self$data, y, fill, add = private$is_htest(), ...) {
      mc <- match.call()
      p <- ggplot(
        data,
        aes(
          x = "", y = {{ y }},
          fill = factor({{ fill }}, levels = unique({{ fill }}))
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

      private$finalize_plot(p, mc, add)
    }
  ),
  private = list(
    finalize_plot = function(p, mc, add_htest = FALSE) {
      cmd_text <- paste(strwrap(deparse(mc, width.cutoff = 500L), width = 60L),
        collapse = "\n"
      )
      p <- p + labs(title = cmd_text)

      if (add_htest && length(self$result) > 0L) {
        prev <- self$result[[length(self$result)]]
        if (inherits(prev, "htest")) {
          p <- p + labs(subtitle = sprintf("p = %.4g", prev$p.value))
        }
      }

      print(p)
      self$command <- private$add_in_list(self$command, mc)
      self$result <- private$add_in_list(self$result, p)
      invisible(self)
    },
    finalize_test = function(htestRes, mc) {
      self$command <- private$add_in_list(self$command, mc)
      self$result <- private$add_in_list(self$result, htestRes)
      invisible(self)
    },
    last = function(x) {
      x[[length(x)]]
    },
    add_in_list = function(x = list(), element) {
      x[[length(x) + 1]] <- element
      x
    },
    is_htest = function(x = self$result) {
      if (length(self$result) == 0L) {
        return(FALSE)
      }
      inherits(private$last(x), "htest")
    }
  )
)
