#' @title Reductor Class
#'
#' @description
#' This is the abstract base class for reduction task including PCA,
#' UMAP and t-SNE.
#' @import ggplot2
#' @export
Artist <- R6::R6Class("Artist",
  public = list(
    #' @description
    #' boxplot
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param fill map to fill
    #' @param ... additional aesthetics properties mapping
    box = function(data, x, y, fill, ...) {
      cols <- c("#CFD8DC", "#90A4AE", "#455A64")

      ggplot(
        data,
        aes(
          x = {{ x }}, y = {{ y }},
          fill = {{ fill }}, ...
        )
      ) +
        geom_boxplot(
          alpha = 0.8, # Fill transparency
          colour = "#474747", # Border col
          outlier.colour = 1
        ) + # Outlier col
        scale_fill_manual(values = cols) # Fill cols
    },
    #' @description
    #' density area plot
    #' @param data data
    #' @param x map to x-axis
    #' @param fill map to fill
    #' @param ... additional aesthetics properties mapping
    density = function(data, x, fill, ...) {
      cols <- c("#F76D5E", "#FFFFBF", "#72D8FF")
      ggplot(data, aes(x = {{ x }}, fill = {{ fill }}, ...)) +
        geom_density(alpha = 0.8, col = NA) +
        scale_fill_manual(values = cols)
    },
    #' @description
    #' violin plot
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param ... additional aesthetics properties mapping
    violin = function(data, x, y, ...) {
      ggplot(
        data,
        aes(
          x = {{ x }}, y = {{ y }},
          fill = {{ x }}, ...
        )
      ) +
        geom_violin(trim = FALSE) +
        geom_boxplot(width = 0.07) +
        scale_fill_manual(values = c("#BCE4D8", "#49A4B9", "#2C5985"))
    },
    #' @description
    #' scatter plot
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param col map to col
    #' @param size map to size
    #' @param ... additional aesthetics properties mapping
    scatter = function(data, x, y, col = NULL, size = NULL, ...) {
      ggplot(data = data, aes(x = {{ x }}, y = {{ y }}, ...)) +
        geom_point(aes(col = {{ col }}, size = {{ size }})) +
        geom_smooth(method = "loess", se = FALSE)
    },
    #' @description
    #' counts plot
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param ... additional aesthetics properties mapping
    counts = function(data, x, y, ...) {
      ggplot(data, aes(x = {{ x }}, y = {{ y }}, ...)) +
        geom_count(col = "tomato3", show.legend = FALSE)
    },
    #' @description
    #' correlation plot
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param fill map to fill
    #' @param ... additional aesthetics properties mapping
    cor = function(data, x, y, fill, ...) {
      ggplot2::ggplot(
        data,
        ggplot2::aes(
          x = {{ x }}, y = {{ y }},
          fill = {{ fill }}, ...
        )
      ) +
        ggplot2::geom_raster() +
        ggplot2::scale_x_discrete(expand = ggplot2::expansion(0)) +
        ggplot2::scale_y_discrete(expand = ggplot2::expansion(0)) +
        ggplot2::scale_fill_gradient2(
          midpoint = 0,
          limit = c(-1, 1)
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.title = ggplot2::element_blank(),
          axis.textx = ggplot2::element_text(angle = 90)
        )
    },
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
    #' heatmap
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param fill map to fill
    #' @param ... additional aesthetics properties mapping
    heatmap = function(data, x, y, fill, ...) {
      ggplot2::ggplot(
        data,
        ggplot2::aes(
          x = {{ x }}, y = {{ y }},
          fill = {{ fill }}, ...
        )
      ) +
        ggplot2::geom_raster() +
        ggplot2::scale_x_discrete(expand = ggplot2::expansion(0)) +
        ggplot2::scale_y_discrete(expand = ggplot2::expansion(0)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.title = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(angle = 90)
        )
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
    #' sankey plot
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param next_x next_x
    #' @param node node
    #' @param next_node node
    #' @param ... additional aesthetics properties mapping
    sankey = function(data, x, y, next_x, node, next_node, ...) {
      ggplot(data, aes(
        x = {{ x }},
        next_x = {{ next_x }},
        node = {{ node }},
        next_node = {{ next_node }},
        fill = factor({{ node }}),
        label = {{ node }},
        ...
      )) +
        ggsankey::geom_sankey(flow.alpha = 0.5, nodecol = 1) +
        ggsankey::geom_sankey_label(size = 3.5, col = 1, fill = "white") +
        scale_fill_viridis_d() +
        theme_sankey(base_size = 16) +
        theme(legend.position = "none")
    },
    #' @description
    #' waterfull plot
    #' @param data data
    waterfall = function(data) {
      waterfalls::waterfall(data)
    },
    #' @description
    #' ridges plot
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param ... additional aesthetics properties mapping
    ridges = function(data, x, y, ...) {
      ggplot(data, aes({{ x }},
        y = {{ y }},
        fill = 0.5 - abs(0.5 - stat(ecdf)),
        ...
      )) +
        ggridges::stat_density_ridges(
          geom = "density_ridges_gradient", calc_ecdf = TRUE
        ) +
        scale_fill_gradient(
          low = "white", high = "#87CEFF",
          name = "Tail prob."
        )
    },
    #' @description
    #' stream plot
    #' @param data data
    #' @param x map to x-axis
    #' @param y map to y-axis
    #' @param fill map to fill
    #' @param ... additional aesthetics properties mapping
    stream = function(data, x, y, fill, ...) {
      cols <- c("#FFB400", "#FFC740", "#C20008", "#FF020D", "#13AFEF")

      ggplot(
        data,
        aes(
          x = {{ x }},
          y = {{ y }}, fill = {{ fill }}, ...
        )
      ) +
        ggstream::geom_stream(extra_span = 0.2) +
        ggstream::geom_stream(
          extra_span = 0.2, true_range = "none",
          alpha = 0.3
        ) +
        scale_fill_manual(values = cols) +
        theme_minimal()
    },
    #' @description
    #' slope plot
    #' @param data data
    #' @param times times
    #' @param measure measure
    #' @param group group
    slope = function(data, times, measure, group) {
      CGPfunctions::newggslopegraph(data,
        times = {{ times }},
        measure = {{ measure }},
        group = {{ group }},
        Title = "GDP evolution",
        SubTitle = "1970-1979",
        Caption = "By R CHARTS",
        ThemeChoice = "tufte"
      )
    },
    #' @description
    #' treemap
    #' @param data data
    #' @param area area
    #' @param lable map to lable
    #' @param fill map to fill
    #' @param subgroup map to subgroup
    #' @param ... additional aesthetics properties mapping
    treemap = function(data, area, fill, lable, subgroup, ...) {
      ggplot(data, aes(
        area = {{ area }}, fill = {{ fill }},
        label = {{ lable }}, subgroup = {{ subgroup }},
        ...
      )) +
        treemapify::geom_treemap() +
        treemapify::geom_treemap_subgroup_border(colour = "white", size = 5) +
        treemapify::geom_treemap_subgroup_text(
          place = "centre", grow = TRUE,
          alpha = 0.25, colour = "black",
          fontface = "italic"
        ) +
        treemapify::geom_treemap_text(
          colour = "white", place = "centre",
          size = 15, grow = TRUE
        )
    },
    #' @description
    #' alluvila plot
    #' @param data data
    #' @param axis1 map to axis
    #' @param axis2 map to axis2
    #' @param y map to y
    #' @param fill map to fill
    #' @param ... additional aesthetics properties mapping
    alluvial = function(data, axis1, axis2, y, fill, ...) {
      ggplot(
        data = {{ data }},
        aes(
          axis1 = {{ axis1 }}, axis2 = {{ axis2 }},
          y = {{ y }}, ...
        )
      ) +
        ggalluvial::geom_alluvium(aes(fill = {{ fill }})) +
        ggalluvial::geom_stratum() +
        geom_text(
          stat = "stratum",
          aes(label = after_stat(stratum))
        ) +
        scale_x_discrete(
          limits = c("Survey", "Response"),
          expand = c(0.15, 0.05)
        ) +
        scale_fill_viridis_d() +
        theme_void()
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
    #' venn plot
    #' @param data data
    venn = function(data, ...) {
      ggvenn::ggvenn(data, ...)
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
