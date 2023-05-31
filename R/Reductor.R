#' @title Reductor Class
#'
#' @description
#' This is the abstract base class for reduction task including PCA,
#' UMAP and t-SNE.
#'
#' @export
#' @examples
#' # reductor <- Reductor$new('tsne')
#' # tune_fit <- reductor$tune(data,
#' # perplexity = c(30, 40, 50, 60),
#' # n_iter = c(1000, 2000, 2500))
#' # reductor$plot(tune_fit)
Reductor <- R6::R6Class("Reductor",
  public = list(
    #' @field method
    #'  which method you choose to do the reduction
    method = NA,

    #' @description
    #' construct a new instance to do the reduction
    #' @param method the method you want to do the reduction
    initialize = function(method) {
      self$method <- method
    },

    #' @description a wrapper for auto-tuning for the reduction
    #' @importFrom data.table CJ
    #' @param input the data you want to do the reduction
    #' @param ... the args space
    #' @param workers number of cores used to do parellel computation
    tune = function(input, ..., workers = parallel::detectCores() - 4) {
      force(input)
      arg_space <- CJ(...) |> asplit(1)

      future::plan(future::multisession, workers = workers)

      results <- furrr::future_map(
        .x = arg_space,
        .f = ~ {
          params <- as.list(.x)
          result <- tryCatch(
            {
              result <- rlang::exec(self$reduction, data = input, !!!params)
              return(list(result = result, params = params))
            },
            error = function(e) {
              logger::log_error("Error occurred:", conditionMessage(e))
              return(NULL)
            }
          )
          if (!is.null(result)) {
            return(result)
          }
        }
      )
      return(results)
    }
  ),
  active = list(
    #' @field reduction a function used to execute the reduction
    reduction = function() {
      private[[self$method]]
    },

    #' @field plot a function used to plot the reduction results
    plot = function() {
      private[[paste0("plot_", self$method)]] #
    }
  ),
  private = list(
    pca = function(data, ...) {
      pc <- stats::prcomp(data,
        center = TRUE,
        scale. = TRUE,
        ...
      )
    },
    plot_pca = function(pc, group) {
      g <- ggbiplot::ggbiplot(pc,
        obs.scale = 1,
        var.scale = 1,
        groups = group,
        ellipse = TRUE,
        circle = TRUE,
        ellipse.prob = 0.68
      )
      g <- g + ggplot2::scale_color_discrete(name = "")
      g <- g + ggplot2::theme(
        legend.direction = "horizontal",
        legend.position = "top"
      )
      print(g)
    },
    tsne = function(data, ...) {
      tsne_fit <- Rtsne::Rtsne(data, ...)
    },
    plot_tsne = function(tsne, group) {
      if (!is.list(tsne[[1]])) {
        data <- purrr::map(list(tsne$Y, group), as.data.frame)
        data <- rlang::exec("cbind", !!!data)
        colnames(data)[1:3] <- c("tSNE1", "tSNE2", "Group")
        ggplot2::ggplot(
          data,
          ggplot2::aes(
            x = tSNE1,
            y = tSNE2,
            color = Group
          )
        ) +
          ggplot2::geom_point() +
          ggplot2::theme(
            legend.position = "bottom",
            legend.title = ggplot2::element_blank()
          )
      } else {
        Y <- purrr::map(
          .x = tsne,
          .f = ~ purrr::pluck(.x, "result", "Y")
        )
        Y <- purrr::map(Y, ~ as.data.frame(.x))
        Y <- purrr::map(Y, ~ {
          `<-`(names(.x), c("tSNE1", "tSNE2"))
          .x
        })

        Y <- purrr::map2(
          list(group <- (as.data.frame(group) |> setNames("Group"))),
          Y, cbind
        )
        purrr::map(
          Y,
          ~ {
            ggplot2::ggplot(
              .x,
              ggplot2::aes(
                x = tSNE1,
                y = tSNE2,
                color = Group
              )
            ) +
              ggplot2::geom_point() +
              ggplot2::theme(
                legend.position = "bottom",
                legend.title = ggplot2::element_blank()
              )
          }
        )
      }
    },
    umap = function(data, ...) {
      umap::umap(data, ...)
    },
    plot_umap = function(umap, group) {
      if (!is.list(umap[[1]])) {
        data <- purrr::map(list(umap$layout, group), as.data.frame)
        data <- rlang::exec("cbind", !!!data)
        colnames(data)[1:3] <- c("umap1", "umap2", "Group")
        ggplot2::ggplot(
          data,
          ggplot2::aes(
            x = umap1,
            y = umap2,
            color = Group
          )
        ) +
          ggplot2::geom_point() +
          ggplot2::theme(
            legend.position = "bottom",
            legend.title = ggplot2::element_blank()
          )
      } else {
        layouts <- purrr::map(.x = umap, .f = ~ purrr::pluck(.x, "result", "layout"))
        layouts <- purrr::map(layouts, ~ as.data.frame(.x))
        layouts <- purrr::map(layouts, ~ {
          `<-`(names(.x), c("umap1", "umap2"))
          .x
        })

        layouts <- purrr::map2(list(group <- (as.data.frame(group) |> setNames("Group"))), layouts, cbind)
        purrr::map(
          layouts,
          ~ {
            ggplot2::ggplot(
              .x,
              ggplot2::aes(
                x = umap1,
                y = umap2,
                color = Group
              )
            ) +
              ggplot2::geom_point() +
              ggplot2::theme(
                legend.position = "bottom",
                legend.title = ggplot2::element_blank()
              )
          }
        )
      }
    }
  )
)
