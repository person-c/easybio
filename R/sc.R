#' Insert character value in specified index  for character in a weird way(Used only for single cell annotation).
#'
#' @param x expression
#' @param na the default value for the NA.
#'
#' @return character
#' @export
#' @examples
#' library(easybio)
#' tmp <- finsert(
#'   x = expression(
#'     c(0, 1, 3) == "Neutrophil",
#'     c(2, 4, 8) == "Macrophage"
#'   ),
#'   na = "Unknown"
#' )
#' names(tmp) <- as.character(0:8)
finsert <- function(x = expression(c(0, 1, 3) == "Neutrophil", c(2, 4, 8) == "Macrophage"), na = "Unknown") {
  x <- eval(substitute(x))
  x <- lapply(x, as.list)
  x <- rapply(x, eval, classes = "call", how = "replace")
  x <- unlist(x, recursive = FALSE)
  itor <- 1
  v <- character()
  while (itor < length(x)) {
    v[x[[itor + 1]] + 1] <- x[[itor + 2]]
    itor <- itor + 3
  }
  v[is.na(v)] <- na
  v
}

#' match markers from cellMarker2
#'
#' @param marker markers from `FindAllMarkers`
#' @param n top number of genes to match
#' @param spc 'Human' or 'Mouse'
#'
#' @return matched cellMarker data.frame
#' @import data.table
#' @export
#' @examples
#' library(easybio)
#' data(exampleMarker)
#' matchCellMarker2(exampleMarker, n = 50, spc = "Human")
matchCellMarker2 <- function(marker, n, spc) {
  . <- markerWith <- NULL
  species <- avg_log2FC <- p_val_adj <- cluster <- gene <- cell_name <- count <- N <- NULL

  marker <- copy(marker)
  setDT(marker)

  cellMarker2 <- cellMarker2[.(spc), on = .(species)]
  marker <- marker[
    avg_log2FC > 0 & p_val_adj < 0.05,
    .SD[order(-avg_log2FC)][1:n],
    keyby = .(cluster)
  ]

  marker <- marker[cellMarker2, on = "gene==marker", nomatch = NULL]
  marker <- marker[, .(markerWith = .(gene), N = .N), by = .(cluster, cell_name)]
  marker <- marker[N > 0, .SD[order(-N)], keyby = .(cluster)]

  topcell <- marker[, head(.SD, 2), keyby = .(cluster)][, unique(cell_name)]
  topmarker <- cellMarker2[.(topcell), .SD, on = .(cell_name)]
  topmarker <- topmarker[!is.na(marker), .(count = .N), by = .(cell_name, marker)]
  topmarker <- topmarker[, head(.SD[order(-count)], 10), by = .(cell_name)]
  setnames(topmarker, old = "marker", new = "Symbol")
  setattr(marker, "topmarker", topmarker)

  marker[, let(uniqueN = sapply(markerWith, FUN = \(x) length(unique(x))))]
  marker[, let(ordered_symbol = lapply(markerWith, FUN = \(x) names(sort(unclass(table(x)), TRUE))))]
  marker[, let(orderN = lapply(markerWith, \(x) as.integer(sort(unclass(table(x)), TRUE))))]
  setcolorder(marker, c("cluster", "cell_name", "uniqueN", "N", "ordered_symbol", "orderN", "markerWith"))
  marker[]
}

#' Check markers from cell returned by `matchCellMarker`
#'
#' @param marker markers from `Seurat::FindAllMarkers`
#' @param n top number of genes to match.
#' @param spc 'Human' or 'Mouse'.
#' @param cl cluster you want to check.
#' @param unique not yet implemented...waiting for.
#'
#' @return matched cellMarker data.frame
#' @import data.table
#' @export
#' @examples
#' library(easybio)
#' data(exampleMarker)
#' check_marker(exampleMarker, n = 50, spc = "Human", cl = c(1, 4, 7))
check_marker <- function(marker, n, spc, cl = c(), unique = TRUE) {
  . <- NULL
  Symbol <- cell_name <- cluster <- count <- NULL

  marker_matched <- matchCellMarker2(marker, n, spc)
  marker_top <- na.omit(attributes(marker_matched)[["topmarker"]])
  marker_top <- unique(marker_top, by = "Symbol")
  marker_top <- marker_top[
    count > 1, .(Symbol = .(Symbol)),
    by = .(cell_name)
  ]

  cell_on_cluster <- marker_matched[.(factor(cl)), head(.SD, 5), by = .(cluster)]
  print(cell_on_cluster)

  cell <- marker_matched[.(factor(cl)), head(.SD, 2), by = .(cluster)][, unique(cell_name)]
  marker_top <- marker_top[.(cell), .SD, on = .(cell_name)]
  cell_features <- setNames(marker_top$Symbol, make.names(marker_top$cell_name))
  cell_features
}


#' Get marker for specific cell in cellMarker2.
#'
#' @param spc 'Human' or 'Mouse'.
#' @param cell cell you want to search.
#' @param number number ot marker you want to get.
#' @param min.count the minimum number the marker reporteded
#'
#' @return marker list
#' @import data.table
#' @export
#' @examples
#' library(easybio)
#' get_marker(spc = "Human", cell = c("Macrophage", "Monocyte"))
get_marker <- function(spc, cell = character(), number = 5, min.count = 1) {
  . <- NULL
  species <- cell_name <- N <- marker <- NULL

  marker <- cellMarker2[.(spc, cell),
    .SD,
    on = .(species, cell_name)
  ][, .N, by = .(cell_name, marker)][
    N > min.count, na.omit(.SD)[order(-N)] |> head(number),
    by = .(cell_name)
  ][, .(marker = .(marker)), by = .(cell_name)]
  marker <- setNames(marker$marker, marker$cell_name)

  marker
}

#' Check dotplot for markers from \code{check_marker}.
#'
#' @param srt Seurat object.
#' @param marker From \code{Seurat::FindAllMarkers}
#' @param n params for \code{check_marker}
#' @param species params for \code{check_marker}
#' @param cls a list contain cluster groups to check
#'
#' @return NULL
#' @import ggplot2
#' @export
plotSeuratDot <- function(srt, marker, n = 50, species, cls) {
  tmpdir <- "tmp_check_marker"
  if (!dir.exists(tmpdir)) {
    message("Creating tmp_check_marker directory...please check all result in ths directory")
    dir.create(tmpdir)
  }

  lapply(cls, \(cl) {
    check_marker(marker, n, species, cl) |>
      Seurat::DotPlot(srt, features = _) +
      scale_x_discrete(
        guide = guide_axis(
          angle = 60,
          theme = theme(text = element_text(size = 4))
        )
      ) +
      theme(
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, angle = 30, vjust = 0.1, hjust = 0)
      )


    ggsave(
      width = 10, height = 8,
      file.path(
        tmpdir,
        paste0(paste0("clusters_", paste0(cl, collapse = "_")), ".pdf")
      )
    )
  })
}
