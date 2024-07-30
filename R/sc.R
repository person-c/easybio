#' Insert character value in specified index  for character in a weird way(Used only for single cell annotation).
#'
#' @param x expression
#' @param na the default value for the NA.
#'
#' @return character
#' @export
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
matchCellMarker2 <- function(marker, n, spc) {
  . <- markerWith <- symbol <- NULL
  species <- avg_log2FC <- p_val_adj <- cluster <- gene <- cell_name <- count <- N <- len <- NULL

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

  marker[, len := sapply(markerWith, FUN = \(x) length(unique(x)))]
  marker[, symbol := lapply(markerWith, FUN = \(x) unique(x))]
  setcolorder(marker, c("cluster", "cell_name", "len", "N", "symbol", "markerWith"))
  marker[]
}

#' Check markers from cell returned by `matchCellMarker`
#'
#' @param marker markers from `Seurat::FindAllMarkers`
#' @param n top number of genes to match.
#' @param spc 'Human' or 'Mouse'.
#' @param cl cluster you want to check.
#'
#' @return matched cellMarker data.frame
#' @import data.table
#' @export
check_marker <- function(marker, n, spc, cl = c()) {
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
