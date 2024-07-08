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
  setDT(marker)

  cellMarker2 <- subset(cellMarker2, species == spc)
  marker <- marker[avg_log2FC > 0 & p_val_adj < 0.05, .SD[order(-avg_log2FC)][1:n], keyby = .(cluster)][cellMarker2, on = "gene==Symbol", nomatch = NULL]
  marker <- marker[, .(markerWith = .(gene), N = .N), by = .(cluster, cell_name)]
  marker <- marker[N > 0, .SD[order(-N)], keyby = .(cluster)]

  topcell <- marker[, head(.SD, 2), keyby = .(cluster)][, unique(cell_name)]
  topmarker <- cellMarker2[.(topcell), let(count = .N), on = .(cell_name), by = .(cell_name, Symbol)][, .SD[order(-count)] |> head(10), by = .(cell_name)]

  setattr(marker, "topmarker", topmarker)
  marker
}
