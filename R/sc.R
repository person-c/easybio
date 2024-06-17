#' match markers from cellMarker2
#'
#' @param marker markers from `FindAllMarkers`
#' @param n top number of genes to match
#'
#' @return matched cellMarker data.frame
#' @export
matchCellMarker2 <- function(marker, n, species) {
  setDT(marker)
  cellMarker2 <- cellMarker2[.(species), .SD, on = .(species)]
  marker <- marker[, .SD[order(-avg_log2FC)][1:n], keyby = .(cluster)][cellMarker2, on = "gene==marker", nomatch = NULL]
  marker <- marker[, .(markerWith = .(gene), N = .N), by = .(cluster, cell_name)]
  marker[N > 0, .SD[order(-N)], keyby = .(cluster)]
}
