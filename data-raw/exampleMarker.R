marker.file <- "./inst/extdata/example-single-cell/exampleMarker.csv"

exampleMarker <- data.table(marker.file)
devtools::use_data(exampleMarker)
