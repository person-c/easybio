## code to prepare `MSigDB` dataset goes here

download.file("http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_Seq.xlsx", "data-raw/Cell_marker_Seq.xlsx")

x <- readxl::read_xlsx("data-raw/Cell_marker_Seq.xlsx")
cellMarker2 <- data.table::setDT(x)
usethis::use_data(cellMarker2, internal = TRUE, overwrite = TRUE)
file.remove("data-raw/Cell_marker_Seq.xlsx")
