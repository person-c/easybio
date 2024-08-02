## code to prepare `MSigDB` dataset goes here

download.file("http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_Seq.xlsx", "data-raw/Cell_marker_Seq.xlsx")

x <- readxl::read_xlsx("data-raw/Cell_marker_Seq.xlsx")
cellMarker2 <- data.table::setDT(x)
cellMarker2[species == "Human", let(marker = toupper(marker))]
cellMarker2[species == "Mouse", let(marker = gsub("(^[[:alpha:]])", "\\U\\1", tolower(marker), perl = TRUE))]
cellMarker2[species == "Human", let(Symbol = toupper(Symbol))]
cellMarker2[species == "Mouse", let(Symbol = gsub("(^[[:alpha:]])", "\\U\\1", tolower(Symbol), perl = TRUE))]
cellMarker2[, let(Symbol = fcoalesce(Symbol, marker))]
setnames(cellMarker2, c("marker", "Symbol"), c("raw_marker", "marker"))
usethis::use_data(cellMarker2, internal = TRUE, overwrite = TRUE)
file.remove("data-raw/Cell_marker_Seq.xlsx")
