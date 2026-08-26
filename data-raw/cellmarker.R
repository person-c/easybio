## code to prepare `cellMarker3` dataset from CellMarker 3.0
## Source: data-raw/single_cell_marker/single_cell_marker.txt
## (raw 142 MB file, git-ignored; see .gitignore)

library(data.table)

x <- data.table::fread("data-raw/single_cell_marker/single_cell_marker.txt")

# Normalize marker: prefer standardised gene symbol, fall back to raw marker
# (empty strings in symbol are treated as missing)
x[, let(marker = fifelse(is.na(symbol) | symbol == "", marker, symbol))]

# Drop rows without an identifiable marker
x <- x[!is.na(marker) & marker != ""]

# Standardise gene name casing
x[species == "Human", let(marker = toupper(marker))]
x[species == "Mouse", let(marker = gsub("(^[[:alpha:]])", "\\U\\1",
  tolower(marker),
  perl = TRUE
))]

# Keep only columns used by the package to reduce sysdata size.
# Rows are NOT deduplicated: repeated marker-cell pairs carry the
# literature-support counts that get_marker(min.count) relies on.
keep_cols <- c("species", "tissue_class", "tissue_type", "cell_name", "marker")
cellMarker3 <- x[, .SD, .SDcols = keep_cols]

# Build index on species for fast lookups
data.table::setindex(cellMarker3, species)

usethis::use_data(cellMarker3, internal = TRUE, overwrite = TRUE)
