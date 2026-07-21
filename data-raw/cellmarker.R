## code to prepare `cellMarker2` dataset from CellMarker 3.0
## Source: data-raw/single_cell_marker/single_cell_marker.txt

x <- data.table::fread("data-raw/single_cell_marker/single_cell_marker.txt")

# Normalize marker: prefer standardised gene symbol, fall back to raw marker
x[, let(marker = fcoalesce(symbol, marker))]

# Standardise gene name casing
x[species == "Human", let(marker = toupper(marker))]
x[species == "Mouse", let(marker = gsub("(^[[:alpha:]])", "\\U\\1",
  tolower(marker),
  perl = TRUE
))]

# Keep only columns used by the package to reduce sysdata size
keep_cols <- c(
  "species", "tissue_class", "tissue_type",
  "cell_name", "marker",
  "uberon_id", "disease", "cellontology_id",
  "gene_id", "gene_type", "uniprot_id",
  "pmid", "journal", "year",
  "marker_source"
)
cellMarker2 <- x[, .SD, .SDcols = intersect(keep_cols, colnames(x))]

# Build index on species for fast lookups
data.table::setindex(cellMarker2, species)

usethis::use_data(cellMarker2, internal = TRUE, overwrite = TRUE)
