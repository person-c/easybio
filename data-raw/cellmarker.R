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

# Some rows store a bare number where a marker belongs (e.g. 2.716302.157),
# they cannot match any reference
x <- x[!grepl("^[0-9]+\\.[0-9]", marker)]

# Standardise gene name casing.
# Human symbols are all uppercase (HGNC)
x[species == "Human", let(marker = toupper(marker))]

# Mouse symbols are sentence case (Cd4, Gapdh). Only uniformly cased names are
# converted, so the ones the rule cannot reproduce stay intact: symbols carrying
# meaningful capitals (H2-K1, RP23-100J14.2), GenBank-style symbols that are
# officially uppercase (AA467197) and RIKEN cDNAs (9130008F23Rik).
sentence_case <- function(v) gsub("(^[[:alpha:]])", "\\U\\1", tolower(v), perl = TRUE)
x[
  species == "Mouse" &
    grepl("^[A-Za-z][A-Za-z0-9]*$", marker) &
    !grepl("^[A-Z]{2}[0-9]{5,6}$|RIK$", marker) &
    (grepl("^[A-Z0-9]+$", marker) | grepl("^[a-z0-9]+$", marker)),
  let(marker = sentence_case(marker))
]

# Mitochondrial genes take a lowercase prefix in mouse (mt-Nd1) where human uses
# uppercase (MT-ND1), rebuild the mouse spelling from whatever casing is stored
x[
  species == "Mouse" & grepl("^[Mm][Tt]-", marker),
  let(marker = paste0("mt-", sentence_case(sub("^[Mm][Tt]-", "", marker))))
]

# Keep only columns used by the package to reduce sysdata size.
# Rows are NOT deduplicated: repeated marker-cell pairs carry the
# literature-support counts that get_marker(min.count) relies on.
keep_cols <- c("species", "tissue_class", "tissue_type", "cell_name", "marker")
cellMarker3 <- x[, .SD, .SDcols = keep_cols]

# Build index on species for fast lookups
data.table::setindex(cellMarker3, species)

# xz rather than use_data()'s bzip2 default: the same object is ~23% smaller
# (1.24 MB against 1.62 MB), and it is decompressed once at install time, so
# the extra compression costs nothing at run time. `version` is pinned to 2,
# the format the file already uses, so a future change of use_data()'s default
# cannot quietly re-serialise it.
usethis::use_data(cellMarker3, internal = TRUE, overwrite = TRUE, compress = "xz", version = 2)
