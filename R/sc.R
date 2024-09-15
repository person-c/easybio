#' Insert Specific Values into a Character Vector at Defined Positions
#'
#' This function constructs a character vector of a specified length, inserting
#' given values at positions determined by numeric indices. It is designed for
#' single cell annotation tasks, where specific annotations need to be placed
#' at certain positions in a vector.
#'
#' @param x An expression defining the value to insert and the positions at which
#'   to insert them. The expression should be a list of logical comparisons,
#'   where the left side is a numeric vector of positions and the right side is
#'   the corresponding character value to insert.
#' @param len The desired length of the output character vector. If the specified
#'   positions exceed this length, the vector will be padded with the `na` value.
#' @param setname A logical value indicating whether to set names for the elements
#'   of the vector. If `TRUE`, names are set as character representations of the
#'   positions from 0 to the length of the vector minus one.
#' @param na The default value to use for positions not specified in `x`. This
#'   value is also used to pad the vector if its length exceeds the positions
#'   specified in `x`.
#'
#' @return A named character vector with the specified values inserted at given
#'   positions and padded with the `na` value if necessary.
#' @export
#'
#' @examples
#' # Example usage:
#' # Insert "Neutrophil" at positions 0, 1, 3 and "Macrophage" at positions 2, 4, 8
#' # in a vector of length 10, with "Unknown" as the default value.
#' library(easybio)
#' annotated_vector <- finsert(
#'   x = expression(
#'     c(0, 1, 3) == "Neutrophil",
#'     c(2, 4, 8) == "Macrophage"
#'   ),
#'   len = 10,
#'   na = "Unknown"
#' )
#' print(annotated_vector)
finsert <- function(
    x = expression(
      c(0, 1, 3) == "Neutrophil",
      c(2, 4, 8) == "Macrophage"
    ),
    len = integer(),
    setname = TRUE,
    na = "Unknown") {
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

  clIdx <- which(sapply(x, is.numeric))
  maxL <- max(unlist(x[clIdx]))

  if (!missing(len) && len > maxL + 1) v <- append(v, rep(NA_character_, len - maxL - 1))
  v[is.na(v)] <- na

  if (setname) names(v) <- as.character(0:(length(v) - 1))
  v
}


#' Retrieve Markers for Specific Cells from cellMarker2
#'
#' This function extracts a list of markers for one or more cell types from the
#' `cellMarker2` dataset. It allows filtering by species, cell type, the number
#' of markers to retrieve, and a minimum count threshold for marker occurrences.
#'
#' @param spc A character string specifying the species, which can be either
#'   'Human' or 'Mouse'.
#' @param cell A character vector of cell types for which to retrieve markers.
#' @param number An integer specifying the number of top markers to return for
#'   each cell type.
#' @param min.count An integer representing the minimum number of times a marker
#'   must have been reported to be included in the results.
#'
#' @return A named list where each name corresponds to a cell type and each
#'   element is a vector of marker names.
#' @import data.table
#' @export
#'
#' @examples
#' # Example usage:
#' # Retrieve the top 5 markers for 'Macrophage' and 'Monocyte' cell types in humans,
#' # with a minimum count of 1.
#' library(easybio)
#' markers <- get_marker(spc = "Human", cell = c("Macrophage", "Monocyte"))
#' print(markers)
get_marker <- function(spc, cell = character(), number = 5, min.count = 1) {
  . <- NULL
  species <- cell_name <- N <- marker <- NULL

  is_exists <- cell %chin% unique(cellMarker2[["cell_name"]])


  sapply(cell[!is_exists], \(x) {
    idx <- grep(
      gsub("\\s+", "", x),
      gsub("\\s+", "", unique(cellMarker2[["cell_name"]])),
      ignore.case = TRUE
    )
    psblCell <- paste0(unique(cellMarker2[["cell_name"]])[idx], collapse = "\n")
    message(sprintf("%s doesn't exist in the CellMarker2 database; Maybe you means?\n%s\n", x, psblCell))
  })

  if (all(!is_exists)) {
    return(NULL)
  }

  marker <- cellMarker2[.(spc, cell), .SD, on = .(species, cell_name)]
  marker <- marker[, .N, by = .(cell_name, marker)]
  marker <- marker[N > min.count, na.omit(.SD)[order(-N)] |> head(number), by = .(cell_name)]
  marker <- marker[, .(marker = .(marker)), by = .(cell_name)]
  marker <- setNames(marker[["marker"]], marker[["cell_name"]])

  marker
}

#' Match Markers with cellMarker2 Dataset
#'
#' This function matches markers from the `FindAllMarkers` output with the
#' `cellMarker2` dataset, filtering by species and selecting the top genes based
#' on their average log2 fold change and adjusted p-values.
#'
#' @param marker A data frame of markers obtained from the `FindAllMarkers`
#'   function, expected to contain columns such as `avg_log2FC`, `p_val_adj`,
#'   and `gene`.
#' @param n An integer specifying the top number of genes to match from the
#'   input markers.
#' @param spc A character string specifying the species, which can be either
#'   'Human' or 'Mouse'.
#'
#' @return A data frame containing matched markers from the `cellMarker2`
#'   dataset, with additional columns indicating the number of matches and
#'   ordered symbols.
#' @import data.table
#' @export
#'
#' @examples
#' # Example usage:
#' # Match the top 50 markers from the pbmc.markers dataset with the Human
#' # species in the cellMarker2 dataset.
#' library(easybio)
#' data(pbmc.markers)
#' matched_markers <- matchCellMarker2(pbmc.markers, n = 50, spc = "Human")
#' print(matched_markers)
matchCellMarker2 <- function(marker, n, spc) {
  . <- markerWith <- NULL
  species <- avg_log2FC <- p_val_adj <- cluster <- gene <- cell_name <- N <- NULL

  marker <- copy(marker)
  setDT(marker)

  cellMarker2 <- cellMarker2[.(spc), on = .(species)]
  marker <- marker[avg_log2FC > 0 & p_val_adj < 0.05, .SD[order(-avg_log2FC)][1:n], keyby = .(cluster)]

  marker <- marker[cellMarker2, on = "gene==marker", nomatch = NULL]
  marker <- marker[, .(markerWith = .(gene), N = .N), by = .(cluster, cell_name)]
  marker <- marker[N > 0, .SD[order(-N)], keyby = .(cluster)]


  marker[, let(uniqueN = sapply(markerWith, FUN = \(x) length(unique(x))))]
  marker[, let(ordered_symbol = lapply(markerWith, FUN = \(x) names(sort(unclass(table(x)), TRUE))))]
  marker[, let(orderN = lapply(markerWith, \(x) as.integer(sort(unclass(table(x)), TRUE))))]
  setcolorder(marker, c("cluster", "cell_name", "uniqueN", "N", "ordered_symbol", "orderN", "markerWith"))
  marker
}

#' Verify Markers for Specific Clusters Using matchCellMarker
#'
#' This function checks the markers for specified clusters returned by the
#' `matchCellMarker2` function. It allows users to filter by species, cluster,
#' and to specify whether to consider cis or trans interactions.
#'
#' @param marker A data frame of markers obtained from `Seurat::FindAllMarkers`.
#' @param n An integer specifying the top number of genes to match from the
#'   input markers.
#' @param spc A character string specifying the species, which can be either
#'   'Human' or 'Mouse'.
#' @param cl An integer or vector of integers specifying the clusters to check.
#' @param topcellN An integer specifying the number of top cells to check for
#'   each cluster.
#' @param cis A logical value indicating whether to check marker directly from the top symbol of
#'   \code{matchCellMarker2} or re-search marker for top cell in cellMarker2.
#'
#' @return A named list where each name corresponds to a cell type and each
#'   element is a vector of marker names.
#' @import data.table
#' @export
#'
#' @examples
#' # Example usage:
#' # Check the top 50 markers for clusters 1, 4, and 7 in the Human species.
#' library(easybio)
#' data(pbmc.markers)
#' verified_markers <- check_marker(pbmc.markers, n = 50, spc = "Human", cl = c(1, 4, 7))
#' print(verified_markers)
check_marker <- function(marker, n, spc, cl = c(), topcellN = 2, cis = FALSE) {
  . <- NULL
  cell_name <- cluster <- NULL

  marker_matched <- matchCellMarker2(marker, n, spc)
  marker_matched <- marker_matched[.(factor(cl)), .SD, on = .(cluster)]

  if (cis) {
    topmarker <- marker_matched[, head(.SD, topcellN), by = .(cluster)]
    topmarker <- setNames(topmarker[["ordered_symbol"]], topmarker[["cell_name"]])
  } else {
    topcell <- marker_matched[, head(.SD, topcellN), keyby = .(cluster)][, unique(cell_name)]
    topmarker <- get_marker(spc, cell = topcell, number = 10, min.count = 0)
  }

  topmarker
}


#' Create Dot Plots for Markers from check_marker
#'
#' This function generates dot plots for the markers obtained from the
#' `check_marker` function for specified cluster groups within a Seurat object.
#' The plots are saved to a temporary directory.
#'
#' @param srt A Seurat object containing the single-cell data.
#' @param cls A list containing cluster groups to check. Each element of the list
#'   should correspond to a cluster or a group of clusters for which to generate
#'   dot plots.
#' @param ... Additional parameters to pass to the `check_marker` function.
#'
#' @return The function returns the temporary directory invisibly.
#' @import ggplot2
#' @export
plotSeuratDot <- function(srt, cls, ...) {
  dotplotList <- lapply(cls, \(cl) {
    features <- check_marker(..., cl = cl)

    if (anyDuplicated(unlist(features)) > 0) {
      features <- unique(list2dt(features), by = "value")
      warning("Duplicated markers are removed!")

      features <- split(features[["value"]], features[["name"]])
    }

    Seurat::DotPlot(srt, features = features) +
      scale_x_discrete(
        guide = guide_axis(
          angle = 60,
          theme = theme(text = element_text(size = 4))
        )
      ) +
      theme(
        axis.text = element_text(size = 4),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, angle = 30, vjust = 0.1, hjust = 0)
      )
  })

  dotplotName <- vapply(cls, \(cl) {
    paste0("clusters_", paste0(cl, collapse = "_"))
  }, "character")

  names(dotplotList) <- dotplotName
  dotplotList
}

#' Plot Distribution of a Marker Across Tissues and Cell Types
#'
#' This function creates a dot plot displaying the distribution of a specified marker across
#' different tissues and cell types, based on data from the CellMarker2.0 database.
#'
#' @param mkr character, the name of the marker to be plotted.
#'
#' @return A ggplot2 object representing the distribution of the marker.
#' @import ggplot2
#' @import data.table
#'
#' @export
#' @examples
#' plotMarkerDistribution("CD14")
plotMarkerDistribution <- function(mkr = character()) {
  . <- cell_name <- tissue_class <- cell_name <- N <- marker <- NULL
  tmp <- cellMarker2[.(mkr), .SD, on = .(marker), by = .(cell_name, tissue_class)]
  tmp <- tmp[, .N, by = .(cell_name, tissue_class)]

  p <- ggplot(tmp, aes(x = cell_name, y = tissue_class)) +
    geom_point(aes(size = N, color = N)) +
    scale_x_discrete(guide = guide_axis(angle = 60)) +
    scale_color_distiller(direction = 1) +
    theme_publication()

  p
}

#' Plot Possible Cell Distribution Based on matchCellMarker2() Results
#'
#' This function creates a dot plot to visualize the distribution of possible cell types
#' based on the results from the `matchCellMarker2()` function, utilizing data from the CellMarker2.0 database.
#'
#' @param marker data.table, the result from the `matchCellMarker2()` function.
#' @param min.uniqueN integer, the minimum number of unique marker genes that must be matched for a cell type to be included in the plot. Default is 2.
#'
#' @return A ggplot2 object representing the distribution of possible cell types.
#' @import ggplot2
#' @import data.table
#' @export
plotPossibleCell <- function(marker, min.uniqueN = 2) {
  cluster <- cell_name <- N <- NULL
  p <- ggplot(marker[uniqueN > min.uniqueN], aes(x = cell_name, y = cluster)) +
    geom_point(aes(size = N, color = N)) +
    scale_x_discrete(guide = guide_axis(angle = 60)) +
    scale_color_distiller(direction = 1) +
    theme_publication()

  p
}



.tuneParameters <- function(srt, resolution, N, spc) {
  cluster <- NULL
  srt <- suppressMessages(Seurat::FindClusters(srt, resolution = resolution))
  srt.markers <- Seurat::FindAllMarkers(srt, only.pos = TRUE)

  markerMatched <- matchCellMarker2(marker = srt.markers, n = N, spc = spc)
  cl2cell <- markerMatched[, head(.SD, 1), by = cluster][, 1:4]
  cl2cell <- setNames(cl2cell[["cell_name"]], as.character(cl2cell[["cluster"]]))
  srt@meta.data[["CellMarker2.0"]] <- cl2cell[as.character(Seurat::Idents(srt))]

  p <- Seurat::DimPlot(srt,
    reduction = "umap",
    label = TRUE, label.size = 1,
    pt.size = 0.6, repel = TRUE,
    group.by = "CellMarker2.0"
  ) +
    labs(title = sprintf("resolution: %s N: %s", resolution, N)) +
    guides(color = guide_legend(override.aes = list(size = 0.5))) +
    theme_publication(base_size = 8)

  p
}
#' Optimize Resolution and Gene Number Parameters for Cell Type Annotation
#'
#' This function tunes the `resolution` parameter in `Seurat::FindClusters()` and the number of top differential genes (`N`) to obtain different cell type annotation results. The function generates UMAP plots for each parameter combination, allowing for a comparison of how different settings affect the clustering and annotation.
#'
#' @param srt Seurat object, the input data object to be analyzed.
#' @param resolution numeric vector, a vector of resolution values to be tested in `Seurat::FindClusters()`.
#' @param N integer vector, a vector of values indicating the number of top differential genes to be used for matching in `matchCellMarker2()`.
#' @param spc character, the species parameter for the `matchCellMarker2()` function, specifying the organism.
#'
#' @return A list of ggplot2 objects, each representing a UMAP plot generated with a different combination of resolution and N parameters.
#' @import ggplot2
#' @import data.table
#' @export
tuneParameters <- function(srt, resolution = numeric(), N = integer(), spc) {
  parameters <- CJ(resolution = resolution, N = N)

  parameterPlot <- Map(
    f = function(x, y) .tuneParameters(srt, x, y, spc),
    x = parameters[["resolution"]],
    y = parameters[["N"]]
  )

  parameterPlot
}

# Used for future
# Assign weight for different markers
.get_marker_weight <- function(spc, cell = character(), min.count = 0, power = 2) {
  . <- NULL
  species <- cell_name <- N <- marker <- NULL

  marker <- cellMarker2[.(spc, cell), .SD, on = .(species, cell_name)]
  marker <- marker[, .(N = .N), by = .(cell_name, marker)]
  marker[, let(weight = N^power / sum(N^power)), by = cell_name]
  marker
}

# Find markers for similar clusters
.findMarkers <- function(SeuratObject, cls = list()) {
  lapply(cls, function(x) {
    Seurat::FindMarkers(SeuratObject, ident.1 = x, group.by = "seurat_clusters")
  })
}
