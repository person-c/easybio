#' Create a Vector from an Index-to-Label Map
#'
#' Constructs a character vector by mapping labels to specified 0-based numeric
#' indices. This is a utility function often used in single-cell analysis to
#' assign cell type annotations to cluster IDs.
#'
#' @param x The mapping of indices to labels. This can be provided in two formats:
#'   \itemize{
#'     \item A \code{list} of formulas, e.g., \code{list(c(0, 1) ~ "LabelA", 2 ~ "LabelB")}.
#'     \item An \code{expression} object, e.g., \code{expression(c(0, 1) == "LabelA", 2 == "LabelB")}.
#'   }
#' @param len An optional integer specifying the minimum length of the output
#'   vector. If the highest index in \code{x} is greater than \code{len}, the
#'   vector will be automatically extended.
#' @param setname A logical value. If \code{TRUE} (the default), the elements of
#'   the output vector are named with their corresponding 0-based index (e.g., "0", "1", "2", ...).
#' @param na The character value used to fill positions that are not specified in
#'   the mapping. Defaults to "Unknown".
#'
#' @return A character vector with the specified labels at the given positions.
#'   The vector is named with 0-based indices if \code{setname} is \code{TRUE}.
#'
#' @export
#'
#' @examples
#' # --- Example 1: Using the default formula list format ---
#' # This is the recommended and default usage.
#' mapping_formula <- list(
#'   c(0, 1, 3) ~ "Neutrophil",
#'   c(2, 4, 8) ~ "Macrophage"
#' )
#' finsert(mapping_formula)
#'
#' # --- Example 2: Using the expression format for backward compatibility ---
#' mapping_expr <- expression(
#'   c(0, 1, 3) == "Neutrophil",
#'   c(2, 4, 8) == "Macrophage"
#' )
#' finsert(mapping_expr, len = 10, na = "Unassigned")
#'
finsert <- function(
    x = list(
      c(0, 1, 3) ~ "Neutrophil",
      c(2, 4, 8) ~ "Macrophage"
    ),
    len = integer(),
    setname = TRUE,
    na = "Unknown") {
  x <- if (is.expression(x)) lapply(x, .exprs2formula) else x

  maxL <- max(unlist(sapply(x, \(.x) eval(.x[[2]]), simplify = FALSE)))
  v <- rep(na, if (!missing(len) && len > (maxL + 1)) len else maxL + 1)
  invisible(lapply(x, \(.x) v[eval(.x[[2]]) + 1] <<- .x[[3]]))

  if (setname) names(v) <- as.character(0:(length(v) - 1))

  return(v)
}

#' Retrieve Available Tissue Classes for a Given Species
#'
#' This function extracts and returns a unique list of available tissue classes
#' from the CellMarker2.0 database for a specified species.
#'
#' @param spc A character string specifying the species (e.g., "Human" or "Mouse").
#'
#' @return A character vector of unique tissue classes available for the given species.
#' If no tissue classes are found, an empty vector is returned.
#'
#' @seealso \code{\link{available_tissue_type}}, \code{\link{get_marker}}
#'
#' @examples
#' # Get all tissue classes for Human
#' available_tissue_class("Human")
#'
#' @export
available_tissue_class <- function(spc) {
  assert_subset(spc, c("Human", "Mouse"), empty.ok = FALSE)

  species <- NULL
  available_ele(cellMarker2, "tissue_class", subset = species == spc)
}

#' Retrieve Available Tissue Types for a Given Species
#'
#' This function extracts and returns a unique list of available tissue types
#' from the CellMarker2.0 database for a specified species.
#'
#' @param spc A character string specifying the species (e.g., "Human" or "Mouse").
#'
#' @return A character vector of unique tissue types available for the given species.
#' If no tissue types are found, an empty vector is returned.
#'
#' @seealso \code{\link{available_tissue_class}}, \code{\link{get_marker}}
#'
#' @examples
#' # Get all tissue types for Human
#' available_tissue_type("Human")
#'
#' @export
available_tissue_type <- function(spc) {
  assert_subset(spc, c("Human", "Mouse"), empty.ok = FALSE)

  species <- NULL
  available_ele(cellMarker2, "tissue_type", subset = species == spc)
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
#' @param tissueClass A character specifying the tissue classes, default `available_tissue_class(spc)`.
#' @param tissueType A character specifying the tissue types, default `available_tissue_type(spc)`.
#'
#' @return A named list where each name corresponds to a cell type and each
#'   element is a vector of marker names.
#' @export
#'
#' @examples
#' # Example usage:
#' # Retrieve the top 5 markers for 'Macrophage' and 'Monocyte' cell types in humans,
#' # with a minimum count of 1.
#' library(easybio)
#' markers <- get_marker(spc = "Human", cell = c("Macrophage", "Monocyte"))
#' print(markers)
#' # Example with a typo in cell name
#' markers_typo <- get_marker(spc = "Human", cell = c("Macrophae", "Monocyte"))
get_marker <- function(
    spc, cell = character(),
    tissueClass = available_tissue_class(spc),
    tissueType = available_tissue_type(spc),
    number = 5, min.count = 1) {
  . <- tissue_type <- tissue_class <- NULL
  species <- cell_name <- N <- marker <- NULL

  all_cell_names <- available_ele(cellMarker2, "cell_name", subset = species == spc)
  is_exists <- cell %chin% all_cell_names

  not_found_cells <- cell[!is_exists]
  if (length(not_found_cells) > 0) {
    suggestions <- vapply(not_found_cells, FUN.VALUE = "character", FUN = function(x) {
      # 1. Fuzzy match with adist for typos
      distances <- adist(x, all_cell_names, ignore.case = TRUE, partial = FALSE)
      min_dist <- min(distances)

      # Heuristic for a "good" match (e.g., distance <= 2)
      if (min_dist <= 2) {
        possible_matches <- all_cell_names[which(distances == min_dist)]
        return(paste(possible_matches, collapse = " or "))
      }

      # 2. Fallback to grep for partial/substring matches
      grep_idx <- grep(x, all_cell_names, ignore.case = TRUE)
      if (length(grep_idx) > 0) {
        return(paste(all_cell_names[grep_idx], collapse = " or "))
      }

      return("") # No suggestion found
    })

    # Format and print message for cells with suggestions
    has_suggestion <- nchar(suggestions) > 0
    if (any(has_suggestion)) {
      msg_lines <- sprintf(
        "- For '%s', did you mean: %s?",
        names(suggestions[has_suggestion]),
        suggestions[has_suggestion]
      )
      message("Some cell types not found. Suggestions:\n", paste(msg_lines, collapse = "\n"))
    }

    # Report cells for which no suggestion could be found
    if (any(!has_suggestion)) {
      message(
        "Could not find any matches for: ",
        paste(names(suggestions[!has_suggestion]), collapse = ", ")
      )
    }
  }

  if (all(!is_exists)) {
    message("No valid cell types provided to fetch markers. Returning NULL.")
    return(NULL)
  }

  # Proceed with only the cell names that exist
  valid_cells <- cell[is_exists]

  cellMarker2_filtered <- cellMarker2[tissue_class %chin% tissueClass & tissue_type %chin% tissueType]
  marker <- cellMarker2_filtered[.(spc, valid_cells), .SD, on = .(species, cell_name), nomatch = NULL]

  if (is.null(marker) || nrow(marker) == 0) {
    return(NULL)
  }

  marker <- marker[, .N, by = .(cell_name, marker)]
  marker <- marker[N >= min.count, na.omit(.SD)[order(-N)] |> head(number), by = .(cell_name)]
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
#' @param tissueClass A character specifying the tissue classes, default `available_tissue_class(spc)`.
#' @param tissueType A character specifying the tissue types, default `available_tissue_type(spc)`.
#'
#' @return A data frame containing matched markers from the `cellMarker2`
#'   dataset, with additional columns indicating the number of matches and
#'   ordered symbols.
#' @export
#'
#' @examples
#' # Example usage:
#' # Match the top 50 differential genes from the pbmc.markers dataset with the Human
#' # species in the cellMarker2 dataset.
#' \dontrun{
#' library(easybio)
#' data(pbmc.markers)
#' matchCellMarker2(pbmc.markers, n = 50, spc = "Human")[]
#' }
matchCellMarker2 <- function(
    marker, n, spc,
    tissueClass = available_tissue_class(spc),
    tissueType = available_tissue_type(spc)) {
  . <- markerWith <- tissue_class <- tissue_type <- NULL
  species <- avg_log2FC <- p_val_adj <- cluster <- gene <- cell_name <- N <- NULL

  marker <- copy(marker)
  setDT(marker)


  cellMarker2 <- cellMarker2[.(spc), .SD, on = .(species), nomatch = NULL]
  cellMarker2 <- cellMarker2[tissue_class %chin% tissueClass & tissue_type %chin% tissueType]

  marker <- marker[avg_log2FC > 0 & p_val_adj < 0.05, .SD[order(-avg_log2FC)][1:n], keyby = .(cluster)]

  marker <- marker[cellMarker2, on = "gene==marker", nomatch = NULL]
  marker <- marker[, .(markerWith = .(gene), N = .N), by = .(cluster, cell_name)]
  marker <- marker[N > 0, .SD[order(-N)], keyby = .(cluster)]


  marker[, let(uniqueN = sapply(markerWith, FUN = \(x) uniqueN(x)))]
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
#' @param tissueClass A character specifying the tissue classes, default `available_tissue_class(spc)`.
#' @param tissueType A character specifying the tissue types, default `available_tissue_type(spc)`.
#' @param cl An integer or vector of integers specifying the clusters to check.
#' @param topcellN An integer specifying the number of top cells to check for
#'   each cluster.
#' @param cis A logical value indicating whether to check marker directly from the top symbol of
#'   \code{matchCellMarker2} or re-search marker for top cell in cellMarker2.
#'
#' @return A named list where each name corresponds to a cell type and each
#'   element is a vector of marker names.
#' @export
#'
#' @examples
#' # Example usage:
#' # Check the top 50 markers for clusters 1, 4, and 7 in the Human species.
#' \dontrun{
#' library(easybio)
#' data(pbmc.markers)
#' check_marker(pbmc.markers, n = 50, spc = "Human", cl = c(1, 4, 7))
#' }
check_marker <- function(
    marker, n, spc,
    tissueClass = available_tissue_class(spc),
    tissueType = available_tissue_type(spc),
    cl = c(), topcellN = 2, cis = FALSE) {
  . <- cell_name <- cluster <- NULL

  marker_matched <- matchCellMarker2(marker, n, spc, tissueClass, tissueType)
  marker_matched <- marker_matched[.(factor(cl)), .SD, on = .(cluster)]

  if (cis) {
    topmarker <- marker_matched[, head(.SD, topcellN), by = .(cluster)]
    topmarker <- setNames(topmarker[["ordered_symbol"]], topmarker[["cell_name"]])
  } else {
    topcell <- marker_matched[, head(.SD, topcellN), keyby = .(cluster)][, unique(cell_name)]
    topmarker <- get_marker(
      spc,
      cell = topcell,
      tissueClass = tissueClass,
      tissueType = tissueType,
      number = 10,
      min.count = 0
    )
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
#' @return A list containing multiple DotPlot.
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
#'
#' @export
#' @examples
#' \dontrun{
#' plotMarkerDistribution("CD14")
#' }
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

.exprs2formula <- function(expr) {
  formula(paste(deparse(expr[[2]]), "~", deparse(expr[[3]])))
}
