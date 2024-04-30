#' Quality control and normalize data
#'
#' @param sc Seurat object
#' @param .subset Quality control expression
#' @param .pattern mitochondrial gene pattern.
#'
#' @return Seurat object
#' @import Seurat
#' @export
dprocess.Seurat <- function(
    sc,
    .subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5,
    .pattern = "^mt-") {
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = .pattern)
  sc <- eval(substitute(subset(sc, subset = .subset)))

  sc <- NormalizeData(sc)
  sc <- FindVariableFeatures(sc)
  sc <- ScaleData(sc, features = rownames(sc))

  options(future.globals.maxSize = 3e+09)
  sc <- SCTransform(sc, vars.to.regress = "percent.mt", verbose = FALSE)
}


#' Find clusters in reduced dimension space.
#'
#' @param sc Seurat object
#' @param resolution higher resolution, more clusters.
#' @param integrate integrate different layers in low dimension space.
#'
#' @return Seurat object
#' @import Seurat
#' @export
dimReduce <- function(sc, resolution, integrate = TRUE) {
  sc <- RunPCA(sc, npcs = 30, verbose = FALSE)

  if (integrate) {
    sc <- IntegrateLayers(
      object = sc,
      method = RPCAIntegration,
      normalization.method = "SCT",
      verbose = FALSE
    )

    sc <- FindNeighbors(sc, dims = 1:30, reduction = "integrated.dr")
    sc <- FindClusters(sc, resolution = resolution)

    sc <- RunUMAP(sc, reduction = "integrated.dr", dims = 1:30, reduction.name = "umap.dr")
    p <- DimPlot(sc, reduction = "umap.dr", label = TRUE, split.by = "orig.ident") +
      scale_colour_brewer(palette = "Paired")
    print(p)

    return(sc)
  }

  sc <- FindNeighbors(sc, dims = 1:30)
  sc <- FindClusters(sc, resolution = resolution)

  sc <- RunUMAP(sc, dims = 1:30)
  p <- DimPlot(sc, reduction = "umap", label = TRUE, split.by = "orig.ident") +
    scale_colour_brewer(palette = "Paired")

  p

  sc
}


#' match markers from cellMarker2
#'
#' @param sc Seurat object
#' @param n top number of genes to match
#'
#' @return matched cellMarker data.frame
#' @import Seurat
#' @export
markerView <- function(sc, n) {
  sc <- PrepSCTFindMarkers(sc)
  marker <- FindAllMarkers(sc, assay = "SCT", only.pos = TRUE)
  setDT(marker)
  marker <- marker[celllMarker2, on = "gene==marker", nomatch = NULL][, .SD[order(-avg_log2FC)][1:n], keyby = .(cluster)]
  marker <- marker[, .SD[duplicated(gene)][, .(markerWith = .(gene), N = .N)], by = .(cluster, cell_name)]
  marker[, .SD[order(-N)], by = .(cluster)][N != 0]
}
