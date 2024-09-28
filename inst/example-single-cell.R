# The data are downloaded from "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"

library(Seurat)
library(easybio)

x <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = x, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

markerTop50Matched <- matchCellMarker2(marker = pbmc.markers, n = 50, spc = "Human")
markerTop50Matched

# You can just use the top matched cell as the annotation
cl2cell <- markerTop50Matched[, head(.SD, 1), by = .(cluster)][, .(cluster, cell_name)]
cl2cell <- setNames(cl2cell[["cell_name"]], cl2cell[["cluster"]])
cl2cell

# or recheck the dot plot for similar clusters
cls <- list(
  c(1, 5, 7),
  c(8),
  c(3),
  c(0, 2, 4, 6)
)

# marker of possible cell for cluster 1, 4, 7
check_marker(pbmc.markers, 50, spc = "Human", cl = c(1, 5, 7))

# Check these markers' distribution of possible cell
srtDotPlot <- plotSeuratDot(srt = pbmc, cls = cls, marker = pbmc.markers, spc = "Human", n = 50)
srtDotPlot[[1]] # view other dotPlot `srtDotPlot[2]`, `srtDotPlot[3]`...

# According to the srtDotplot
cl2cell <- finsert(
  expression(
    c(3) == "B cell",
    c(8) == "Megakaryocyte",
    c(7) == "DC",
    c(1, 5) == "Monocyte",
    c(0, 2, 4) == "Naive CD8+ T cell",
    c(6) == "Natural killer cell"
  ),
  len = 9
)
cl2cell
pbmc@meta.data[["anno"]] <- cl2cell[as.character(Idents(pbmc))]
DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "anno")
