setwd("data-raw")

fn <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
download.file(fn, "pbmc.tar.gz")
untar("pbmc.tar.gz", list = TRUE)
untar("pbmc.tar.gz")


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
ggplot2::ggsave("UMAP_Raw.png", width = 4.62, height = 3.26)


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
data.table::fwrite(pbmc.markers, "pbmc.markers.csv")

markerTop50Matched <- matchCellMarker2(marker = pbmc.markers, n = 50, spc = "Human")
cls <- list(
  c(1, 5, 7),
  c(8),
  c(3),
  c(0, 2, 4, 6)
)

dotplotList <- plotSeuratDot(srt = pbmc, cls = cls, marker = pbmc.markers, spc = "Human", n = 50)

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

pbmc@meta.data[["anno"]] <- cl2cell[as.character(Idents(pbmc))]
DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "anno")
ggplot2::ggsave("UMAP_Anno.png", width = 6.35, height = 4.51)

usethis::use_data(pbmc.markers, compress = "xz", overwrite = TRUE)
fdel <- list.files(pattern = "tar", full.names = TRUE)
file.remove(fdel)
unlink("filtered_gene_bc_matrices/", recursive = TRUE)
