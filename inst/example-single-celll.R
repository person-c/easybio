library(Seurat)
savedir <- "extdata/example-single-cell"
# download file from here
# download.file("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz")

theme_set(theme_classic())
bf <- theme(
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  axis.text.x = element_text(angle = 60, size = 8),
  text = element_text(size = 8),
  legend.position = "none"
)

pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)
ggplot2::ggsave(file.path(savedir, "umap.png"), width = 5, height = 4)

marker <- FindAllMarkers(pbmc)
fwrite(marker, file.path(savedir, "exampleMarker.csv"))

## easybio single cell workflow to annotate cell
library(easybio)
check_marker(marker, 50, "Human", c(1, 5, 7)) |> DotPlot(pbmc, features = _) + scale_x_discrete(guide = guide_axis(angle = 60)) + bf
ggsave(file.path(savedir, "dotplot_1_5_7.png"), width = 10, height = 4)
check_marker(marker, 50, "Human", c(8)) |> DotPlot(pbmc, features = _) + scale_x_discrete(guide = guide_axis(angle = 60)) + bf
ggsave(file.path(savedir, "dotplot_8.png"), width = 4, height = 3)
check_marker(marker, 50, "Human", c(3)) |> DotPlot(pbmc, features = _) + scale_x_discrete(guide = guide_axis(angle = 60)) + bf
ggsave(file.path(savedir, "dotplot_3.png"), width = 4, height = 3)
check_marker(marker, 50, "Human", c(0, 2, 4, 6)) |> DotPlot(pbmc, features = _) + scale_x_discrete(guide = guide_axis(angle = 60)) + bf
ggsave(file.path(savedir, "dotplot_0_2_4_6.png"), width = 10, height = 4)


# According to the dotplot; We can annotate the clusters as follows
cluster2cell <- easybio::finsert(expression(
  c(1, 5) == "Monocyte",
  c(7) == "DC",
  c(8) == "megakaryocyte",
  c(3) == "B.cell",
  c(0, 2) == "Naive.CD8.T.cell",
  c(4) == "Cytotoxic.T.Cell",
  c(6) == "Natural.killer.cell"
))
cluster2cell <- setNames(cluster2cell, as.character(0:8))

# add annotation to metadata
pbmc@meta.data$anno <- cluster2cell[as.character(Idents(pbmc))]
DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "anno")
ggplot2::ggsave(file.path(savedir, "umap_anno.png"), width = 9, height = 6)


## You can get marker for a specific cell using easybio
easybio::get_marker(spc = "Human", cell = c("Monocyte", "Neutrophil"))
