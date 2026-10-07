# ==============================================================================
# easybio Single-Cell RNA-Seq Annotation Workflow Example
# ==============================================================================
#
# This script demonstrates the complete workflow for annotating single-cell clusters
# using the `easybio` package. The process starts with a standard Seurat
# analysis pipeline to obtain cell clusters and their marker genes. Then, it
# uses `easybio`'s functions to perform automated annotation, interactive
# verification, and final manual curation.
#
# The data used here is the 3k PBMC dataset from 10x Genomics.
# Download URL: "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
#
# To run this example, ensure the data is downloaded and unzipped in your working directory.
# You can see the raw code of this example by running:
# file.show(system.file(package = "easybio", "example-single-cell.R"))
# ------------------------------------------------------------------------------

# Load necessary libraries
library(Seurat)
library(easybio)

# ---
# Step 1: Standard Seurat Pre-processing and Clustering
# ---
# This section covers the standard Seurat workflow. The goal is to normalize the
# data, find cell clusters, and identify marker genes for each cluster. These
# markers will be the input for the `easybio` annotation workflow.

# Load the 10x Genomics dataset
setwd("data-raw")
x <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

# Create the Seurat object with initial filtering
pbmc <- CreateSeuratObject(
  counts = x,
  project = "pbmc3k",
  min.cells = 3,
  min.features = 200
)

# Quality control: calculate mitochondrial DNA percentage and filter out low-quality cells
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
)

# Normalize and scale the data, and find variable features
pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 1e4
)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Perform dimensionality reduction (PCA and UMAP) and clustering
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Visualize the initial, unannotated clusters
tmp <- DimPlot(pbmc, reduction = "umap", label = TRUE)
print(tmp)

# Find marker genes for each cluster. This is the crucial input for match_ref.
# Keep Seurat's default min.pct here: it gates which genes are *tested* in either
# population, so a gene detected in only a handful of cells of the cluster can
# still come out as a positive marker. The detection rate of every marker is
# carried in the `pct.1` column instead, and match_ref() can gate on that from
# its own `min_pct` argument.
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)


# ---
# Step 2: Automated Annotation with `match_ref`
# ---
# This is the first core step of the `easybio` workflow. We use the marker genes
# found in the previous step to query the CellMarker 3.0 database and get a list
# of potential cell types for each cluster.

# The search is bounded to the tissues this sample could have come from. Every
# database entry carries a tissue_class (the organ or system) and a tissue_type
# (the sample the original study described), and match_ref() ANDs the two.
# Bounding the class is what matters here: the default searches every organ, and
# bounding the type as well removes evidence rather than correcting for
# anything. PBMC are peripheral blood mononuclear cells, so blood plus the bone
# marrow they are made in. available_tissue_type("Human", tissue_class = "Blood")
# lists what a class actually holds.
marker <- match_ref(
  marker = pbmc.markers,
  n = 50,
  spc = "Human",
  tissue_class = c("Blood", "Bone marrow")
)
# Let's look at the top results. Within each cluster the candidates are ranked by
# `uniqueN` (how many of our marker genes they matched), with `N` (the total
# number of matching database entries) as the tie-breaker. `pct_with` reports the
# detection rate of each matched gene, aligned with `ordered_symbol`, so you can
# see how much of the evidence is actually expressed.
marker |> head()

# For a quick first look, you can extract the top matched cell type for each cluster.
# This provides a preliminary, automated annotation.
cl2cell <- marker[, head(.SD, 1), by = .(cluster)][, .(cluster, cell_name)]
cl2cell <- setNames(cl2cell[["cell_name"]], cl2cell[["cluster"]])
print("Initial automated annotation based on top hits:")
cl2cell

pbmc@meta.data[["anno_initial"]] <- cl2cell[as.character(Idents(pbmc))]
tmp <- DimPlot(
  pbmc,
  reduction = "umap",
  label = TRUE,
  group.by = "anno_initial"
)
print(tmp)

# For a global view of every candidate, use `plot_possible_cell()`: one point per
# cluster and candidate, sized and coloured by the strength of the evidence.
tmp <- plot_possible_cell(
  marker[, head(.SD), by = .(cluster)],
  min_unique_n = 2
)
print(tmp)

# The same view can be coloured by how well the evidence is actually detected.
# With `value = "pct"` the fill shows the share of a candidate's matched markers
# whose detection rate reaches `min_pct`, which is the quickest way to spot an
# annotation that rests on barely expressed genes.
tmp <- plot_possible_cell(
  marker[, head(.SD), by = .(cluster)],
  min_unique_n = 2,
  value = "pct"
)
print(tmp)


# ---
# Step 3: Verification and Exploration with `check_marker` and `plot_seurat_dot`
# ---
# The automated annotation is just a starting point. This step is crucial for
# verifying the results and gaining confidence in the annotations.

# Question 1: "Why did the algorithm make these annotations?"
# Use `cis = TRUE` to see which of OUR marker genes matched the database,
# providing the evidence for the annotation.
tmp <- check_marker(marker, cl = c(1, 5, 7), top_cell_n = 2, cis = TRUE)

# Question 2: "Are these annotations correct?"
# Use `cis = FALSE` to retrieve the CANONICAL markers for the suggested cell types
# from the CellMarker 3.0 database. We can then check if these canonical markers
# are actually expressed in our clusters.
tmp <- check_marker(marker, cl = c(1, 5, 7), top_cell_n = 2, cis = FALSE)

# Now, let's visually confirm the expression of the supporting markers using a Dot Plot.
# We can pipe the results from `check_marker` directly into `plot_seurat_dot`.
# This plot shows the expression of the genes that led to the annotation (from cis = TRUE).
tmp <- check_marker(marker, cl = c(1, 5, 7), top_cell_n = 2, cis = TRUE) |>
  plot_seurat_dot(srt = pbmc)

# We can systematically check all interesting cluster groups.
# Here, we create a list of cluster groups based on the UMAP plot.
cls <- list(
  c(1, 5, 7),
  c(8),
  c(3),
  c(0, 2, 4, 6)
)

# Loop through the list and generate a dot plot for each group to inspect the evidence.
tmp <- lapply(cls, \(cl) {
  check_marker(marker, cl = cl, top_cell_n = 2, cis = TRUE) |>
    plot_seurat_dot(srt = pbmc) +
    ggplot2::ggtitle(
      paste0("Cluster ", paste(cl, collapse = ","), "'s possible cell types")
    )
})

print(tmp[[1]]) # Show the first plot as an example

# The entire workflow from annotation to visualization can be done in a single pipe:
tmp <- match_ref(marker = pbmc.markers, n = 50, spc = "Human") |>
  check_marker(cl = c(1, 5, 7), top_cell_n = 2, cis = TRUE) |>
  plot_seurat_dot(srt = pbmc)

print(tmp)
# ---
# Step 4: Final Manual Curation and Annotation
# ---
# After reviewing the dot plots and the evidence from `check_marker`, we can make
# an informed, final decision on the cell type for each cluster.

# Based on the dot plots, we create a final mapping.
# For example, we might conclude that cluster 3 is B cells, cluster 8 is Megakaryocytes, etc.
# `finsert` is a convenient helper function to create this mapping vector.
cl2cell <- finsert(
  list(
    c(3) ~ "B cell",
    c(8) ~ "Megakaryocyte",
    c(7) ~ "DC",
    c(1, 5) ~ "Monocyte",
    c(0, 2, 4) ~ "Naive CD8+ T cell",
    c(6) ~ "Natural killer cell"
  ),
  len = 9 # Ensure the vector has length for all clusters (0-8)
)
print("Final curated annotation:")
cl2cell

# Add the final, curated annotations to the Seurat object's metadata.
pbmc@meta.data[["anno"]] <- cl2cell[as.character(Idents(pbmc))]

# Visualize the final annotated UMAP.
tmp <- DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "anno")
print(tmp)

# ---
# Additional `easybio` Utility Functions
# ---
# The package also includes functions for direct queries.

# `get_marker`: Retrieve canonical markers for specific cell types directly.
get_marker(
  spc = "Human",
  cell = c("Monocyte", "Neutrophil"),
  number = 5,
  min_count = 1
)

# `plot_marker_distribution`: Visualize how a single marker is distributed across
# all cell types and tissues in the CellMarker 3.0 database.
plot_marker_distribution(mkr = "CD68")
