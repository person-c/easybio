# Gene expression aligned against hg38
library(TCGAbiolinks)
library(easybio)
library(data.table)

# data
query <- GDCquery(
  project = "TCGA-CHOL",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification"
)
GDCdownload(query = query)
data <- GDCprepare(query = query)

# data
lt <- prepare_tcga(data)
lt$all$sampleInfo[["group"]] <- fifelse(lt$all$sampleInfo$sample_type %ilike% "Tumor", "Tumor", "Normal")

# limma-voom workflow
x <- dgeList(lt$all$exprCount, lt$all$sampleInfo, lt$all$featuresInfo)
x <- dprocess_dgeList(x, "group", 10)
efit <- limmaFit(x, "group")

get_attr(efit, "contrast")

CHOL_DEGs <- limma::topTable(fit = efit, coef = 1, number = Inf)
setDT(CHOL_DEGs, keep.rownames = "rid")
head(CHOL_DEGs)

CHOL_DEGs[, let(
  tumor_vs_normal = fcase(
    adj.P.Val < 0.05 & logFC > 2, "Up",
    adj.P.Val < 0.05 & logFC < -2, "Down",
    default = "Not-Significant"
  )
)]

plotVolcano(
  data = CHOL_DEGs,
  x = logFC,
  y = -log10(adj.P.Val),
  color = tumor_vs_normal
)

# Over Presentation Analysis(ORA)
pathwayGO <- r4msigdb::query("Hs", pathway = "^GO(BP|CC|MF)_")
pathwayGO <- setNames(pathwayGO$symbol, pathwayGO$standard_name)

oraRes <- fgsea::fora(
  pathways = pathwayGO,
  genes = CHOL_DEGs[.("Up"), gene_name, on = .(tumor_vs_normal)],
  universe = unique(CHOL_DEGs$gene_name)
)
oraRes[, let(
  category = fcase(
    pathway %like% "GOBP", "BP",
    pathway %like% "GOMF", "MF",
    pathway %like% "GOCC", "CC"
  )
)]

oraRes <- oraRes[, .SD[order(padj)], by = .(category)]
oraRes[, let(pathwayGO = factor(pathway, levels = rev(pathway)))]
plotORA(
  data = oraRes[, head(.SD, 5), by = category],
  x = -log10(padj),
  y = pathwayGO,
  size = log10(overlap),
  fill = category
)
