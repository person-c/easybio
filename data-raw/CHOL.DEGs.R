# Gene expression aligned against hg38
setwd("data-raw")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)

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
## limma workflow
x <- dgeList(lt$all$exprCount, lt$all$sampleInfo, lt$all$featuresInfo)

x <- dprocess_dgeList(x, "group", 10)
efit <- limmaFit(x, "group")

CHOL.DEGs <- limma::topTable(fit = efit, coef = 1, number = Inf)
plotVolcano(data = degs, x = logFC, y = -log10(adj.P.Val))
usethis::use_data(CHOL.DEGs, compress = "xz", overwrite = TRUE)
