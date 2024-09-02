#' @title Download and Process GEO Data
#'
#' @description
#' This function downloads gene expression data from the Gene Expression Omnibus (GEO) database.
#' It retrieves either the expression matrix or the supplementary tabular data if the expression data is not available.
#' The function also allows for the conversion of probe identifiers to gene symbols and can combine multiple probes into a single symbol.
#'
#' @param geo A character string specifying the GEO Series ID (e.g., "GSE12345").
#' @param dir A character string specifying the directory where files should be downloaded. Default is the current working directory (`"."`).
#' @param combine A logical value indicating whether to combine multiple probes into a single gene symbol. Default is `TRUE`.
#' @param method A character string specifying the method to use for combining probes into a single gene symbol. Options are `"max"` (take the maximum value) or `"mean"` (compute the average). Default is `"max"`.
#' @param filter_regex A character string specifying a regular expression to filter the supplementary files for download. Default is `NULL`, which means no filtering is applied.
#'
#' @return A list containing:
#' \item{data}{A data frame of the expression matrix.}
#' \item{sample}{A data frame of the sample metadata.}
#' \item{feature}{A data frame of the feature metadata, which includes gene symbols if combining probes.}
#'
#' @importFrom utils download.file
#' @import data.table
#' @export
prepare_geo <- function(geo, dir = ".", combine = TRUE, method = "max", filter_regex = NULL) {
  . <- ID <- symbol <- gene_assignment <- NULL

  eset <- GEOquery::getGEO(GEO = geo, destdir = dir, getGPL = FALSE)
  if (length(eset) > 1) warning("There are more than one geo dataset;only the first one extracted")
  exp <- as.data.frame(eset[[1]]@assayData$exprs)
  pd <- eset[[1]]@phenoData@data


  if (nrow(exp) == 0L) {
    GEOquery::getGEOSuppFiles(geo)
    return(eset[[1]]@phenoData@data)
  }

  gpl <- GEOquery::getGEO(eset[[1]]@annotation, destdir = ".")
  gpl <- GEOquery::Table(gpl)
  setDT(gpl)
  gpl <- gpl[.(rownames(exp)), on = .(ID)]

  if (!combine) {
    gpl <- setDF(gpl, gpl[[1]])
    return(list(data = exp, sample = pd, feature = gpl))
  }

  gpl2 <- copy(gpl)
  if ("gene_assignment" %chin% colnames(gpl)) {
    gpl2[, symbol := sapply(strsplit(x = gene_assignment, "//"), `[`, 2)]
  }

  if (any(colnames(gpl) %ilike% "symbol|genename")) {
    gpl2[["symbol"]] <- gpl[[which(colnames(gpl) %ilike% "symbol|genename")]]
  }

  exp2 <- as.data.table(exp)
  exp2[, symbol := gpl2[, symbol]]
  exp2 <- exp2[symbol != ""]

  if (method == "max") {
    exp2 <- exp2[, .SD[which.max(rowMeans(.SD, na.rm = TRUE))], by = symbol, .SDcols = is.numeric]
  }
  if (method == "mean") {
    exp2 <- exp2[, lapply(.SD, function(x) sum(x) / length(x)), by = symbol, .SDcols = is.numeric]
  }

  exp2 <- setDF(exp2, exp2$symbol)
  exp2$symbol <- NULL
  gpl2 <- gpl2[.(rownames(exp2)), on = .(symbol), mult = "first"]
  gpl2 <- setDF(gpl2, gpl2$symbol)
  gpl2$symbol <- NULL
  return(list(data = exp2, sample = pd, feature = gpl2))
}

#' Prepare TCGA Data for Analysis
#'
#' This function prepares TCGA data for downstream analyses such as differential expression analysis with `limma` or survival analysis.
#' It extracts and processes the necessary information from the TCGA data object, separating tumor and non-tumor samples.
#'
#' @param data A `SummarizedExperiment` object containing TCGA data, typically obtained from R package `TCGABiolinks`.
#'
#' @return A list.
#' @import data.table
#' @export
prepare_tcga <- function(data) {
  sampleInfo <- as.data.frame(data@colData)
  sampleInfo[["OS"]] <- fcoalesce(sampleInfo[["days_to_death"]], sampleInfo[["days_to_last_follow_up"]])

  featuresInfo <- as.data.frame(data@rowRanges)
  expr <- as.data.frame(data@assays@data$unstranded, row.names = rownames(featuresInfo))
  colnames(expr) <- rownames(sampleInfo)

  # tumor smaple data
  tumorIdx <- sampleInfo[["sample_type"]] %ilike% "Tumor"

  expr2 <- as.data.frame(data@assays@data$fpkm_unstrand, row.names = rownames(featuresInfo))
  colnames(expr2) <- rownames(sampleInfo)
  expr2 <- expr2[, tumorIdx]

  sampleInfo2 <- sampleInfo[, tumorIdx]

  structure(list(
    all = list(
      exprCount = expr,
      featuresInfo = featuresInfo,
      sampleInfo = sampleInfo
    ),
    tumor = list(
      exprFpkm = expr2,
      featuresInfo = featuresInfo,
      sampleInfo = sampleInfo2
    )
  ))
}
