#' @title Download and Process GEO Data
#'
#' @description
#' This function downloads gene expression data from the Gene Expression
#' Omnibus (GEO) database. It retrieves either the expression matrix or the
#' supplementary tabular data if the expression data is not available.
#' The function also allows for the conversion of probe identifiers to gene
#' symbols and can combine multiple probes into a single symbol.
#'
#' @param geo A character string specifying the GEO Series ID (e.g., "GSE12345").
#' @param dir A character string specifying the directory where files should be
#'   downloaded. Default is the current working directory (`"."`).
#' @param combine A logical value indicating whether to combine multiple probes
#'   into a single gene symbol. Default is `TRUE`.
#' @param method A character string specifying the method to use for combining
#'   probes into a single gene symbol. Options are `"max"` (take the maximum
#'   value) or `"mean"` (compute the average). Default is `"max"`.
#'
#' @return A list containing:
#' \item{data}{A data frame of the expression matrix, or `NULL` if not available.}
#' \item{sample}{A data frame of the sample metadata.}
#' \item{feature}{A data frame of the feature metadata, or `NULL` if not available.}
#' \item{status}{A character string indicating the data source:
#'   `"expression_matrix"`, `"supplementary_files"`, or `"no_data"`.}
#' \item{supplementary}{Only present when `status` is `"supplementary_files"`.
#'   A named list of `data.table` objects parsed from supplementary files.}
#'
#' @importFrom utils download.file
#' @export
prepare_geo <- function(geo, dir = ".", combine = TRUE, method = "max") {
  . <- ID <- symbol <- gene_assignment <- NULL # nolint: object_name_linter.

  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop(
      "To get GEO datasets, prepare_geo() requires 'GEOquery' package which ",
      "cannot be found. Please install 'GEOquery' using 'BiocManager::install('GEOquery')'."
    )
  }

  eset <- GEOquery::getGEO(GEO = geo, destdir = dir, getGPL = FALSE)
  if (length(eset) > 1) warning("There are more than one geo dataset;only the first one will be extracted")
  exp <- as.data.frame(eset[[1]]@assayData$exprs)
  pd <- eset[[1]]@phenoData@data


  if (nrow(exp) == 0L) {
    warning("No expression data is retrieved in series matrix; try to check the supplementary file")
    # code from GEOquery::getGEOSuppFiles()
    stub <- gsub("\\d{1,3}$", "nnn", geo, perl = TRUE)
    url <- sprintf(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/",
      stub, geo
    )

    fnames <- try(
      {
        a <- xml2::read_html(url)
        grep("^G", xml2::xml_text(xml2::xml_find_all(a, "//a/@href")),
          value = TRUE
        )
      },
      silent = TRUE
    )
    f_idx <- grep(pattern = "(count)|(fpkm)|(tpm)", x = fnames, ignore.case = TRUE)

    if (inherits(fnames, "try-error") || length(f_idx) == 0L) {
      message(sprintf("No potential expression data is detected in supplementary files"))
      message("Check URL manually if in doubt")
      message(url)

      return(list(
        data = NULL, sample = pd, feature = NULL,
        status = "no_data"
      ))
    }

    message("detect potential expression data: \n", paste0(fnames[f_idx], "\n"))
    message("read potential expression data in supplementary files...")
    res <- lapply(f_idx, \(idx) fread(paste0(url, fnames[[idx]])))
    names(res) <- make.names(fnames[[f_idx]])

    return(list(
      data = NULL, sample = pd, feature = NULL,
      status = "supplementary_files",
      supplementary = res
    ))
  }

  gpl <- GEOquery::getGEO(eset[[1]]@annotation, destdir = ".")
  gpl <- GEOquery::Table(gpl)
  setDT(gpl)
  if (!is.character(gpl[["ID"]])) {
    warning("The gpl annotation data's ID column is not character; Please check the gpl data carefully!")
    gpl[, let(ID = as.character(ID))]
  }
  gpl <- gpl[.(rownames(exp)), on = .(ID)]

  if (!combine) {
    gpl <- setDF(gpl, gpl[[1]])
    return(list(data = exp, sample = pd, feature = gpl, status = "expression_matrix"))
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
  return(list(data = exp2, sample = pd, feature = gpl2, status = "expression_matrix"))
}

#' Prepare TCGA Data for Analysis
#'
#' This function prepares TCGA data for downstream analyses such as
#' differential expression analysis with `limma` or survival analysis.
#' It extracts and processes the necessary information from the TCGA data
#' object, separating tumor and non-tumor samples.
#'
#' @param data A `SummarizedExperiment` object containing TCGA data, typically obtained from R package `TCGABiolinks`.
#'
#' @return A list.
#' @export
prepare_tcga <- function(data) {
  sample_info <- as.data.frame(data@colData)
  sample_info[["OS"]] <- fcoalesce(sample_info[["days_to_death"]], sample_info[["days_to_last_follow_up"]])

  features_info <- as.data.frame(data@rowRanges)
  rownames(features_info) <- data@rowRanges@ranges@NAMES
  expr <- as.data.frame(data@assays@data$unstranded, row.names = rownames(features_info))
  colnames(expr) <- rownames(sample_info)

  # tumor smaple data
  tumor_idx <- sample_info[["sample_type"]] %ilike% "Tumor"

  expr2 <- as.data.frame(data@assays@data$fpkm_unstrand, row.names = rownames(features_info))
  colnames(expr2) <- rownames(sample_info)
  expr2 <- expr2[, tumor_idx]

  sample_info2 <- sample_info[tumor_idx, ]

  structure(list(
    all = list(
      exprCount = expr,
      features_info = features_info,
      sample_info = sample_info
    ),
    tumor = list(
      exprFpkm = expr2,
      features_info = features_info,
      sample_info = sample_info2
    )
  ))
}
