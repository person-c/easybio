#' download GEO data
#'
#' A function to download GEO data and convert it to a standard expression if existed;
#' It will try to download the supplementary tabular data if no expression data is checked.
#'
#' @param geo GSE ID.
#' @param dir download directory.
#' @param combine whether to combine probes.
#' @param method ways to process many probes to one symbol(max or mean).
#' @param filter_regex regex keywords to decide what kind of supplementary file
#' to download.
#'
#' @return a list containing expression matrix and metadata.
#' @importFrom utils download.file
#' @import data.table
#' @export
download_geo <- function(geo, dir = ".", combine = TRUE, method = "max", filter_regex = NULL) {
  eset <- GEOquery::getGEO(GEO = geo, destdir = dir, getGPL = FALSE)
  if (length(eset) > 1) warning("There are more than one geo dataset;only the first one extracted")
  exp <- as.data.frame(eset[[1]]@assayData$exprs)
  pd <- eset[[1]]@phenoData@data


  if (nrow(exp) == 0L) {
    GEOquery::getGEOSuppFiles(geo)
    return(eset[[1]]@phenoData@data)
  }

  # if (nrow(exp) == 0L) {
  #   message(sprintf("%s seems not array data; Trying to download supplementary files.", geo))

  #   stub <- gsub("\\d{1,3}$", "nnn", geo, perl = TRUE)
  #   mirror <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", stub, geo)
  #   h <- rvest::read_html(mirror)
  #   nodes <- rvest::html_nodes(h, "a")
  #   a <- rvest::html_attr(nodes, "href")
  #   index <- ifelse(!is.null(filter_regex), grep(pattern = filter_regex, a), grep(pattern = ".gz|xls|csv", a))
  #   destfile_url <- paste0(mirror, a[index])
  #   Map(download.file, destfile_url, a[index])

  #   return(eset[[1]]@phenoData@data)
  # }

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
