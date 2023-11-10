#' download GEO data
#'
#' A function to download GEO data and convert it to a standard expression if existed;
#' It will try to download the supplementary tabular data if no expression data is checked.
#'
#' @param geo GSE ID.
#' @param dir download directory.
#' @param method ways to process many probes to one symbol(max or mean).
#' @param filter_regex regex keywords to decide what kind of supplementary file
#' to download.
#'
#' @return a list containing expression matrix and metadata.
#' @importFrom utils download.file
#' @export
#' @examples
#' # ADD_EXAMPLES_HERE
download_geo <- function(geo, dir = ".", method = "max", filter_regex = NULL) {
  eset <- GEOquery::getGEO(GEO = geo, destdir = dir, getGPL = FALSE)
  exp <- as.data.frame(eset[[1]]@assayData$exprs)
  pd <- eset[[1]]@phenoData@data

  if (length(exp) == 0L) {
    message(sprintf("%s seems not array data; Trying to download supplementary files.", geo))

    stub <- gsub("\\d{1,3}$", "nnn", geo, perl = TRUE)
    mirror <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", stub, geo)
    h <- rvest::read_html(mirror)
    nodes <- rvest::html_nodes(h, "a")
    a <- rvest::html_attr(nodes, "href")
    index <- ifelse(!is.null(filter_regex), grep(pattern = filter_regex, a), grep(pattern = ".gz|xls|csv", a))
    destfile_url <- paste0(mirror, a[index])
    purrr::walk2(destfile_url, a[index], download.file)

    return(eset[[1]]@phenoData@data)
  }

  gpl <- GEOquery::getGEO(eset[[1]]@annotation, destdir = ".")
  gpl <- GEOquery::Table(gpl)

  if ("gene_assignment" %in% colnames(gpl)) {
    gene_symbol <- strsplit(x = gpl$gene_assignment, "//")
    gene_symbol <- purrr::map_chr(gene_symbol, ~ `[`(.x, 2))
    ids <- data.frame(probe_id = gpl[[1]], symbol = gene_symbol)
    ids <- stats::na.omit(ids)
  } else if (any(grepl(x = colnames(gpl), pattern = "symbol|genename", ignore.case = TRUE))) {
    ids <- gpl[, c(1, grep(colnames(gpl), pattern = "symbol|genename", ignore.case = TRUE))]
    ids <- stats::na.omit(ids)
    names(ids) <- c("probe_id", "symbol")
  } else {
    exp <- data.table::setDT(exp, keep.rownames = TRUE)
    warning(sprintf(
      "%s: It seems that no columns containing symbol information
      in the GPL annotation; try to check the GPL annotation's column names",
      geo
    ))

    return(list(exprMatrix = exp, metaData = pd))
  }

  ids <- data.table::setDT(ids)
  exp <- data.table::setDT(exp, keep.rownames = TRUE)
  exp <- merge(ids, exp, by.x = "probe_id", by.y = "rn")

  if (method == "max") {
    exp <- exp[, .SD[which.max(rowMeans(.SD, na.rm = TRUE))], keyby = symbol, .SDcols = is.numeric]
  }
  if (method == "mean") {
    exp <- exp[, lapply(.SD, function(x) sum(x) / length(x)), keyby = symbol, .SDcols = is.numeric]
  }

  exp <- exp[, symbol := gsub("\\s.+", "", symbol)][!grepl("^$", symbol)]

  return(list(exprMatrix = exp, metaData = pd))
}
