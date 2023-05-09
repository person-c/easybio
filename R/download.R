#' download GEO data
#'
#' A function to download GEO data and convert it to a standard expression
#' if existed; It will try to download the supplementary tabular data if no
#' expression data checked.
#'
#' @param geo GSE ID.
#' @param dir download directory.
#' @param method ways to process many probes to one symbol(max or mean).
#' @param filter_regex regex keywords to decide what kind of supplementary file
#' to download.
#'
#' @return a list containing expression matrix and metadata.
#' @importFrom utils download.file
#' @importFrom stats na.omit
#' @export
#' @examples
#' # ADD_EXAMPLES_HERE
download_geo <- function(geo, dir, method = "max", filter_regex = NULL) {

eset <- GEOquery::getGEO(GEO = geo, destdir = dir, getGPL = FALSE)
exp <- eset[[1]]@assayData$exprs

# check whether the GSE data are high throught data
if (dim(exp)[[1]] == 0) {
  message(sprintf("%s seems not array data; Trying to download supplementary files.", geo))

  stub <- gsub("\\d{1,3}$", "nnn", geo, perl = TRUE)
  mirror <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", stub, geo)
  h <- rvest::read_html(mirror)

  nodes <- rvest::html_nodes(h, "a")
  a <- rvest::html_attr(nodes, "href")

  if (!is.null(filter_regex)) {
    index <- which(grepl(pattern = filter_regex, a) == TRUE)
  } else {
    index <- which(grepl(pattern = ".gz|xls|csv", a) == TRUE)
  }

  destfile_url <- paste(mirror, a[index], sep = "")
  purrr::walk2(destfile_url, a[index], download.file)

  return(eset[[1]]@phenoData@data)
}
# download GPL annotation and convert probe_id to symbol
gpl_annotation <- GEOquery::getGEO(eset[[1]]@annotation, destdir = ".")
gpl_annotation <- GEOquery::Table(gpl_annotation)

if (purrr::some( #check gene_assignment;if existence, get symbol from assignment
      colnames(gpl_annotation),
      ~ stringr::str_detect(string = .x, pattern = "gene_assignment"))) {

    gene_symbol <- strsplit(x = gpl_annotation$gene_assignment, "//")
    gene_symbol <- purrr::map_chr(gene_symbol, ~ `[`(.x, 2))
    ids <- data.frame(probe_id = gpl_annotation[[1]], symbol = gene_symbol)
    ids <- na.omit(ids)

  } else {
    if (purrr::some( # get symbol column
      colnames(gpl_annotation),
      ~ stringr::str_detect(string = .x, pattern = "[Ss]ymbol"))) {

      ids <- gpl_annotation[,
        c(1, which(stringr::str_detect(colnames(gpl_annotation),
          pattern = "[Ss]ymbol") == TRUE))]
      ids <- na.omit(ids)
      names(ids) <- c("probe_id", "symbol")

      } else {
        warning(sprintf("%s: Thers are no columns containing string [Ss]ymbol
        in the GPL annotation; try to check the GPL annotation's column names",
        geo)
        )

        return(list(exprMatrix = exp, metaData = eset[[1]]@phenoData@data))
      }
  }

exp <- exp[rownames(exp) %in% ids$probe_id, ]
ids <- ids[match(rownames(exp), ids$probe_id), ]

# a symbol corresponds to many probe
# choose the max row mean probe or get the mean vale of these probe
if (method == "max") {
  tmp <- by(exp, ids$symbol,
    function(x) rownames(x)[which.max(rowMeans(x))])
  probes <- as.character(tmp)
  exp <-  exp[rownames(exp) %in% probes, ]
  rownames(exp) <- ids[match(rownames(exp), ids$probe_id), 2]
  # get sample informaiton
  pd <- eset[[1]]@phenoData@data

  return(list(exprMatrix = exp, metaData = pd))
}

if (method == "mean") {
  tmp <- by(exp, ids$symbol, function(x) colSums(x) / nrow(x))

  exp <- rlang::exec("rbind", !!!tmp)
  pd <- eset[[1]]@phenoData@data

  return(list(exprMatrix = exp, metaData = pd))
}
}
