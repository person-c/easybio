#' use limma to do differential analysis
#'
#' do differential analysis accroding to smaple names
#'
#' @param data expression matrix
#' @param group_key key to distinguish group
#' @param data_type RNA-seq or array.
#'
#' @return list containing differential analysis and data
#' @export
#' @examples
#' # ADD_EXAMPLES_HERE
analysis_limma <- function(data, group_key, data_type) {

# array data check
if (data_type == "array") {
  data_max <-  summary(data)[6, ] |>
  as.character() |>
  strsplit(split = ":") |>
  purrr::map_chr(~ `[`(.x, 2)) |>
  as.numeric()

if (purrr::some(data_max, ~ .x > 100)) {
  warning("You'd better do log transformation before diff-analysis")
}
}

# group accroding to regex match
group_list <- ifelse(
  grepl(group_key, colnames(data)) == TRUE, "control", "treat")

data <- data[, c(which(group_list == "control"), which(group_list == "treat"))]
group_list <- sort(group_list)

# design matrix
design <- stats::model.matrix(~0 + group_list)
colnames(design) <- gsub("group_list", "", colnames(design))
rownames(design) <- colnames(data)
design

# contrast matrix
contrast_matrix <- limma::makeContrasts(
  paste0(c("treat", "control"), collapse = "-"), levels = design)
contrast_matrix

# RNA-Seq diff

if (data_type == "RNA-Seq") {
  data <- edgeR::DGEList(data)

  keep_exprs <- edgeR::filterByExpr(data, group = group_list)
  data <- data[keep_exprs, , keep.lib.sizes = FALSE]
  data <- edgeR::calcNormFactors(data, method = "TMM")

  v <- limma::voom(data, design, plot = TRUE)
  vfit <- limma::lmFit(v, design)
  vfit <- limma::contrasts.fit(vfit, contrasts = contrast_matrix)
  efit <- limma::eBayes(vfit)
  limma::plotSA(efit, main = "Final model: Mean-variance trend")

  result <- limma::topTable(efit, coef = 1, n = Inf)
  return(
    result <- list(diff = result, design_matrix = design,
      contrast = contrast_matrix, diff_input = data)
    )
    class(result) <- "limma"
}

if (data_type == "array") {
  fit <- limma::lmFit(
    limma::normalizeBetweenArrays(data, method = "cyclicloess"), design)

  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  efit <- limma::eBayes(fit2)  ## default no trend

  result <- limma::topTable(efit, coef = 1, n = Inf)

  return(
    result <- list(diff = result, design_matrix = design,
      contrast = contrast_matrix, diff_input = data)
    )
    class(result) <- "limma"
}
}
