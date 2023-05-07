#' a wrapper function to bio generic function
#'
#' assign a task(class) to the data to execute task assignment
#' to bio generic function
#'
#' @param object data used to analysis
#' @param type task you want to execute
#' @param ... additional arguments
#'
#' @return depending on the certain task
#' @export
#' @examples
#' # ADD_EXAMPLES_HERE
analysis <- function(object, type, ...) {
  class(object) <- c(type, class(object))
  bio(object, ...)

}

#' bio generic function
#'
#' assign task for bio generic function
#'
#' @param object data used to analysis
#' @param ... additional arguments for specific bio. function
#'
#' @return depending on the specific task
#' @examples
#' # ADD_EXAMPLES_HERE
bio <- function(object, ...) {
  UseMethod("bio")
}

#' use limma to do differential analysis
#'
#' do differential analysis accroding to smaple names
#'
#' @param data expression matrix
#' @param group_key key to distinguish group
#' @param data_type RNA-seq or array.
#' @param ... aditional arguments
#'
#' @return list containing differential analysis and data
#' @examples
#' # ADD_EXAMPLES_HERE
bio.limma <- function(data, group_key, data_type, ...) {

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
print(colnames(data))
group_list <- ifelse(
  grepl(group_key, colnames(data)) == TRUE, "control", "treat")

data <- data[, c(which(group_list == "control"), which(group_list == "treat"))]
group_list <- sort(group_list)

# design matrix
design <- stats::model.matrix(~0 + group_list)
colnames(design) <- gsub("group_list", "", colnames(design))
rownames(design) <- colnames(data)

contrast_matrix <- limma::makeContrasts(
  paste0(c("treat", "control"), collapse = "-"), levels = design)

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
  result <- list(diff = result, design_matrix = design,
    contrast = contrast_matrix, diff_input = data)
  class(result) <- c("limma", class(result))

  result

}

if (data_type == "array") {
  fit <- limma::lmFit(
    limma::normalizeBetweenArrays(data, method = "cyclicloess"), design)

  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  efit <- limma::eBayes(fit2)  ## default no trend

  result <- limma::topTable(efit, coef = 1, n = Inf)

  result <- list(diff = result, design_matrix = design,
      contrast = contrast_matrix, diff_input = data)
  class(result) <- c("limma", class(result))
  result

}
}



#' use clusterProfiler
#'
#' GO(Gene ontology) enrichment analysis
#'
#' @param .gene gene name used to do the enrichment analysis
#' @param .universe the background gene set list
#' @param .org_db the species information database
#' @param from_type original gene name type
#' @param to_type modified gene name type
#' @param ... additional arguments
#'
#' @importFrom clusterProfiler enrichGO
#'
#' @return list
#' @examples
#' # ADD_EXAMPLES_HERE
bio.go <- function(.gene, .universe = NULL,
  .org_db, from_type, to_type, ...) {

  dot_lists <- list(...) # nolint

  if (is.null(.universe)) {
    target_gene <- clusterProfiler::bitr(
      .gene, from_type, to_type, .org_db, drop = TRUE)

    arg_lists <- list(gene = target_gene[[2]],
    OrgDb = .org_db, keyType = to_type)
  } else {
    target_gene <- lapply(list(.gene, .universe),
    function(x) {
      clusterProfiler::bitr(x, from_type, to_type, .org_db, drop = TRUE)
    })

    arg_lists <- list(gene = target_gene[[1]][[2]],
    universe = target_gene[[2]][[2]],
    OrgDb = .org_db, keyType = to_type)
  }

  result <- do.call("enrichGO", args = c(dot_lists, arg_lists))
  class(result) <- "go"

  return(result)
}



#' use clusterProfiler
#'
#' kegg(Gene ontology) enrichment analysis
#'
#' @param .gene gene name used to do the enrichment analysis
#' @param .universe the background gene set list
#' @param .org_db the species information database
#' @param from_type original gene name type
#' @param to_type modified gene name type
#' @param ... additional arguments
#'
#' @importFrom clusterProfiler enrichKEGG
#'
#' @return list
#' @examples
#' # ADD_EXAMPLES_HERE
bio.kegg <- function(.gene, .universe = NULL,
  .org_db, from_type, to_type, ...) {

  dot_lists <- list(...)

  if (is.null(.universe)) {
    target_gene <- clusterProfiler::bitr(
      .gene, from_type, to_type, .org_db, drop = TRUE)

    arg_lists <- list(gene = target_gene[[2]])
  } else {
    target_gene <- lapply(list(.gene, .universe),
    function(x) {
      clusterProfiler::bitr(x, from_type, to_type, .org_db, drop = TRUE)
    })

    arg_lists <- list(gene = target_gene[[1]][[2]],
    universe = target_gene[[2]][[2]])
  }

  result <- do.call("enrichKEGG", args = c(dot_lists, arg_lists))
  class(result) <- "kegg"

  return(result)
}


#' use fgsea to do gsea enrich analysis
#'
#' gsea analysis
#' more on here
#' \url{https://www.bioconductor.org/packages/
#' devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html}
#' @param  pathways list containing pathways and its associated genes
#' @param ranks ranked gene list according to the logFC
#'
#' @return list
#' @examples
#' # ADD_EXAMPLES_HERE
bio.gsea <- function(pathways, ranks) {
  result <- fgsea::fgsea(pathways = pathways,
    stats = ranks,
    eps  = 0.0,
    minSize  = 15,
    maxSize  = 500)

  class(result) <- "gsea"
  return(result)
}


#' use surv to do survival analysis
#'
#' survival analysis
#' @param  data survival data
#' @param form survival form
#'
#' @return list
#' @examples
#' # ADD_EXAMPLES_HERE
bio.surv <- function(data, form) {
  force(data)
  arg <- rlang::enexpr(form)

  fit <- rlang::expr(survival::survfit(!!arg, data))
  result <- eval(fit)
  class(result) <- "surv"
  return(result)
}
