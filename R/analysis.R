#' A wrapper function to bio generic function
#'
#' Assign a task(class) to the data to execute task assignment.
#'
#' @param object data.
#' @param task task you want to execute.
#' @param ... additional arguments.
#'
#' @return Depending on the certain task.
#' @export
analyze <- function(object, task, ...) {
  class(object) <- c(task, class(object))
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

#' Use limma to do differential analysis
#'
#' limma is a R package used to do differential analysis for RNAseq and Array data.
#'
#' @param object A numeric data.frame with rownames.
#' @param pattern Regex pattern to divide the data into two groups.
#' @param data_type rnaSeq or array.
#' @param ... aditional arguments.
#'
#' @return list containing differential analysis and data
bio.limma <- function(object, pattern, data_type, ...) {
  index_control <- grep(pattern, colnames(object))
  index_treat <- grep(pattern, colnames(object), invert = TRUE)
  data <- object[, c(index_control, index_treat)]
  group <- rep(c("control", "treat"), c(length(index_control), length(index_treat)))

  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  rownames(design) <- colnames(data)
  contrast_matrix <- limma::makeContrasts(paste0(c("treat", "control"), collapse = "-"), levels = design)


  if (data_type == "rnaSeq") {
    data <- edgeR::DGEList(data)

    keep_exprs <- edgeR::filterByExpr(data, group)
    data <- data[keep_exprs, , keep.lib.sizes = FALSE]
    data <- edgeR::calcNormFactors(data, method = "TMM")

    v <- limma::voom(data, design, plot = TRUE)
    vfit <- limma::lmFit(v, design)
    vfit <- limma::contrasts.fit(vfit, contrasts = contrast_matrix)
    efit <- limma::eBayes(vfit)
    limma::plotSA(efit, main = "Final model: Mean-variance trend")

    result <- limma::topTable(efit, coef = 1, n = Inf)
    result <- list(diff = result, design_matrix = design, contrast = contrast_matrix)
    class(result) <- c("limma", class(result))

    return(result)
  }

  if (data_type == "array") {
    col_max <- sapply(X = object, FUN = max)
    if (any(col_max > 25)) warning("You'd better do log transformation before diff-analysis")

    fit <- limma::lmFit(limma::normalizeBetweenArrays(data, method = "cyclicloess"), design)
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    efit <- limma::eBayes(fit2)
    result <- limma::topTable(efit, coef = 1, n = Inf)

    result <- list(diff = result, design_matrix = design, contrast = contrast_matrix)

    class(result) <- c("limma", class(result))
    return(result)
  }
}

#' ClusterProfiler.
#'
#' GO(Gene ontology) enrichment analysis.
#'
#' @param object A gene vector used to do the enrichment analysis.
#' @param background Background gene set.
#' @param db The species information database.
#' @param from Original gene name task.
#' @param to Modified gene name task.
#' @param ... additional arguments
#'
#' @return list
bio.go <- function(
    object, background = NULL,
    db = org.Hs.eg.db, from = "SYMBOL", to = "ENSEMBL", ...) {
  dot_lists <- list(...)
  if (is.null(background)) {
    gene <- clusterProfiler::bitr(object, from, to, db, drop = TRUE)
    arg_lists <- list(gene = gene[[2]], OrgDb = db, keyType = to)
  } else {
    gene <- lapply(
      list(object, background),
      function(x) {
        clusterProfiler::bitr(x, from, to, db, drop = TRUE)
      }
    )

    arg_lists <- list(
      gene = gene[[1]][[2]],
      backgroud = gene[[2]][[2]],
      OrgDb = db,
      keytype = to
    )
  }

  args <- c(dot_lists, arg_lists)
  expr <- rlang::expr(clusterProfiler::enrichGO(!!!args))

  result <- eval(expr)
  class(result) <- "go"
  return(result)
}



#' use clusterProfiler
#'
#' kegg(Gene ontology) enrichment analysis
#'
#' @param object gene name used to do the enrichment analysis
#' @param backgroud the background gene set list
#' @param db the species information database
#' @param from original gene name task
#' @param to modified gene name task
#' @param ... additional arguments
#'
#' @return list
bio.kegg <- function(
    object, backgroud = NULL,
    db = org.Hs.eg.db, from = "SYMBOL", to = "ENTREZID", ...) {
  dot_lists <- list(...)
  if (is.null(backgroud)) {
    gene <- clusterProfiler::bitr(object, from, to, db, drop = TRUE)
    arg_lists <- list(gene = gene[[2]])
  } else {
    gene <- lapply(
      list(object, backgroud),
      function(x) {
        clusterProfiler::bitr(x, from, to, db, drop = TRUE)
      }
    )

    arg_lists <- list(gene = gene[[1]][[2]], backgroud = gene[[2]][[2]])
  }

  args <- c(dot_lists, arg_lists)
  expr <- rlang::expr(clusterProfiler::enrichKEGG(!!!args))

  result <- eval(expr)
  class(result) <- "kegg"
  return(result)
}


#' Use fgsea to do gsea enrich analysis
#'
#' GSEA
#' more on here
#' \url{https://www.bioconductor.org/packages/
#' devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html}
#' @param object Ranked gene vector according to the logFC.
#' @param  pathways List containing pathways and its associated genes.
#' @param ... additional arguments
#'
#' @return List
bio.gsea <- function(object, pathways, ...) {
  fgseaRes <- fgsea::fgsea(
    pathways = pathways,
    stats = object,
    minSize = 15,
    maxSize = 500
  )

  result <- list(fgseaRes, ranks = object, pathways = pathways)
  class(result) <- "gsea"
  return(result)
}


#' use surv to do survival analysis
#'
#' survival analysis
#' @param  object survival data
#' @param form survival form
#' @param ... additional parameters
#'
#' @return list
bio.surv <- function(object, form, ...) {
  arg <- rlang::enexpr(form)
  fit <- rlang::expr(survival::survfit(!!arg, object))
  fit <- eval(fit)

  result <- list(fit = fit, data = object)
  class(result) <- c("surv", class(result))

  return(result)
}

#' use surv to do cox analysis
#'
#' cox analysis
#' @param  object survival data
#' @param form cox form
#' @param ... additional parameters
#'
#' @return list
bio.cox <- function(object, form, ...) {
  arg <- rlang::enexpr(form)
  fit <- rlang::expr(survival::coxph(!!arg, object))
  fit <- eval(fit)

  result <- list(fit = fit, data = object)
  class(result) <- c("cox", class(result))
  return(result)
}
