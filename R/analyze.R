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
#' limma is a R package used to do differential analysis
#' for RNAseq and Array data.
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
  group <- rep(
    c("control", "treat"),
    c(length(index_control), length(index_treat))
  )

  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  rownames(design) <- colnames(data)
  contrast_matrix <- limma::makeContrasts(
    paste0(c("treat", "control"), collapse = "-"),
    levels = design
  )


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
    result <- list(
      diff = result,
      design_matrix = design,
      contrast = contrast_matrix
    )
    class(result) <- c("limma", class(result))

    return(result)
  }

  if (data_type == "array") {
    if (any(sapply(object, max) > 20)) {
      warning("You'd better do log transformation before diff-analysis")
    }

    fit <- limma::lmFit(data, design)
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    efit <- limma::eBayes(fit2)
    result <- limma::topTable(efit, coef = 1, n = Inf)

    result <- list(
      diff = result,
      design_matrix = design,
      contrast = contrast_matrix
    )

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
  fgseaRes <- clusterProfiler::GSEA(object, TERM2GENE = pathways, ...)
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

#' WGCNA - gene cluster analysis
#'
#' Use WGCNA to do gene cluster analysis
#' @param object a numeric matrix or data frame,
#' with genes(feature) in columns and samples in rows.
#' @param trait a numeric matrix or data frame with
#' clinical features in columns.
#' All classification features should be converted to numeric.
#' @param  powers a numeric atomic vector to be used for the powers pick.
#' @param power the picked power from powers.
#' @param ... additional arguments
#'
#' @return depending on the arguments you supply.

bio.wgcna <- function(
    object,
    trait = NULL,
    powers = NULL,
    power = NULL,
    ...) {
  # check missings
  gsg <- WGCNA::goodSamplesGenes(object, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
      dynamicTreeCut::printFlush(
        paste(
          "Removing genes:",
          paste(names(object)[!gsg$goodGenes], collapse = ", ")
        )
      )
    }
    if (sum(!gsg$goodSamples) > 0) {
      dynamicTreeCut::printFlush(
        paste(
          "Removing samples:",
          paste(rownames(object)[!gsg$goodSamples], collapse = ", ")
        )
      )
    }
    object <- object[gsg$goodSamples, gsg$goodGenes]
  }

  nsamples <- nrow(object)
  # check outliers of samples
  sptree <- fastcluster::hclust(dist(object), method = "average")

  if (is.null(trait)) {
    png("sample-cluster.png")
    WGCNA::sizeGrWindow(12, 9)
    par(cex = 0.6)
    par(mar = c(0, 4, 2, 0))
    plot(sptree,
      main = "Sample clustering to detect outliers",
      sub = "", xlab = "", cex.lab = 1.5,
      cex.axis = 1.5, cex.main = 2
    )
  } else {
    trait2colors <- WGCNA::numbers2colors(trait, signed = FALSE)
    # Plot the sample dendrogram and the colors underneath.
    WGCNA::sizeGrWindow(12, 9)
    p <- WGCNA::plotDendroAndColors(sptree, trait2colors,
      groupLabels = names(trait),
      main = "Sample dendrogram and trait heatmap"
    )

    p
  }


  # pick soft-traitColors
  if (!is.null(powers)) {
    WGCNA::enableWGCNAThreads()
    sft <- WGCNA::pickSoftThreshold(object, powerVector = powers, verbose = 5)
    # Plot the results:
    png("power-tunning.png")
    WGCNA::sizeGrWindow(9, 5)
    par(mfrow = c(1, 2))
    cex1 <- 0.9
    # Scale-free topology fit index as a function of the soft-thresholding power
    plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
      xlab = "Soft Threshold (power)",
      ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
      main = paste("Scale independence")
    )
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
      labels = powers, cex = cex1, col = "red"
    )
    # this line corresponds to using an R^2 cut-off of h
    abline(h = 0.90, col = "red")
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
      xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
      main = paste("Mean connectivity")
    )
    text(sft$fitIndices[, 1], sft$fitIndices[, 5],
      labels = powers, cex = cex1, col = "red"
    )
    return(sft)
  }

  # network construction based on topologial overlap matrix
  if (!is.null(power)) {
    cor <- WGCNA::cor
    WGCNA::enableWGCNAThreads()
    net <- WGCNA::blockwiseModules(object,
      power = power,
      TOMType = "unsigned", minModuleSize = 30,
      reassignThreshold = 0, mergeCutHeight = 0.25,
      numericLabels = TRUE, pamRespectsDendro = FALSE,
      maxBlockSize = 60000,
      saveTOMs = TRUE,
      verbose = 3
    )
    cor <- stats::cor

    mcolors <- WGCNA::labels2colors(net$colors)
    p <- WGCNA::plotDendroAndColors(net$dendrograms[[1]],
      mcolors[net$blockGenes[[1]]],
      "Module colors",
      dendroLabels = FALSE, hang = 0.03,
      addGuide = TRUE, guideHang = 0.05
    )
    print(p)
    gcluster <- data.table::data.table(
      gene = names(net$colors),
      colors = WGCNA::labels2colors(net$colors)
    )


    meg <- WGCNA::moduleEigengenes(object, mcolors)$eigengenes
    meg <- WGCNA::orderMEs(meg)
    data.table::fwrite(meg, "eigenvector.csv")
    data.table::fwrite(x = gcluster, "module.csv")

    if (!is.null(trait)) {
      cor5mt <- WGCNA::cor(meg, trait, use = "p")
      p5cor5mt <- WGCNA::corPvalueStudent(cor5mt, nsamples)

      # view module-traits relationship heatmap
      tmatrix <- paste(
        signif(cor5mt, 2),
        "\n(", signif(p5cor5mt, 1), ")",
        sep = ""
      )
      dim(tmatrix) <- dim(cor5mt)
      par(mar = c(6, 8.5, 3, 3))
      # Display the correlation values within a heatmap plot
      p <- WGCNA::labeledHeatmap(
        Matrix = cor5mt,
        xLabels = names(trait),
        yLabels = names(meg),
        ySymbols = names(meg),
        colorLabels = FALSE,
        colors = WGCNA::blueWhiteRed(50),
        textMatrix = tmatrix,
        setStdMargins = FALSE,
        cex.text = 0.5,
        zlim = c(-1, 1),
        main = paste("Module-trait relationships")
      )

      print(p)
    }

    return(net)
  }
}


#' Use STRINGdb to do PPI analysis
#'
#' PPI analysis
#' @param object a data.frame with one column of gene
#' @param ... additional arguments
#'
#' @return data.table
bio.ppi <- function(object, ...) {
  dt <- object
  st <- STRINGdb::STRINGdb$new(
    version = "12.0",
    species = 9606, score_threshold = 400
  )

  stid <- st$map(dt, colnames(dt)[[1]], removeUnmappedRows = TRUE)
  st$plot_network(stid$STRING_id)

  inter <- st$get_interactions(stid$STRING_id)
  id2symbol <- stid$candidant
  names(id2symbol) <- stid$STRING_id

  data.table::setDT(inter)[, from := id2symbol[from]][, to := id2symbol[to]]
  data.table::setnames(inter, old = "combined_score", "weight")

  class(inter) <- c("ppi", class(inter))
}
