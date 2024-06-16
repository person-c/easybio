#' CoreAlg performs nu-regression using support vector machines (SVM) and calculates weights, root mean squared error (RMSE), and correlation coefficient (R).
#'
#' @param X Input matrix or data frame containing the predictor variables.
#' @param y Numeric vector containing the response variable.
#' @param absolute Logical value indicating whether to use absolute space or relative space for the weights.
#' @param abs_method String indicating the method to calculate the weights in absolute space. Can be either 'sig.score' or 'no.sumto1'.
CoreAlg <- function(X, y, absolute, abs_method) {
  # try different values of nu
  svn_itor <- 3

  res <- function(i) {
    if (i == 1) {
      nus <- 0.25
    }
    if (i == 2) {
      nus <- 0.5
    }
    if (i == 3) {
      nus <- 0.75
    }
    model <- svm(X, y, type = "nu-regression", kernel = "linear", nu = nus, scale = F)
    model
  }

  if (Sys.info()["sysname"] == "Windows") {
    out <- mclapply(1:svn_itor, res, mc.cores = 1)
  } else {
    out <- mclapply(1:svn_itor, res, mc.cores = svn_itor)
  }

  nusvm <- rep(0, svn_itor)
  corrv <- rep(0, svn_itor)

  # do cibersort
  t <- 1
  while (t <= svn_itor) {
    weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights < 0)] <- 0
    w <- weights / sum(weights)
    u <- sweep(X, MARGIN = 2, w, "*")
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  # pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  # get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q < 0)] <- 0
  if (!absolute || abs_method == "sig.score") w <- (q / sum(q)) # relative space (returns fractions)
  if (absolute && abs_method == "no.sumto1") w <- q # absolute space (returns scores)

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
}

#' doPerm performs permutation-based sampling and runs the CoreAlg function iteratively.
#'
#' @param perm Number of permutations to perform.
#' @param X Input matrix or data frame containing the predictor variables.
#' @param Y Numeric vector containing the response variable.
#' @param absolute Logical value indicating whether to use absolute space or relative space for the weights.
#' @param abs_method String indicating the method to calculate the weights in absolute space. Can be either 'sig.score' or 'no.sumto1'.
doPerm <- function(perm, X, Y, absolute = FALSE, abs_method = "sig.score") {
  itor <- 1
  ylist <- as.list(data.matrix(Y))
  dist <- numeric(length = perm)

  while (itor <= perm) {
    yr <- as.numeric(ylist[sample(length(ylist), nrow(X))])
    yr <- (yr - mean(yr)) / sd(yr)
    result <- CoreAlg(X, yr, absolute, abs_method)
    mix_r <- result$mix_r

    dist[[itor]] <- mix_r

    itor <- itor + 1
  }
  res <- list(dist = dist)
  res
}

#' Estimation of the abundances of member cell types.
#'
#' @param sig_matrix  Cell type GEP barcode matrix: row 1 = sample labels; column 1 = gene symbols; no missing values;
#' @param mixture_file  GEP matrix: row 1 = sample labels; column 1 = gene symbols; no missing values
#' @param perm Set permutations for statistical analysis (above 100 permutations recommended).
#' @param QN Quantile normalization of input mixture (default = TRUE)
#' @param absolute  Run CIBERSORT in absolute mode (default = FALSE)
#' @param abs_method  if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'
#' @author Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' @return cibersrot with immune cell fractions
#' @export
CIBERSORT <- function(
    sig_matrix, mixture_file, perm, QN = FALSE, absolute = FALSE,
    abs_method = c("sig.score", "no.sumto1")) {
  if (length(intersect(rownames(mixture_file), rownames(sig_matrix))) == 0L) {
    stop("None identical gene between eset and reference had been found.")
  }
  if (absolute) match.arg(abs_method)


  X <- data.matrix(sig_matrix)
  Y <- data.matrix(mixture_file)


  # order
  X <- X[order(rownames(X)), ]
  Y <- Y[order(rownames(Y)), ]

  P <- perm # number of permutations

  # anti-log if max < 50 in mixture file
  if (max(Y) < 50) {
    Y <- 2^Y
  }

  # quantile normalization of mixture file
  if (QN) {
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- limma::normalizeBetweenArrays(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  # store original mixtures
  Yorig <- Y
  Ymedian <- max(median(Yorig), 1)

  # intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX, ]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY, ]

  # standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  # empirical null distribution of correlation coefficients
  if (P > 0) {
    nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)
  }

  header <- c("Mixture", colnames(X), "P-value", "Correlation", "RMSE")
  if (absolute) header <- c(header, paste("Absolute score (", abs_method, ")", sep = ""))

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  # iterate through mixtures
  while (itor <= mixtures) {
    y <- Y[, itor]

    # standardize mixture
    y <- (y - mean(y)) / sd(y)

    # run SVR core algorithm
    result <- CoreAlg(X, y, absolute, abs_method)

    # get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    if (absolute && abs_method == "sig.score") {
      w <- w * median(Y[, itor]) / Ymedian
    }

    # calculate p-value
    if (P > 0) {
      pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
    }

    # print output
    out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
    if (absolute) out <- c(out, sum(w))
    if (itor == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }

    itor <- itor + 1
  }

  obj <- rbind(header, output)
  obj <- obj[-1, -1]
  obj <- matrix(as.numeric(unlist(obj)), nrow = nrow(obj))
  rownames(obj) <- colnames(Y)
  if (!absolute) {
    colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE")
  } else {
    colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE", paste("Absolute score (", abs_method, ")", sep = ""))
  }
  obj
}
