#' Calculate a linear model
#'
#' \code{calculate_lm} describes the linear model of the last listed variable
#' in \code{means} and \code{covs} as a function of all other variables in
#' \code{means} and \code{covs}.
#'
#' @param means a vector of means of all model predictors and the response with
#'   the last element the response mean.
#' @param covs a matrix of the covariance of all model predictors and the
#'   response with the order of rows/columns corresponding to the order of
#'   \code{means}.
#' @param n sample size
#' @param add_intercept logical. If \code{TRUE} adds an intercept to the model.
#' @param cl call
#' @param terms terms
#'
#' @importFrom stats pt
#'
#' @references Wolf, J.M., Barnard, M., Xueting, X., Ryder, N., Westra, J., and
#'   Tintle, N.  (2020). Computationally efficient, exact, covariate-adjusted
#'   genetic principal component analysis by leveraging individual marker
#'   summary statistics from large biobanks. \emph{Pacific Symposium on
#'   Biocomputing}, 25, 719-730.
#'
#'   Gasdaska A., Friend D., Chen R., Westra J., Zawistowski M., Lindsey W. and
#'   Tintle N. (2019) Leveraging summary statistics to make inferences about
#'   complex phenotypes in large biobanks. \emph{Pacific Symposium on
#'   Biocomputing}, 24, 391-402.
#'
calculate_lm <- function(means, covs, n, add_intercept = FALSE, cl = NULL, terms = NULL) {
  p <- ncol(covs) - 1
  if (add_intercept) {
    p <- p + 1
    covs <- rbind(0, cbind(0, covs))
    means <- c(1, means)
    if (!is.null(dimnames(covs)) & !is.null(dimnames(covs))) {
      names(means)[1] <- "(Intercept)"
      colnames(covs)[1] <- "(Intercept)"
      rownames(covs)[1] <- "(Intercept)"
    }
  }

  covX <- covs[-(p + 1), -(p + 1)]
  meanX <- means[-(p + 1)]
  XtX <- (n - 1) * covX + n * meanX %*% t(meanX)
  Xty <- (n - 1) * covs[p + 1, -(p + 1)] + n * meanX * means[p + 1]
  yty <- (n - 1) * covs[p + 1, p + 1] + n * means[p + 1]^2

  beta <- drop(solve(XtX) %*% Xty)

  sigma2 <- drop((yty - t(beta) %*% Xty) / (n - p))

  var_beta <- diag(sigma2 * solve(XtX))
  sd_beta <- sqrt(var_beta)

  t_stat <- beta / sd_beta
  p_val <- 2 * pt(abs(t_stat), df = n - p, lower.tail = F)

  coefficients <- cbind(beta, sd_beta, t_stat, p_val)
  colnames(coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

  aliased <- rep(FALSE, p)
  names(aliased) <- rownames(coefficients)

  sigma <- sqrt(sigma2)

  df <- c(p, n - p, p)


  SSE <- drop(yty - t(Xty) %*% solve(XtX) %*% Xty)
  SST <- (n - 1) * covs[p + 1, p + 1]
  SSR <- SST - SSE

  SS <- c(SSR = SSR, SSE = SSE, SST = SST)

  MSR <- SSR / (p - 1)
  MSE <- SSE / (n - p)

  r.squared <- 1 - (SSE / SST)
  adj.r.squared <- 1 - (1 - r.squared) * (n - 1) / (n - p)

  fstatistic <- c(value = MSR / MSE, numdf = p - 1, numdf = n - p)

  cov.unscaled <- solve(XtX)

  re <- list(
    call = cl, terms = terms, residuals = NULL,
    coefficients = coefficients,
    aliased = aliased,
    sigma = sigma,
    df = df,
    r.squared = r.squared,
    adj.r.squared = adj.r.squared,
    fstatistic = fstatistic,
    cov.unscaled = cov.unscaled,
    `Sum Sq` = SS
  )
  
  class(re) <- "pcsslm"

  return(re)
}
