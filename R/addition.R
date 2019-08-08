## Functions to add linear combinations of phenotypes


#' Adjust a linear model for covariates
#' @param coefs 2 x p Matrix of model coefficients for y ~ x_1, ... y ~ x_p
#'   The first and second rows correspond to intercept and slope coefficients
#'   for y ~ x, respectively.
#' @param xMeans Vector of length p of means of all covariates
#' @param xCov Covariance matrix of all covariates
#' @param yMean Response mean
#' @param yVar Response variance
#' @param n Sample size
#'
#' @examples
#' data(mtcars)
#'
#' dispMod <- lm(mpg ~ disp, data = mtcars)
#' cylMod  <- lm(mpg ~ cyl,  data = mtcars)
#'
#' coefs <- cbind(coefficients(dispMod), coefficients(cylMod))
#'
#' xMeans <- c(mean(mtcars$disp), mean(mtcars$cyl))
#'
#' xCov <- cov(mtcars)[c('disp', 'cyl'), c('disp', 'cyl')]
#'
#' yMean <- mean(mtcars$mpg)
#' yVar <- var(mtcars$mpg)
#'
#' n <- nrow(mtcars)
#'
#' out <- lm_adjust(coefs = coefs,
#'                  xMeans = xMeans,
#'                  xCov = xCov,
#'                  yMean = yMean,
#'                  yVar = yVar,
#'                  n = n)
#'
#' out
#' summary(lm(mpg ~ disp + cyl, data = mtcars))$coefficients


lm_adjust <- function(coefs,
                      xMeans,
                      xCov,
                      yMean,
                      yVar,
                      n) {
  ## nCovar := Number of covariates. (Including intercept)
  nCovar <- ncol(coefs) + 1
  
  ## xMeans := column means of design matrix
  xMeans <- c(1, xMeans)
  
  ## varX := Covariance of design matrix
  varX <- cbind(rep(0, nCovar),
                rbind(rep(0, nCovar - 1),
                      xCov))
  ## XtX := X'X
  ## To do: re-write code to avoid loop and use matrix algebra
  XtX <- matrix(NA, nrow = nrow(varX), ncol = ncol(varX))
  for (i in 1:nrow(XtX)) {
    for (j in 1:nrow(XtX)) {
      XtX[i, j] <- varX[i, j] * (n - 1) + xMeans[i] * xMeans[j] * n
    }
  }
  
  ## XtY := X'y
  XtY <- matrix(NA, nrow = nCovar,  ncol = 1)
  ## To do: re-write code to avoid loop and use matrix algebra
  XtY[1, 1] <- yMean * n
  for (i in 2:nCovar) {
    XtY[i, 1] <-
      coefs[2, i - 1] * varX[i, i] * (n - 1) + xMeans[i] * yMean * n
  }
  
  coefs <- solve(XtX) %*% XtY
  
  YtY <- yVar * (n - 1) + yMean ^ 2 * n
  
  varCoefs <-
    as.double((YtY - t(coefs) %*% XtX %*% coefs) / (n - nCovar)) * solve(XtX)
  
  seCoefs <- sqrt(diag(varCoefs))
  
  out <- list()
  out$coefs <- coefs
  out$seCoefs <- seCoefs
  return(out)
}

#' Model a linear combination of phenotypes
#'
#' @param coefs (p + 1) * m matrix of model coefficients for
#'   y_1 ~ x_1 + ... + x_p, ..., y_m ~ x_1 + ... + x_p.
#'   Columns correspond to responses and rows correspond to covariates.
#' @param xMeans Vector of length p of all means of covariates
#' @param xCov Covariance matrix of all covariates
#' @param yMeans Vector of length m of all means of responses
#' @param yCov Covariance matrix of all responses
#' @param n Sample size
#' @param weights Vector of length m of phenotype weights. Defaults to
#'   equal weights.
#' @examples 
#' mpgMod <- lm(mpg ~ disp + cyl, data = mtcars)
#' hpMod  <- lm(hp  ~ disp + cyl, data = mtcars)
#' 
#' coefs <- cbind(coefficients(mpgMod), coefficients(hpMod))
#' xMeans <- c(mean(mtcars$disp), mean(mtcars$cyl))
#' xCov <- cov(mtcars)[c('disp', 'cyl'), c('disp', 'cyl')]
#' 
#' yMeans <- c(mean(mtcars$mpg), mean(mtcars$hp))
#' yCov   <- cov(mtcars)[c('mpg', 'hp'), c('mpg', 'hp')]
#' 
#' n <- nrow(mtcars)
#' 
#' out <- lm_combine(coefs = coefs,
#'                   xMeans = xMeans,
#'                   xCov = xCov,
#'                   yMeans = yMeans,
#'                   yCov = yCov,
#'                   n = n,
#'                   weights = 1)
#' 
#' out
#' summary(lm((mtcars$mpg + mtcars$hp) ~ disp + cyl, data = mtcars))$coefficients

lm_combine <- function(coefs,
                       xMeans,
                       xCov,
                       yMeans,
                       yCov,
                       n,
                       weights = 1) {
  out <- list()
  
  nResp <- ncol(coefs)
  nCovar <- nrow(coefs)
  
  if (length(weights) == 1) {
    if (weights == 1) {
      weights <- rep(1, nResp)
    }
  }
  
  ## Calcualte model coefficients
  out$coefs <- apply(
    coefs,
    MARGIN = 1,
    FUN = function(rowi) {
      sum(rowi * weights)
    }
  )
  
  ## Calculate coefficient standard errors

  
  ## xMeans := column means of design matrix
  xMeans <- c(1, xMeans)
  
  ## varX := Covariance of design matrix
  varX <- cbind(rep(0, nCovar),
                rbind(rep(0, nCovar - 1),
                      xCov))
  ## XtX := X'X
  ## To do: re-write code to avoid loop and use matrix algebra
  XtX <- matrix(NA, nrow = nrow(varX), ncol = ncol(varX))
  for (i in 1:nrow(XtX)) {
    for (j in 1:nrow(XtX)) {
      XtX[i, j] <- varX[i, j] * (n - 1) + xMeans[i] * xMeans[j] * n
    }
  }
  ## YtY := y'y
  ## To do: re-write code to avoid loop and use matrix algebra
  YtY <- 0
  for (i in 1 : nResp){
    for (j in 1 : nResp){
      YtY <- YtY + yCov[i, j] * (n - 1) + yMeans[i] * yMeans[j] * n
    }
  }
  
  betas <- as.matrix(out$coefs, ncol = 1)
  varCoefs <- 
    as.double((YtY - t(betas) %*% XtX %*% betas) / (n - nCovar))  * solve(XtX)
  
  out$seCoefs <- sqrt(diag(varCoefs))
  
  return(out)
}
