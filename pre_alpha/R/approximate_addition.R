#' Approximate the sum of two variables based on summary statistics.
#'
#' @param slopes Vector of slopes of linear regressions of each
#'  phenotype as a function of the genotype
#' @param stdErrors Vector of standard errors for each slope coefficient
#' @param intercepts Vector of intercepts of linear regressions of each
#'  phenotype as a function of the genotype
#' @param weights vector of weights for each phenotype (defaults to 1)
#' @param varX
#' @param n Sample size
#' @param covPhenos Variance-covariance matrix of all phenotypes
#'
#' @export
#'

approx_addition <- function(slopes, stdErrors, intercepts,
                                 weights = 1, varX, n, covPhenos){

  out <- list()

  predSlope <- sum(slopes * weights)
  predIntercept <- sum(intercepts * weights)
  
  out$predIntercept <- predIntercept
  out$predSlope <- predSlope
  
  
  ## Current bug if default weights ... problem with == if weights is vector of
  ## length > 1
  if (length(weights) == 1 & weights == 1){
    weights <- rep(1, length(slopes))
  }

  m <- length(slopes)
  if (m == 2){
    ## If only 2 phenotypes, use eqn 7 from Gadaska 2019
    se <- sqrt(sum(weights ^ 2 * stdErrors ^ 2) + 2 * prod(weights) / (n - 2) *
               (as.double(covPhenos) / varX - prod(slopes)))
  } else {
    ## If > 2 phenotypes ... 

    ## calculate \sum_{q=1}^{m-1} \sum_{r = q+1}^m c_q c_r cov(y_q, y_r)
    ## (See eq. 10 in Gadaska 2019)
    covSum <- 0
    for (q in 1 : (m - 1)){
      print(paste('q = ', q))
      for (r in (q + 1) : m){
      print(paste('r =', r))
        covSum <- covSum +
          (weights[q] * weights[r] * covPhenos[q, r])
      }
    }
    
    ## calculate \sum_{q=1}^{m-1} \sum_{r = q+1}^m c_q c_r beta_q beta_r
    betaPairSum <- 0
    for (q in 1 : (m - 1)){
      for (r in (q + 1) : m){
        betaPairSum <- betaPairSum +
          (weights[q] * weights[r] * slopes[q] * slopes[r])
      }
    }
    
    se <- sqrt(sum(weights ^ 2 * stdErrors ^ 2) +
      2 / (n - 2) * ((covSum / varX) - betaPairSum))
  }

  out$se <- se
  
  return(out)
}

#' Approx the sum of two variables based on summary statistics with covariates.
#'
#' @param coefs A matrix of coefficients. Each column represents intercept,
#'  and slope coefficients for a linear regression of a phenotype as a function
#'  of a genotype and other covariates.
#' @param covariateMeans Vector of means of covariates
#' @param covariateVars Variance-covariance matrix of all covariates
#' @param weights Vector of weights for each phenotype (defaults to 1)
#' @param respVars Variance-covariance matrix of all phenotypes
#' @param respMeans Vector of means of all phenotypes
#' @param n Sample size
#'
#'

approx_addition.covars <-
  function(coefs, covariateMeans, covariateVars, weights = 1,
           respVars, respMeans, n){
  
  tryCatch({
  out <- list()

  out$coefs <- apply(coefs, MARGIN = 1, FUN = function(rowi){
                     if(length(rowi) != length(weights)){
                     print('length rowi != weights')
                     print(rowi)
                     print(weights)
                     }
    sum(rowi * weights)
  })
  
  betas <- as.matrix(out$coefs, ncol = 1)


  nResp <- ncol(coefs)
  nCovars <- nrow(coefs)
  
  if (length(covariateMeans) != nCovars){
    covariateMeans <- c(1, covariateMeans)
  }
  
  
  
  varC <- cbind(rep(0, nCovars), rbind(rep(0, nCovars - 1), covariateVars))
  
  XtX <- matrix(NA, nrow = nrow(varC), ncol = ncol(varC))
  for (i in 1:nrow(XtX)){
    for (j in 1:nrow(XtX)){
      XtX[i, j] <- varC[i, j] * (n - 1) + covariateMeans[i] * covariateMeans[j] * n
    }
  }
  
  
  YtY <- 0
  for (i in 1:nResp){
    for (j in 1:nResp){
      YtY <- YtY + weights[i] * weights[j] *
      (respVars[i, j] * (n - 1) + respMeans[i] * respMeans[j] * n)
    }
  }
  
  varCoefs <- as.double((YtY - t(betas) %*% XtX %*% betas) / (n - nCovars))  * solve(XtX)
  
  out$seBeta <- sqrt(diag(varCoefs))
  }, error = function(e){out <- list(NA, NA)})
  return(out)
  
}


#' Adjust for covariates from simple linear regressions.
#'
#' @param coefs A matrix of coefficients. Each column represents intercept,
#'  and slope coefficient for a linear regression of a common phenotype as a 
#'  function of a single covariates.
#' @param covariateVars Variance-covariance matrix of all covariates
#' @param respVar Variance of the phenotype
#' @param respMean Mean of the phenotype
#' @param n Sample size
#'

covar_adjust <- function(coefs, covariateVars, respMean, respVar, n){
  out <- list()
  tryCatch({

  nResp <- 1
  nCovars <- ncol(coefs)
  
  # column means of design matrix
  xMeans <- c(1, apply(coefs, MARGIN = 2, FUN = 
                       function(x){(respMean - x[1]) / x[2]}
                       )) 

  # variance-covariance of design matrix
  varC <- cbind(rep(0, nCovars + 1), rbind(rep(0, nCovars), covariateVars))
  
  XtX <- matrix(NA, nrow = nrow(varC), ncol = ncol(varC))
  for (i in 1:nrow(XtX)){
    for (j in 1:nrow(XtX)){
      XtX[i, j] <- varC[i, j] * (n - 1) + xMeans[i] * xMeans[j] * n
    }
  }

  XtY <- matrix(NA, nrow = nCovars + 1,  ncol = 1) 
  XtY[1, 1] <- respMean * n
  for (i in 2 : (nCovars + 1)){
      XtY[i, 1] <- coefs[2, i - 1] * varC[i, i] * (n - 1) + xMeans[i] * respMean * n 
  }

  betas <- solve(XtX) %*% XtY
  
  out$coefs <- betas

  YtY <- respVar * (n - 1) + respMean ^ 2 * n

  #varCoefs <- as.double(YtY - t(betas) %*% XtY) * solve(XtX)
  varCoefs <- as.double((YtY - t(betas) %*% XtX %*% betas) / (n - (nCovars + 1)))  * solve(XtX)
  out$seBeta <- sqrt(diag(varCoefs))
  # out$xMeans <- xMeans
  # out$varC <- varC
  # out$XtX <- XtX
  # out$XtY <- XtY
  }, error = function(e){out <- list(NA, NA)})
  
  return(out)



}

## revised function that takes xBar as an input instead as a workaround
## to issues with incorrect xBar calculations from summary stats due
## to missing data...
## Note: xMeans := column means of design matrix.
##      -> first element of xMeans should be '1' (mean of c(1, 1, ... ,1))
covar_adjust2 <- function(coefs, covariateVars, respMean, respVar, n, xMeans){
  out <- list()
  tryCatch({

  nResp <- 1
  nCovars <- ncol(coefs)
  
  # variance-covariance of design matrix
  varC <- cbind(rep(0, nCovars + 1), rbind(rep(0, nCovars), covariateVars))
  
  XtX <- matrix(NA, nrow = nrow(varC), ncol = ncol(varC))
  for (i in 1:nrow(XtX)){
    for (j in 1:nrow(XtX)){
      XtX[i, j] <- varC[i, j] * (n - 1) + xMeans[i] * xMeans[j] * n
    }
  }

  XtY <- matrix(NA, nrow = nCovars + 1,  ncol = 1) 
  XtY[1, 1] <- respMean * n
  for (i in 2 : (nCovars + 1)){
      XtY[i, 1] <- coefs[2, i - 1] * varC[i, i] * (n - 1) + xMeans[i] * respMean * n 
  }

  betas <- solve(XtX) %*% XtY
  
  out$coefs <- betas

  YtY <- respVar * (n - 1) + respMean ^ 2 * n

  varCoefs <- as.double((YtY - t(betas) %*% XtX %*% betas) / (n - (nCovars + 1)))  * solve(XtX)
  out$seBeta <- sqrt(diag(varCoefs))
  # out$xMeans <- xMeans
  # out$varC <- varC
  # out$XtX <- XtX
  # out$XtY <- XtY
  }, error = function(e){out <- list(NA, NA)})
  
  return(out)



}
