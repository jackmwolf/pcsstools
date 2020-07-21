#' Approximate the relationship between a genotype and
#' the principal components of a set of phenotypes.
#'
#' @param slopes Vector of slopes of linear regressions of each
#'  phenotype as a function of the genotype
#' @param stdErrors Vector of standard errors for each slope coefficient
#' @param intercepts Vector of intercepts of linear regressions of each
#'  phenotype as a function of the genotype
#' @param varX Variance of X, the independent variable (genotype)
#' @param covPhenos variance-covariance matrix of all phenotypes
#' @param n Sample size
#' @param scaled Should summary statistics be converted to correspond to
#'  scaled response variables? (Or, are they already scaled?)
#' @param centered Center response variables? (Or, already centered?)
#' @param meanPhenos Required for \code{centered = T}. Vector of means
#'  of each phenotype.
#' @export

approx_pca <- function(slopes, stdErrors, intercepts, varX, covPhenos, n, scaled = T,
                      meanPhenos = NULL){
  
  centered <- TRUE ## Locked argument for now. Doesn't work else due to SVD != eigen PCA
  
  ## Calculate PC weights -----------------------------------------------------
  if (scaled == T){
    # Calculate correlation matrix
    corPhenos <- cov2cor(covPhenos)
    pcaWeights <- eigen(corPhenos)$vectors
  } else if (scaled == F){
    pcaWeights <- eigen(covPhenos)$vectors
  }
  
  ## Calculate scaled coefficients --------------------------------------------
  if (scaled == T){
    slopes <- slopes / sqrt(diag(covPhenos))
    stdErrors <- stdErrors / sqrt(diag(covPhenos))
  }
  if (centered == T) {
    intercepts <- (intercepts - meanPhenos)
    if (scaled == T){
      intercepts <- intercepts / sqrt(diag(covPhenos))
    }
  }
  
  
  
  
  ## For each set of weights (ith principal component),
  ## carry out approx_addition()
  out <- apply(pcaWeights, MARGIN = 2, FUN = function(x){
    approx_addition(weights = x,
                    slopes = slopes, intercepts = intercepts,
                    stdErrors = stdErrors, varX = varX,
                    covPhenos = covPhenos,
                    n = n)
  })
  
  ## also return pca weights
  out$pcaWeights <- pcaWeights
  
  return(out)
  
}
