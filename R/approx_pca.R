#' Approximate the relationship between a genotype and
#' the principal components of a set of phenotypes.
#'
#' @param slopes Vector of slopes of linear regressions of each
#'  phenotype as a function of the genotype
#' @param stdErrors Vector of standard errors for each slope coefficient
#' @param intercepts Vector of intercepts of linear regressions of each
#'  phenotype as a function of the genotype
#' @param varX Variance of X, the independent variable (genotype)
#' @param covPhenos (Optional) Variance-covariance matrix of all phenotypes
#' @param n Sample size
#'
#' @export

approx_pca <- function(slopes, stdErrors, intercepts, varX, covPhenos, n){
  
  ## Question: prcomp() returns -1 * eigen(), assume that
  ## it gives us what we wants so I multiply by -1 for now...
  pcaWeights <- -1 * eigen(covPhenos)$vectors
  
  ## For each set of weights (ith principal component),
  ## carry out approx_addition()
  out <- apply(pcaWeights, MARGIN = 1, FUN = function(x){
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
