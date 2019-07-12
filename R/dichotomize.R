#' Approximate coefficients for a dichotomous version of a continuous variable
#'
#' '\code{dichotmize} approximates the relationship between a genotype and
#' a dichotomous version of a phenotype given information about the relation
#' between the genotype and a linear version of the phenotype. E.g. given
#' information about the relationship between X and blood pressure, approximate
#' the model between X and wether a patient has high blood pressure
#' (BP > 140mm Hg)
#'
#' @param slope Slope coefficient for the original continuous response as
#'   a function of the genoptype
#' @param intercept Intercept coefficient for the original continuous
#'   response as a function of the genoptype
#' @param varY Variance of the original continuous response
#' @param varX Variance of the genotype
#' @param maf Minor allele frequency
#' @param cutoff Cutoff point for dichotomization
#' @param n Sample size

approxDichotomous <- function(slope, intercept, varY, varX, maf, cutoff,
                              n, inheritance){
  
  residVar <- (1 - (slope * sqrt(varX) / sqrt(varY)) ^ 2) * varY
  
  out <- list()
  
  ## Calculate slope and intercept ---------------------------------------------
  predIntercept <- pnorm(q = cutoff, mean = intercept, sd = sqrt(residVar),
                         lower.tail = F)
  predSlope <- pnorm(q = cutoff, mean = intercept + slope,
                     sd = sqrt(residVar), lower.tail = F) - predIntercept
  
  out$predIntercept <- predIntercept
  out$predSlope <- predSlope
  
  
  ## Calculate standard error -------------------------------------------------
  ## Currently only possible for dominant inheritance
  if (inheritance == 'dominant'){
  ## Expected number of Y|x=0
  nY1X0 <- n * predIntercept * 1 - maf
  
  ## Expected number of 1-Y|x=0
  nY0X0 <- n * (1 - predIntercept) * 1 - maf
  
  ## Expected number of Y|X=1
  nY1X1 <- n * (predSlope + predIntercept) * maf
  
  ## Expected number of 1-Y|X=1
  nY0X1 <- n * (1 - predSlope + predIntercept) * maf
  
  predSE <- sqrt((nY1X0 * (1 - predIntercept) ^ 2 + nY1X1 * (1 - predSlope - predIntercept) ^ 2 +
                  nY0X0 * (0 - predIntercept) ^ 2 + nY0X1 * (0 - predSlope - predIntercept) ^ 2) /
                 ((n - 2) * varX * (n - 1)))
  
  out$predSE <- predSE
  }
  
  return(out)
  
}
