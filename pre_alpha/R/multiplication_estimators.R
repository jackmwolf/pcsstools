#' Approximate coefficients for a product of phenotypes
#'
#' @param slopes A vector of slope coeffients for each phenotype as a function
#'  of a common genotype
#' @param intercepts A vector of intercepts for each phenotype as a function of
#'  a common genotype
#' @param inheritance Phenotype inheritance method. Either 'domiant' or
#'  'additive'
#' @param maf Minor allele frequency
#' @param varCovar Variance-covariance matrix of all phenotypes
#' @param yBar Vector of means of phenotypes
#' @param n Sample size
#'
#' @seealso \code{\link{getMultCoefs}} for coeffient estimation and
#'  \code{\link{getMultSE}} for coefficient standard error estimation.
#'
#' @export


approxPhenoProd <- function(slopes, intercepts, maf, inheritance = 'dominant',
                            varCovar, yBar, n){
  coefs <- getMultCoefs(slopes = slopes, intercepts = intercepts,
                        maf = maf, inheritance = inheritance)
  se <- getMultSE(slopes = slopes, intercepts = intercepts,
                  maf = maf, varCovar = varCovar, yBar = yBar,
                  n = n, predBeta = coefs$predSlope,
                  inheritance = inheritance)
  
  out <- list()
  out$coefs <- coefs
  out$se <- se
  
  return(out)
}



#' Approximate coefficients for a product of phenotypes
#'
#' @param slopes A vector of slope coeffients for each phenotype as a function
#'  of a common genotype
#' @param intercepts A vector of intercepts for each phenotype as a function of
#'  a common genotype
#' @param maf Minor allele frequency. Required for inheritance = 'additive'
#' @param inheritance Phenotype inheritance method. Either 'domiant' or
#'  'additive'

getMultCoefs <- function(slopes, intercepts, maf = NULL,
                                    inheritance = 'dominant'){
  
  coefs <- cbind(slopes, intercepts)
  
  predIntercept <- prod(coefs[, 'intercepts'])

  
  ## Dominant inheritance =====================================================
  if (inheritance == 'dominant'){
    predSlope <- prod(rowSums(coefs)) - prod(coefs[, 'intercepts'])
    
  }
    
  ## Additive inheritance =====================================================
  if (inheritance == 'additive'){
    ## Use maf to calculate expected proportions of x = 0, 1, and 2
    
    p0 <- (1 - maf) ^ 2
    p1 <- 2 * maf * (1 - maf)
    p2 <- maf ^ 2
    
    ## Estimated sample means at each value of x
    mu0 <- prod(coefs[, 'intercepts'])
    mu1 <- prod(rowSums(coefs))
    mu2 <- prod(rowSums(coefs %*% matrix(c(1, 0, 0, 2), ncol = 2)))
    
    out$predSlope <- ((p0 + p1 + p2) * (p1 * mu1 + 2 * p2 * mu2) -
                      (p1 + 2 * p2) * (p0 * mu0 + p1 * mu1 + p2 * mu2)) /
                     ((p0 + p1 + p2) * (p1 + 4 * p2) - (p1 + 2 * p2) ^ 2)
    
  }
  
  ## Use maf to calculate expected proportions of x = 0, 1, and 2
  
  
  out <- list()
  
  out$predSlope <- predSlope
  out$predIntercept <- predIntercept
  
  return(out)
  
}

#' Approximate slope standard error for a product of phenotypes
#'
#' @param slopes A vector of slope coeffients for each phenotype as a function
#'  of a common genotype
#' @param intercepts A vector of intercepts for each phenotype as a function of
#'  a common genotype
#' @param maf Minor allele frequency
#' @param varCovar Variance-covariance matrix of all phenotypes
#' @param yBar Vector of means of phenotypes
#' @param predBeta Predicted slope from getCoefs()
#' @param n Sample size
#' @param inheritance Phenotype inheritance method. Either 'domiant' or
#'  'additive'

getMultSE <- function(slopes, intercepts, maf, varCovar, yBar, n, predBeta,
                                 inheritance){
  
  # Vector to hold slopes for y_1 * y_2 * ... * y_n-1 ~ y_n
  recBeta <- vector(length = length(slopes))
  recBeta[1] <- 0
  
  if (inheritance == 'dominant'){
    xBar <- maf
  } else if (inheritance == 'additive'){
    xBar <- maf * 2
  }
  
  ## Dominant inheritance =====================================================
  for (i in 2:length(slopes)){
    pairWiseBetas <- varCovar[i, 1 : (i - 1)] / varCovar[i, i]
    
    recBeta[i] <- getMultCoefs(slopes = pairWiseBetas,
                                       intercepts = yBar[1 : i - 1] - pairWiseBetas * xBar,
                                       maf = maf,
                                       inheritance = 'dominant' # dominant inherntance == binaryexplanatory
                                       )$predSlope
  }
  
  
  # Estimate the sum of y_1,i * ... * y_n,i -----------------------------------
  # product_{k = i+1}^j yBar_k
  yBarProdTail <- vector(length = length(slopes))
  yBarProdTail[length(slopes)] <- 1
  for (i in 1 : (length(slopes) - 1)){
    yBarProdTail[i] <- prod(yBar[i + 1 : length(yBar)], na.rm = T)
  }
  
  # est sum of y_1i * y_2i * ... y_ni
  sumYprod <- prod(yBar) * n + (n - 1) * sum(recBeta * diag(varCovar) * yBarProdTail)
  
  # variance of x
  if (inheritance == 'dominant'){
    varX <- maf * (1 - maf)
  } else if (inheritance == 'additive'){
    varX <- maf * (1 - maf) * 2
  }
  
  SE <- sqrt((sumYprod - (sumYprod ^ 2) / n - predBeta ^ 2 * (n - 1) * varX) /
             ((n - 2) * varX * (n - 1)))
 
 
 SE <- SE
 #out$sumYprod <- sumYprod
 
 return(SE)
}




