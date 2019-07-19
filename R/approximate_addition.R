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
#' @param covPhenos (Optional) Variance-covariance matrix of all phenotypes
#'
#' @export
#'

approx_addition <- function(slopes, stdErrors, intercepts,
                                 weights = 1, varX, n, covPhenos = NULL){

  out <- list()

  predSlope <- sum(slopes * weights)
  predIntercept <- sum(intercepts * weights)
  
  out$predSlope <- predSlope
  out$predIntercept <- predIntercept
  
  ## Estimate covariance of phenotypes if not provided --------------
  ## Or ... make separate function for before hand. Looks like it
  ## assumes other information is known that would be annoying to
  ## also use in here...
  if (is.null(covPhenos)){
    ## TO DO ========================================================
  }
  
  
  ## Current bug if default weights ... problem with == if weights is vector of
  ## length > 1
  #if (weights == 1){
  #  weights <- rep(1, length(slopes))
  #}

  ## calculate \sum_{q=1}^{m-1} \sum_{r = q+1}^m c_q c_r cov(y_q, y_r)
  ## (See eq. 10 in Gadaska 2019)
  m <- length(slopes)
  covSum <- 0
  for (q in 1 : (m - 1)){
    #print(paste('q = ', q))
    for (r in (q + 1) : m){
    #print(paste('r =', r))
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
  
  out$se <- se
  
  return(out)
}
