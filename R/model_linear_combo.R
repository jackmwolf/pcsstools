#' Model the principal component score of a set of phenotypes using PCSS
#'
#' \code{model_prcomp} calculates the linear model for the mth principal
#'   component score of a set of phenotypes as a function of a set of
#'   predictors.
#'
#' @param formula an object of class \code{formula} whose dependent variable is
#'   a series of variables joined by \code{+} operators. \code{model_prcomp}
#'   will treat a principal component score of those variables as the actual
#'   dependent variable. All model terms must be accounted for in \code{means}
#'   and \code{covs}.
#'
#' @param comp integer indicating which principal component score to analyze.
#'   Must be less than or equal to the total number of phenotypes.
#' @param n sample size.
#' @param means named vector of predictor and response means.
#' @param covs named matrix of the covariance of all model predictors and the
#'   responses.
#' @param center logical. Should the dependent variables be centered before
#'   principal components are calculated?
#' @param standardize logical. Should the dependent variables be standardized 
#'   before principal components are calculated? 
#' @param ... additional arguments
#' 
#' @inherit pcsslm return
#'
#' @references{
#'   
#'   \insertRef{wolf_computationally_2020}{pcsstools}
#'   
#' }
#'
#' @examples
#' ex_data <- pcsstools_example[c("g1", "x1", "x2", "y1", "y2", "y3")]
#' head(ex_data)
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)
#'
#' model_prcomp(
#'   y1 + y2 + y3 ~ g1 + x1 + x2,
#'   comp = 1, n = n, means = means, covs = covs
#' )
#' @export
model_prcomp <- function(formula, comp = 1, n, means, covs, 
                         center = FALSE, standardize = FALSE, ...) {
  cl <- match.call()
  terms <- terms(formula)

  xterms <- extract_predictors(formula)
  yterms <- parse_sum(extract_response(formula))
  
  check_terms_combo(xterms$predictors, yterms, means, covs)

  # Re-arrange means, covs, and predictors to match given formula
  means0 <- means[c(xterms$predictors, yterms)]
  covs0 <- covs[c(xterms$predictors, yterms), c(xterms$predictors, yterms)]
  add_intercept <- xterms$add_intercept
  
  # Adjust pcss if response is centered, standardized
  if (center) {
    means0[yterms] <- 0
  }
  if (standardize) {
    var0 <-  c(rep(1, length(xterms$predictors)), diag(covs0)[yterms])
    
    covs0 <- (1 / sqrt(var0)) * t((1 / sqrt(var0)) *  covs0)

  }

  # Calculate weights for PCA
  ysigma <- covs0[yterms, yterms]
  phi <- eigen(ysigma)$vectors[, comp]

  re <- calculate_lm_combo(
    means = means0, covs = covs0, n = n, phi = phi, 
    add_intercept = add_intercept, terms = terms, ...
  )
  re$call <- cl
  class(re) <- "pcsslm"

  return(re)
}


#' Model a linear combination of a set of phenotypes using PCSS
#'
#' \code{model_combo} calculates the linear model for a linear combination of
#'   phenotypes as a function of a set of predictors.
#'
#' @param formula an object of class \code{formula} whose dependent variable is
#'   a series of variables joined by \code{+} operators. \code{model_combo}
#'   will treat a principal component score of those variables as the actual
#'   dependent variable. All model terms must be accounted for in \code{means}
#'   and \code{covs}.
#'
#' @param phi named vector of linear weights for each variable in the
#'   dependent variable in \code{formula}.
#' @param n sample size.
#' @param means named vector of predictor and response means.
#' @param covs named matrix of the covariance of all model predictors and the
#'   responses.
#' @param ... additional arguments
#'
#' @references{
#' 
#'   \insertRef{wolf_computationally_2020}{pcsstools}
#'   
#'   \insertRef{gasdaska_leveraging_2019}{pcsstools}
#' 
#' }
#'
#' @inherit pcsslm return
#' @examples
#' ex_data <- pcsstools_example[c("g1", "x1", "x2", "x3", "y1", "y2", "y3")]
#' head(ex_data)
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)
#' phi <- c("y1" = 1, "y2" = -1, "y3" = 0.5)
#'
#' model_combo(
#'   y1 + y2 + y3 ~ g1 + x1 + x2 + x3, 
#'   phi = phi, n = n, means = means, covs = covs
#' )
#'
#' summary(lm(y1 - y2 + 0.5 * y3 ~ g1 + x1 + x2 + x3, data = ex_data))
#' @export
model_combo <- function(formula, phi, n, means, covs, ...) {
  
  cl <- match.call()
  terms <- terms(formula)

  xterms <- extract_predictors(formula)
  yterms <- parse_sum(extract_response(formula))

  check_terms_combo(xterms$predictors, yterms, means, covs)
  
  # Re-arrange means, covs, and predictors to match given formula
  means0 <- means[c(xterms$predictors, yterms)]
  covs0 <- covs[c(xterms$predictors, yterms), c(xterms$predictors, yterms)]
  add_intercept <- xterms$add_intercept

  phi0 <- phi[yterms]

  re <- calculate_lm_combo(
    means = means0, covs = covs0, n = n, phi = phi, 
    add_intercept = add_intercept, terms = terms, ...
  )
  re$call <- cl
  class(re) <- "pcsslm"

  return(re)
}


#' Model an individual phenotype using PCSS
#'
#' \code{model_singular} calculates the linear model for a singular
#'   phenotype as a function of a set of predictors.
#'
#' @param formula an object of class \code{formula} whose dependent variable is
#'   only variable. All model terms must be accounted for in \code{means}
#'   and \code{covs}.
#' @param n sample size.
#' @param means named vector of predictor and response means.
#' @param covs named matrix of the covariance of all model predictors and the
#'   responses.
#' @param ... additional arguments
#'
#' @references{
#' 
#'   \insertRef{wolf_computationally_2020}{pcsstools}
#'   
#' }
#'
#' @inherit pcsslm return
#'
#' @export
#' @examples 
#' ex_data <- pcsstools_example[c("g1", "x1", "y1")]
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)
#'
#' model_singular(
#'   y1 ~ g1 + x1,
#'   n = n, means = means, covs = covs
#' )
#' summary(lm(y1 ~ g1 + x1, data = ex_data))
model_singular <- function(formula, n, means, covs, ...) {
  cl <- match.call()
  terms <- terms(formula)
  
  xterms <- extract_predictors(formula)
  yterms <- parse_sum(extract_response(formula))
  
  check_terms_combo(xterms$predictors, yterms, means, covs)
  
  # Re-arrange means, covs, and predictors to match given formula
  means0 <- means[c(xterms$predictors, yterms)]
  covs0 <- covs[c(xterms$predictors, yterms), c(xterms$predictors, yterms)]
  add_intercept <- xterms$add_intercept
  
  re <- calculate_lm(
    means = means0, covs = covs0, n = n, add_intercept = add_intercept, 
    terms = terms, ...
  )
  re$call <- cl
  class(re) <- "pcsslm"
  
  return(re)
}




#' Calculate a linear model for a linear combination of responses
#'
#' \code{calculate_lm_combo} describes the linear model for a linear combination
#'   of responses as a function of a set of predictors.
#'
#' @param means a vector of means of all model predictors and the response with
#'   the last \code{m} elements the response means (with order corresponding to
#'   the order of weights in \code{phi}).
#' @param covs a matrix of the covariance of all model predictors and the
#'   responses with the order of rows/columns corresponding to the order of
#'   \code{means}.
#' @param n sample size.
#' @param m number of responses to combine. Defaults to \code{length(weighs)}.
#' @param phi vector of linear combination weights with one entry per response
#'   variable.
#' @param add_intercept logical. If \code{TRUE} adds an intercept to the model.
#' @param ... additional arguments
#'
#' @references{
#' 
#'   \insertRef{wolf_computationally_2020}{pcsstools}
#'   
#'   \insertRef{gasdaska_leveraging_2019}{pcsstools}
#' 
#' }
#' @inherit pcsslm return
calculate_lm_combo <- function(means, covs, n, phi, m = length(phi), 
                               add_intercept, ...) {
  p <- length(means) - m

  # Covariances with linear combo and variance/mean of the linear combo
  new_covs <- covs[1:p, (p + 1):(p + m)] %*% phi
  new_var <- drop(t(phi) %*% covs[(p + 1):(p + m), (p + 1):(p + m)] %*% phi)
  new_mean <- sum(phi * means[(p + 1):(p + m)])

  means0 <- c(means[1:p], new_mean)
  covs0 <- rbind(cbind(covs[1:p, 1:p], new_covs), c(t(new_covs), new_var))
  
  colnames(covs0) <- c(names(means)[1:p], NA)
  rownames(covs0) <- c(names(means)[1:p], NA)

  calculate_lm(means = means0, covs = covs0, n = n, 
               add_intercept = add_intercept, ...)
}
