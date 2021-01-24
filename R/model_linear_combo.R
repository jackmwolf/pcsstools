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
#' @examples 
#' ex_data <- cont_data[c("g", "x", "y1", "y2", "y3")]
#' head(ex_data)
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)   
#' 
#' model_prcomp(
#'   y1 + y2 + y3 ~ g + x, comp = 1, n = n, means = means, covs = covs
#' )
#' 
#' 
#' @export
model_prcomp <- function(formula, comp = 1, n, means, covs) {
  all_vars <- names(means)
  
  xterms <- extract_predictors(formula, all_vars)
  yterms <- parse_sum(extract_response(formula), all_vars)
  
  # Re-arrange means, covs, and predictors to match given formula
  means0 <- means[c(xterms$predictors, yterms)]
  covs0  <- covs[c(xterms$predictors, yterms), c(xterms$predictors, yterms)]
  add_intercept <- xterms$add_intercept
  
  # Calculate weights for PCA
  ysigma <- covs0[yterms, yterms]
  phi <- eigen(ysigma)$vectors[, comp]
  
  model <- calculate_lm_combo(
    means = means0, covs = covs0, n = n, phi = phi, add_intercept = add_intercept
  )
  
  return(model)
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
#' @examples 
#' ex_data <- cont_data[c("g", "x", "y1", "y2", "y3")]
#' head(ex_data)
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)   
#' phi <- c("y1" = 1, "y2" = -1, "y3" = 0.5)
#' 
#' model_combo(
#'   y1 + y2 + y3 ~ g + x, phi = phi, n = n, means = means, covs = covs
#' )
#' 
#'
#' coef(summary(lm(y1 - y2 + 0.5 * y3 ~ g + x, data = ex_data)))
#' 
#' @export
model_combo <- function(formula, phi, n, means, covs) {
  all_vars <- names(means)
  
  xterms <- extract_predictors(formula, all_vars)
  yterms <- parse_sum(extract_response(formula), all_vars)
  
  # Re-arrange means, covs, and predictors to match given formula
  means0 <- means[c(xterms$predictors, yterms)]
  covs0  <- covs[c(xterms$predictors, yterms), c(xterms$predictors, yterms)]
  add_intercept <- xterms$add_intercept
  
  phi0 <- phi[yterms]
  
  model <- calculate_lm_combo(
    means = means0, covs = covs0, n = n, phi = phi, add_intercept = add_intercept
  )
  
  return(model)
}