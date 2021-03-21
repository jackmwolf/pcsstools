#' Approximate a linear model for a series of logical AND statements using PCSS
#'
#' \code{model_and} approximates the linear model for the conjunction
#'   of m phenotypes as a function of a set of predictors.
#'
#' @param formula an object of class \code{formula} whose dependent variable is
#'   a combination of variables and logical \code{&} operators. All model terms
#'   must be accounted for in \code{means} and \code{covs}.
#' @param n sample size.
#' @param means named vector of predictor and response means.
#' @param covs named matrix of the covariance of all model predictors and the
#'   responses.
#' @param predictors named list of objects of class \code{predictor}.
#' @param ... additional arguments
#'
#' @inherit pcsslm return
#' 
#' @references{
#' 
#'   \insertRef{wolf_using_2021}{pcsstools}
#'   
#' }
#' @examples
#' ex_data <- pcsstools_example[c("g1", "x1", "y4", "y5")]
#' head(ex_data)
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)
#' predictors <- list(
#'   g1 = new_predictor_snp(maf = mean(ex_data$g1) / 2),
#'   x1 = new_predictor_normal(mean = mean(ex_data$x1), sd = sd(ex_data$x1))
#' )
#'
#' model_and(
#'   y4 & y5 ~ g1 + x1,
#'   means = means, covs = covs, n = n, predictors = predictors
#' )
#' summary(lm(y4 & y5 ~ g1 + x1, data = ex_data))
#' @export
#'
model_and <- function(formula, n, means, covs, predictors, ...) {
  cl <- match.call()
  terms <- terms(formula)

  xterms <- extract_predictors(formula)
  yterms <- parse_and(extract_response(formula))
  
  check_terms_logical(xterms$predictors, yterms, means, covs, predictors, ...)
  

  # Re-arrange means, covs, and predictors to match given formula
  means0 <- means[c(xterms$predictors, yterms)]
  covs0 <- covs[c(xterms$predictors, yterms), c(xterms$predictors, yterms)]
  predictors0 <- predictors[xterms$predictors]
  add_intercept <- xterms$add_intercept

  re <- approx_and(
    means = means0, covs = covs0, n = n, predictors = predictors0,
    add_intercept = add_intercept, terms = terms
  )
  re$call <- cl
  class(re) <- "pcsslm"

  return(re)
}

#' Approximate a linear model for a series of logical OR statements using PCSS
#'
#' \code{model_or} approximates the linear model for the a disjunction
#'   of m phenotypes as a function of a set of predictors.
#'
#' @param formula an object of class \code{formula} whose dependent variable is
#'   a combination of variables and logical \code{|} operators. All model terms
#'   must be accounted for in \code{means} and \code{covs}.
#' @param n sample size.
#' @param means named vector of predictor and response means.
#' @param covs named matrix of the covariance of all model predictors and the
#'   responses.
#' @param predictors named list of objects of class \code{predictor}.
#' @param ... additional arguments
#' 
#' @inherit pcsslm return
#'
#' @references{
#' 
#'   \insertRef{wolf_using_2021}{pcsstools}
#'   
#' }
#' 
#' @examples
#' ex_data <- pcsstools_example[c("g1", "x1", "y4", "y5")]
#' head(ex_data)
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)
#' predictors <- list(
#'   g1 = new_predictor_snp(maf = mean(ex_data$g1) / 2),
#'   x1 = new_predictor_normal(mean = mean(ex_data$x1), sd = sd(ex_data$x1))
#' )
#'
#' model_or(
#'   y4 | y5 ~ g1 + x1,
#'   means = means, covs = covs, n = n, predictors = predictors
#' )
#' summary(lm(y4 | y5 ~ g1 + x1, data = ex_data))
#' @export
#'
model_or <- function(formula, n, means, covs, predictors, ...) {
  cl <- match.call()
  terms <- terms(formula)

  xterms <- extract_predictors(formula)
  yterms <- parse_or(extract_response(formula))
  
  check_terms_logical(xterms$predictors, yterms, means, covs, predictors, ...)

  # Re-arrange means, covs, and predictors to match given formula
  means0 <- means[c(xterms$predictors, yterms)]
  covs0 <- covs[c(xterms$predictors, yterms), c(xterms$predictors, yterms)]
  predictors0 <- predictors[xterms$predictors]
  add_intercept <- xterms$add_intercept

  re <- approx_or(
    means = means0, covs = covs0, n = n, predictors = predictors0,
    add_intercept = add_intercept,
    terms = terms, ...
  )
  re$call <- cl
  class(re) <- "pcsslm"

  return(re)
}


#' Approximate a linear model for a series of logical AND statements
#'
#' \code{approx_and} approximates the linear model for the a conjunction
#'   of m phenotypes as a function of a set of predictors.
#'
#' @param means vector of predictor and response means with the last \code{m}
#'   means being the means of \code{m} binary responses to combine in a
#'   logical and statement.
#' @param covs a matrix of the covariance of all model predictors and the
#'   responses with the order of rows/columns corresponding to the order of
#'   \code{means}.
#' @param n sample size.
#' @param predictors list of objects of class \code{predictor} corresponding
#'   to the order of the predictors in \code{means}.
#' @param add_intercept logical. Should the linear model add an intercept term?
#' @param response_assumption character. Either \code{"binary"} or
#'   \code{"continuous"}. If \code{"binary"}, specific calculations will be done
#'   to estimate product means and variances.
#' @param verbose should output be printed to console?
#' @param ... additional arguments
#'
#' @references{
#' 
#'   \insertRef{wolf_using_2021}{pcsstools}
#'   
#' }
approx_and <- function(means, covs, n, predictors, add_intercept = TRUE,
                       verbose = FALSE, response_assumption = "binary",
                       ...) {
  m <- length(means) - length(predictors)
  p <- length(means) - m

  # Generate responses' pmfs
  r_means <- means[(p + 1):(p + m)]
  responses <- lapply(r_means, new_predictor_binary)

  approx0 <- approx_mult_prod(
    means = means, covs = covs, n = n,
    response = response_assumption, responses = responses,
    predictors = predictors, verbose = verbose
  )


  model <- calculate_lm(
    means = approx0$means, covs = approx0$covs, add_intercept = add_intercept,
    n = n, ...
  )
  return(model)
}

#' Approximate a linear model for a series of logical OR statements
#'
#' \code{approx_or} approximates the linear model for a disjunction of m
#'   phenotypes as a function of a set of predictors.
#'
#' @param means vector of predictor and response means with the last m
#'   means being the means of m binary responses to combine in a
#'   logical OR statement.
#' @param covs a matrix of the covariance of all model predictors and the
#'   responses with the order of rows/columns corresponding to the order of
#'   \code{means}.
#' @param n sample size.
#' @param predictors list of objects of class \code{predictor} corresponding
#'   to the order of the predictors in \code{means}.
#' @param add_intercept logical. Should the linear model add an intercept term?
#' @param response_assumption character. Either \code{"binary"} or
#'   \code{"continuous"}. If \code{"binary"}, specific calculations will be done
#'   to estimate product means and variances.
#' @param verbose should output be printed to console?
#' @param ... additional arguments
#'
#' @references{
#' 
#'   \insertRef{wolf_using_2021}{pcsstools}
#'   
#' }
#' 
approx_or <- function(means, covs, n, predictors, add_intercept = TRUE,
                      verbose = FALSE, response_assumption = "binary", ...) {
  # Model "y1 or y2 or ..." via "not(not y1 and not y2 and ...)"
  m <- length(means) - length(predictors)
  p <- length(means) - m

  not_means <- c(means[1:p], 1 - means[(p + 1):(p + m)])
  phi <- c(rep(1, p), rep(-1, m))
  not_covs <- t(covs * phi) * phi

  # Generate responses' pmfs
  r_means <- not_means[(p + 1):(p + m)]
  responses <- lapply(r_means, new_predictor_binary)

  approx_not_and <- approx_mult_prod(
    means = not_means, covs = not_covs, n = n,
    response = response_assumption, responses = responses,
    predictors = predictors,
    verbose = verbose
  )

  tau <- c(rep(1, p), -1)
  out_covs <- t(approx_not_and$covs * tau) * tau
  out_means <- approx_not_and$means * tau + c(rep(0, p), 1)

  model <- calculate_lm(
    means = out_means, covs = out_covs, n = n, 
    add_intercept = add_intercept, ...
  )
  return(model)
}

