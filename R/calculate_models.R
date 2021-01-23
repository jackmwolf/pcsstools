#' Calculate a linear model
#'
#' \code{calculate_lm} describes the linear model of the last listed variable
#' in \code{means} and \code{covs} as a function of all other variables in
#' \code{means} and \code{covs}.
#'
#' @param means a vector of means of all model predictors and the response with
#'   the last element the response mean.
#' @param covs a matrix of the covariance of all model predictors and the
#'   response with the order of rows/columns corresponding to the order of
#'   \code{means}.
#' @param n sample size
#' @param add_intercept logical. If \code{TRUE} adds an intercept to the model.
#' 
#' @importFrom stats pt
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
#' ex_data <- cont_data[c("g", "x", "y1")]
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)
#' calculate_lm(means = means, covs = covs, n = n, add_intercept = TRUE)
#' # Compare results to...
#' mod <- lm(y1 ~ 1 + g + x, data = ex_data)
#' summary(mod)$coef
#' summary(mod)$sigma^2
#'
#' @export
calculate_lm <- function(means, covs, n, add_intercept = FALSE) {
  p <- ncol(covs) - 1
  if (add_intercept) {
    p <- p + 1
    covs <- rbind(0, cbind(0, covs))
    means <- c(1, means)
    if (!is.null(dimnames(covs)) & !is.null(dimnames(covs))) {
      names(means)[1] <- "(Intercept)"
      colnames(covs)[1] <- "(Intercept)"
      rownames(covs)[1] <- "(Intercept)" 
    }
  }
  covX <- covs[-(p + 1), -(p + 1)]
  meanX <- means[-(p + 1)]
  XtX <- (n - 1) * covX + n * meanX %*% t(meanX)
  Xty <- (n - 1) * covs[p + 1, -(p + 1)] + n * meanX * means[p + 1]
  yty <- (n - 1) * covs[p + 1, p + 1] + n * means[p + 1]^2

  beta <- drop(solve(XtX) %*% Xty)
  
  sigma2 <- drop((yty - t(beta) %*% Xty) / (n - p))

  var_beta <- diag(sigma2 * solve(XtX))
  sd_beta <- sqrt(var_beta)

  t_stat <- beta / sd_beta
  p_val <- 2 * pt(abs(t_stat), df = n - p, lower.tail = F)

  return(list(beta = beta, sd_beta = sd_beta,
              t_stat = t_stat, p_val = p_val,
              sigma2 = sigma2))
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
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)
#' phi <- c(1, 0.5, 2)
#' calculate_lm_combo(means = means, covs = covs, n = n, phi = phi, m = 3, add_intercept = TRUE)
#' # Compare results to...
#' mod <- lm(y1 + 0.5 * y2 + 2 * y3 ~ 1 + g + x, data = ex_data)
#' summary(mod)$coef
#'
#' @export
calculate_lm_combo <- function(means, covs, n, phi, m = length(phi), add_intercept) {
  p <- length(means) - m

  # Covariances with linear combo and variance/mean of the linear combo
  new_covs <- covs[1:p, (p + 1):(p + m)] %*% phi
  new_var  <- drop(t(phi) %*% covs[(p + 1):(p + m), (p + 1):(p + m)] %*% phi)
  new_mean <- sum(phi * means[(p + 1):(p + m)])

  means0 <- c(means[1:p], new_mean)
  covs0 <- rbind(cbind(covs[1:p, 1:p], new_covs), c(t(new_covs), new_var))

  calculate_lm(means = means0, covs = covs0, n = n, add_intercept = add_intercept)
}

#' Extract independent variables from a formula
#' @param formula an object of class \code{formula}.
#' @importFrom stats terms
extract_predictors <- function(formula = formula()) {
  terms0 <- terms(formula)
  
  model_terms <- attributes(terms0)$term.labels
  add_intercept <- attributes(terms0)$intercept
  
  re <- list(
    predictors = model_terms,
    add_intercept = add_intercept
  )
  re
}
#' Extract dependent variables from a formula as a string
#' @param formula an object of class \code{formula}.
#' @importFrom stats terms
extract_response <- function(formula = formula()) {
  terms0 <- terms(formula)
  response <- as.character(attributes(terms0)$variables)[[2]]
  return(response)
}


parse_response <- function(response = character(), split = character()) {
  terms0 <- strsplit(response, split = split)[[1]]
  # Trim white space at start and end
  terms0 <- sapply(terms0, trimws, simplify = T)
  return(terms0)
}

parse_product <- function(response = character()) {
  terms0 <- unname(parse_response(response = response, split = "\\*"))
  return(terms0)
}

parse_or <- function(response = character()) {
  terms0 <- unname(parse_response(response = response, split = "\\|"))
  return(terms0)
}

parse_and <- function(response = character()) {
  terms0 <- unname(parse_response(response = response, split = "\\&"))
  return(terms0)
}


