#' Approximate the covariance of a set of predictors and a product of responses
#'
#' \code{approx_mult_prod} recursively estimates the covariances and means of a
#'   set of responses
#'
#' @param means a vector of predictor and response means with all response means
#'   at the end of the vector.
#' @param covs covariance matrix of all predictors and responses with column
#'   and row order corresponding to the order of \code{means}.
#' @param n sample size (an integer).
#' @param response a string. Currently supports \code{"binary"} or
#'   \code{"continuous"}.
#' @param predictors,responses lists of objects of class \code{predictor} where
#'   each entry corresponds to one predictor/response variable.
#' @param verbose logical.
#' 
#' @importFrom stats median integrate cov2cor
#'
#' @examples
#' # 3 Binary Phenotypes
#' ex_data <- bin_data[c("g", "x", "y1", "y2", "y3")]
#' head(ex_data)
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)
#' predictors <- list(
#'   new_predictor_snp(maf = mean(ex_data$g) / 2),
#'   new_predictor_normal(mean = mean(ex_data$x), sd = sd(ex_data$x))
#' )
#' responses <- lapply(means[3:length(means)], new_predictor_binary)
#'
#' approx_mult_prod(means, covs, n, response = "binary",
#'   predictors = predictors, responses = responses, verbose = TRUE)
#' # Compare to
#' w3 <- with(ex_data, y1 * y2 * y3)
#' colMeans(cbind(ex_data[1:2], w3))
#' cov(cbind(ex_data[1:2], w3))
#'
#' @export
approx_mult_prod <- function(means, covs, n, response, predictors, responses, verbose = FALSE) {
  # Number of responses
  m <- length(means) - length(predictors)
  p <- length(means) - m

  if (m > 3) {
    warning("Estimation with products of >3 responses has not yet been tested for accuracy.")
  }

  perms <- lapply(make_permutations(m), rev)
  approx <-
    lapply(perms, function(order0) {

      # Re-arrange order of means and covs to match order0
      means0 <- means[c(1:p, p + order0)]
      covs0  <- covs[c(1:p, p + order0), c(1:p, p + order0)]
      predictors0 <- predictors[order0]

      # Print estimation order
      if (verbose & m > 2) {
        resp_names <- names(means0)[(p + 1):(p + m)]
        if (any(is.null(resp_names)) | "" %in% resp_names) {
          order_chr <- paste(rev(order0), collapse = " * ")
          cat(paste("Approximating with response columns order: "), order_chr, "\n")
        } else {
          y_names <- paste(rev(resp_names), collapse = " * ")
          cat(paste("Approximating with responses ordered as: ", y_names, "\n"))
        }
      }

      # Estimate covariances and means recursively according to order0
      approx0 <-
        approx_prod_recursive(means = means0, covs = covs0, n = n,
                              response = response, predictors = predictors,
                              responses = responses[order0])
    })

  # Take pairwise medians of estimates
  pred_mean <- median(
    sapply(approx, function(.) .$means[p + 1])
  )
  pred_covs <- apply(
    sapply(approx, function(.) .$covs[, p + 1]),
    MARGIN = 1, FUN = median
  )


  means_out <- c(means[1:p], pred_mean)
  covs_out  <- rbind(cbind(covs[1:p, 1:p], pred_covs[1:p]), pred_covs)

  # Prepare names
  prod_name <- paste(names(means)[(p + 1) : (p + m)], collapse = "")
  names(means_out)[p + 1]  <- prod_name
  colnames(covs_out)[p + 1] <- prod_name
  rownames(covs_out)[p + 1] <- prod_name

  return(list(means = means_out, covs = covs_out))
}


#' Recursively approximate covariances with a product of responses
#'
approx_prod_recursive <- function(means, covs, n, response, predictors, responses) {
  # Number of responses and predictors
  m <- length(means) - length(predictors)
  p <- length(means) - m

  if (m == 2) {
    cov_out <- approx_prod_stats(means = means, covs = covs, n = n,
                                 response = response, predictors = predictors)
    return(cov_out)
  } else if (m > 2) {
    # mean of y_m, covariance with each predictor, and variance
    mean_m <- means[p + 1]
    covs_m  <- covs[1:(p + 1), p + 1]

    # mean of w_(m-1), covariance with each predictor, and variance
    ids2 <- -(p  + 1)
    approx_w <- approx_prod_recursive(
      means = means[ids2],
      covs = covs[ids2, ids2],
      n = n, response = response, predictors = predictors,
      responses = responses[-1]
    )



    # Covariance of y_m and w_(m-1)
    # Location of the responses in means and covs
    resp_ids <- (p + 1):(p + m)
    w_cov <- approx_response_cov_recursive(
      ids = 1:m, r_means = means[resp_ids], r_covs = covs[resp_ids, resp_ids],
      n = n, response = response, responses = responses)


    means_in <- c(approx_w$means, mean_m)
    cov_in <- c(covs_m[1:p], w_cov[["cov"]], covs_m[p+1])
    covs_in <- rbind(cbind(approx_w$covs, cov_in[1:(length(cov_in) - 1)]),
                     cov_in
    )

    cov_out <- approx_prod_stats(means = means_in, covs = covs_in, n = n,
                                 response = response, predictors = predictors)

    return(cov_out)
  } else {
    stop("Enter at least 2 responses.")
  }
}





#' Approximate the covariance of one response with an arbitrary product of
#' responses.
#' @param ids Column ids of responses to use. First is taken alone while 2nd to
#'   last are to be multiplied
#' @param r_covs Response covariance matrix
#' @param r_means Response means (vector)
#' @param n Sample size
#' @param response Character, Either "binary" or "continuous"
#' @param responses List of lists with elements of class predictor
#'
#' @return A vector with the approximated covariance, and approximated mean and
#'   variance of the product
#'
#' @examples
#' Y <- bin_data[c("y1", "y2", "y3")]
#' r_means <- colMeans(Y)
#' r_covs <- cov(Y)
#' n <- nrow(Y)
#' responses <- lapply(r_means, new_predictor_binary)
#' approx_response_cov_recursive(ids = c(3, 2, 1), r_covs, r_means, n, responses, response = "binary", verbose = T)
#' w2 <- with(Y, y1 * y2)
#' cov(Y$y3, w2)
#' mean(w2)
#' var(w2)
#'
#'
approx_response_cov_recursive <- function(ids, r_covs, r_means, n, responses,
                                          response, verbose = FALSE) {
  if (length(ids) == 3) {
    if (verbose) {
      if (!any(is.null(names(r_means)))) {
        resp_names <- names(r_means)
      } else {
        resp_names <- 1:(length(r_means))
      }

      message(
        paste("Estimating covariance of", resp_names[ids[1]], "and",
              paste(resp_names[ids[2:length(ids)]], collapse = "")))
    }
    covs0 <- r_covs[ids, ids]
    means0 <- r_means[ids]
    cov_out <- do.call(approx_cov,
                       c(list(means = means0, covs = covs0, n = n, response = response),
                         responses[[ids[1]]])
    )
    return(cov_out)
  } else if (length(ids) > 3) {
    if (verbose) {
      message(paste("More than 3 ids input, approximating recursively"))
    }
    cov_12 <- r_covs[ids[1], ids[2]]

    mean_1 <- r_means[ids[1]]
    var_1  <- r_covs[ids[1], ids[1]]

    mean_2 <- r_means[ids[2]]
    var_2  <- r_covs[ids[2], ids[2]]

    # Approx cov of 1st id and product of 3rd to END
    ids0 <- c(ids[1], ids[3:length(ids)])
    approx_13 <- approx_response_cov_recursive(
      ids = ids0, r_covs = r_covs, r_means = r_means, n = n, responses = responses,
      response = response, verbose = verbose)
    cov_13 <- approx_13[["cov"]]
    # Approx cov of 2nd id and product of 3rd to END and mean of product of 3rd
    # to END
    ids0 <- c(ids[2], ids[3:length(ids)])
    approx_23 <- approx_response_cov_recursive(
      ids = ids0, r_covs = r_covs, r_means = r_means, n = n, responses = responses,
      response = response, verbose = verbose)
    cov_23 <- approx_23[["cov"]]

    var_3  <- median(c(approx_13[["var"]],  approx_23[["var"]]))
    mean_3 <- median(c(approx_13[["mean"]], approx_23[["mean"]]))

    means0 <- c(mean_1, mean_2, mean_3)
    covs0  <- matrix(c(
      var_1,   cov_12, cov_13,
      cov_12, var_2,  cov_23,
      cov_13, cov_23, var_3
    ), ncol = 3)

    cov_out <- do.call(approx_cov,
                       c(list(means = means0, covs = covs0, n = n, response = response),
                         responses[[ids[1]]]))
    return(cov_out)
  } else {
    stop("At least three id values are required.")
  }
}




#' Approximate summary statistics for a product of phenotypes and a set of predictors
#'
#' @param means Vector of means of predictors and the two phenotypes to be multiplied
#' @param covs Covariance matrix of all predictors and the two phenotypes
#' @param n Sample size
#' @param response character. Either "binary" or "continuous".
#' @param predictors a list of elements of class predictor
#' @return A list with the predicted covariance matrix of all predictors and
#'   the product and the means of all predictors and the product.
#'
#' @examples
#' ex_data <- cont_data[c("g", "x", "y3", "y5")]
#' predictors <- list(
#'   new_predictor_snp(maf = mean(ex_data$g) / 2),
#'   new_predictor_normal(mean = mean(ex_data$x), sd = sd(ex_data$x))
#'   )
#' approx_prod_stats(colMeans(ex_data), covs = cov(ex_data), n = nrow(ex_data),
#'   response = "continuous", predictors = predictors)
#' with(ex_data, cov(cbind(g, x, y3y5 = y3 * y5)))
#'
#' @export
approx_prod_stats <- function(means, covs, n, response, predictors) {
  n_preds <- length(means) - 2
  preds <- list()
  for (j in 1:n_preds) {
    params <- list(
      means = means[c(j, n_preds + 1:2)],
      covs = covs[c(j, n_preds + 1:2), c(j, n_preds + 1:2)],
      response = response,
      n = n
    )
    preds[[j]] <- do.call(approx_cov, c(params, predictors[[j]]))
  }
  prod_covs  <- sapply(preds, function(.) .[["cov"]])
  prod_vars  <- sapply(preds, function(.) .[["var"]])
  prod_means <- sapply(preds, function(.) .[["mean"]])

  prod_var <- median(prod_vars)
  prod_mean <- prod_means[[1]]

  means_out <- c(means[1:n_preds], prod_mean)
  cov_out <- rbind(
    cbind(covs[1:n_preds, 1:n_preds], unname(prod_covs)),
    c(prod_covs, prod_var)
  )

  # Put names on
  prod_name <- paste(names(means)[n_preds + 1:2], collapse = "")
  names(means_out)[n_preds + 1]  <- prod_name
  colnames(cov_out)[n_preds + 1] <- prod_name
  rownames(cov_out)[n_preds + 1] <- prod_name

  return(list(means = means_out, covs = cov_out))
}

#' Approximate the mean of Y conditional on X
#' @param means Vector of the mean of X and the mean of Y
#' @param covs Matrix of covariances for X and Y
#' @param response Character. If "binary" truncates means to interval [0, 1].
#'   If "continuous" does not restrict.
#' @param n Sample size
#' @return A list of length 2 consisting of 2 functions that give the
#'   estimated conditional mean and conditional variance of Y as a function of X
#'
approx_conditional <- function(means, covs, response, n) {
  # OLS coefficients for Y = a + bX
  b <- covs[1, 2] / covs[1, 1]
  a <- means[2] - b * means[1]

  if (response == "binary") {
    p_mu <- function(x0) {
      p <- a + b * x0
      p[p > 1] <- 1
      p[p < 0] <- 0
      return(p)
    }

    p_var <- function(x0) p_mu(x0) * (1 - p_mu(x0))

  } else if (response == "continuous") {
    p_mu <- function(x0) {
      p <- a + b * x0
      return(p)
    }
    # Approximate the conditional variance of a continuous Y given X under the
    # OLS assumption of homoscedasticity.
    p_s2 <- (n * means[2]^2 + (n - 1) * covs[2, 2] - a * n * means[2] -
               b * (n * means[1] * means[2] + (n - 1) * covs[1, 2])) / (n - 2)
    p_var <- function(x0) p_s2
  }

  return(list(c_mu = p_mu, c_var = p_var))
}

#' Approximate the partial correlation of Y and Z given X
#' @param covs Covariance matrix of X, Y, and Z.
#' @param cors Correlation matrix of X, Y, and Z.
#' @return Approximated partial correlation of the later two terms given the first
get_pcor <- function(covs, cors = cov2cor(covs)) {
  if (isTRUE(all.equal(cors[1, 2], 1)) | isTRUE(all.equal(cors[1, 3], 1))) {
    rho <- 0
  } else {
    rho <- (cors[2, 3] - cors[1, 2] * cors[1, 3]) /
      (sqrt(1 - cors[1, 2]^2) * sqrt(1 - cors[1, 3]^2))
  }
  return(rho)
}

#' Approximate the covariance of X and Y*Z as well as the mean and variance of Y * Z
#'
#' @examples
#' ex_data <- bin_data[c("g", "y4", "y5")]
#' approx_cov(colMeans(ex_data), cov(ex_data), predictor_type = "discrete",
#'            response = "binary", n = nrow(ex_data),
#'            f = function(x0) dbinom(x0, size = 2, prob = mean(ex_data$g) / 2),
#'            support = 0:2)
#' with(ex_data, cov(g, y4 * y5))
#'
approx_cov <- function(means, covs, predictor_type, response, n, f, ...) {
  # MEAN ##
  pred_mean <- get_mean(means = means[2:3], covs = covs[2:3, 2:3], n = n)

  # COVARIANCE ##
  # Conditional means / variances for phenotype 1 and 2
  c_1 <- approx_conditional(means = means[c(1, 2)], covs = covs[c(1, 2), c(1, 2)],
                            response = response, n = n)

  c_2 <- approx_conditional(means = means[c(1, 3)], covs = covs[c(1, 3), c(1, 3)],
                            response = response, n = n)

  # Partial correlation
  rho <- get_pcor(covs)

  # Predicted product conditional mean
  c_prod_mean <- function(x0) {
    c_1$c_mu(x0) * c_2$c_mu(x0) + rho * sqrt(c_1$c_var(x0) * c_2$c_var(x0)) * (n - 1) / n
  }

  if (predictor_type == "discrete") {
    # Estimated product covariance
    pred_cov <- approx_cov_discrete(c_prod_mean = c_prod_mean,
                                    predictor_mean = means[1],
                                    f = f, ...)
  } else if (predictor_type == "continuous") {
    pred_cov <- approx_cov_continuous(c_prod_mean = c_prod_mean,
                                      predictor_mean = means[1],
                                      f = f, ...)
  } else {
    stop("Invalid predictor_type argument to approx_cov")
  }

  ## VARIANCE ##
  if (response == "binary") {
    pred_var <- pred_mean * (1 - pred_mean) *  n / (n - 1)

  } else if (response == "continuous") {
    # Estimate the conditional variance
    c_prod_var <- function(x0) {
      (c_1$c_mu(x0)^2/c_1$c_var(x0) + c_2$c_mu(x0)^2/c_2$c_var(x0) +
         2 * rho *  c_1$c_mu(x0) * c_2$c_mu(x0) / sqrt(c_1$c_var(x0) * c_2$c_var(x0)) +
         1 + rho^2) * c_1$c_var(x0) * c_2$c_var(x0)
    }

    # Take the expectation(?) of the conditional variance of
    if (predictor_type == "discrete") {
      pred_var <- approx_var_discrete(
        c_prod_mean = c_prod_mean, c_prod_var = c_prod_var,
        prod_mean = pred_mean, n = n, f = f, ...)
    } else if (predictor_type == "continuous") {
      pred_var <- approx_var_continuous(
        c_prod_mean = c_prod_mean, c_prod_var = c_prod_var,
        prod_mean = pred_mean, n = n, f = f, ...)
    }
  }
  return(c(cov = unname(pred_cov), mean = unname(pred_mean), var = unname(pred_var)))
}

approx_cov_discrete <- function(c_prod_mean, predictor_mean, f, support) {
  g <- function(x0) f(x0) * (x0 - predictor_mean) * c_prod_mean(x0)
  pred_cov <- sum(g(support))
  return(pred_cov)
}

approx_cov_continuous <- function(c_prod_mean, predictor_mean, f, lb, ub) {
  g <- function(x0) f(x0) * (x0 - predictor_mean) * c_prod_mean(x0)
  pred_cov <- integrate(g, lower = lb, upper = ub, stop.on.error = FALSE)$value
  return(pred_cov)
}

#' Calculate the mean of the product of X and Y
get_mean <- function(means, covs, n) {
  prod_mean <- means[1] * means[2] + covs[1, 2] * (n - 1) / n
  return(prod_mean)
}

approx_var_discrete <- function(c_prod_mean, c_prod_var, prod_mean, n, f, support) {
  g <- function(x0) {
    (n * f(x0) - 1) * c_prod_var(x0) + n * f(x0) * (c_prod_mean(x0) - prod_mean)^2
  }
  pred_var <- sum(g(support))/(n - 1)
  return(pred_var)
}

approx_var_continuous <- function(c_prod_mean, c_prod_var, prod_mean, n, f, lb, ub) {
  g <- function(x0) {
    # f(x0) * c_prod_var(x0) + f(x0) * c_prod_mean(x0) - prod_mean^2
    (n * f(x0) - 1) * c_prod_var(x0) + n * f(x0) * (c_prod_mean(x0) - prod_mean)^2
  }
  pred_var <- integrate(g, lower = lb, upper = ub, stop.on.error = FALSE)$value / (n - 1)
  return(pred_var)
}


make_permutations <- function(m) {
  # All possible permutations. We only need half of these as the order of the
  # first two columns does not matter.
  if (m < 3) {
    re_perms <- list(1:m)
    return(re_perms)
  }
  all_perms <- gtools::permutations(m, m)
  dups <- which(duplicated(as.matrix(all_perms[, 3:ncol(all_perms)])))
  re_perms <- all_perms[-dups, ]
  re_perms <- split(re_perms, row(re_perms))
  return(re_perms)
}
