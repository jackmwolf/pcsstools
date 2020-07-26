
#' Approximate a linear model for a series of logical AND statements
#'
#' \code{approx_and} approximates the linear model for the a conjunction
#'   of m phenotypes as a function of a set of predictors.
#'
#'  @param means vector of predictor and response means with the last \code{m}
#'   means being the means of \code{m} binary responses to combine in a
#'   logical and statement.
#' @param covs a matrix of the covariance of all model predictors and the
#'   responses with the order of rows/columns corresponding to the order of
#'   \code{means}.
#' @param n sample size.
#' @param predictors list of objects of class \code{predictor} corresponding
#'   to the order of the predictors in \code{means}.
#' @param add_intercept logical. Should the linear model add an intercept term?
#'
#' @examples
#' ex_data <- bin_data[c("g", "x", "y1", "y2")]
#' head(ex_data)
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)
#' predictors <- list(
#'   new_predictor_snp(maf = mean(ex_data$g) / 2),
#'   new_predictor_normal(mean = mean(ex_data$x), sd = sd(ex_data$x))
#' )
#'
#' approx_and(means = means, covs = covs, n = n, predictors = predictors,
#'   add_intercept = TRUE)
#' y1_and_y2 <- with(ex_data, y1 & y2)
#' coef(summary(lm(y1_and_y2 ~ g + x + 1, data = ex_data)))
#'
#' @export
approx_and <- function(means, covs, n, predictors, add_intercept = TRUE, verbose = FALSE) {
  m <- length(means) - length(predictors)
  p <- length(means) - m

  # Generate responses' pmfs
  r_means <- means[(p + 1) : (p + m)]
  responses <- lapply(r_means, new_predictor_binary)

  approx0 <- approx_mult_prod(means = means, covs = covs, n = n,
                              response = "binary", responses = responses,
                              predictors = predictors, verbose = verbose)
  do.call(calculate_lm, c(approx0, n = n, add_intercept = TRUE))
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
#'
#' @examples
#' ex_data <- bin_data[c("g", "x", "y1", "y2")]
#' head(ex_data)
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)
#' predictors <- list(
#'   new_predictor_snp(maf = mean(ex_data$g) / 2),
#'   new_predictor_normal(mean = mean(ex_data$x), sd = sd(ex_data$x))
#' )
#'
#' approx_or(means = means, covs = covs, n = n, predictors = predictors,
#'   add_intercept = TRUE)
#' coef(summary(lm(y1 | y2 ~ 1 + g + x, data = ex_data)))
approx_or <- function(means, covs, n, predictors, add_intercept = TRUE, verbose = FALSE) {
  # Model "y1 or y2 or ..." via "not(not y1 and not y2 and ...)"
  m <- length(means) - length(predictors)
  p <- length(means) - m
  
  not_means <- c(means[1:p], 1 - means[(p + 1):(p + m)])
  phi <- c(rep(1, p), rep(-1, m))
  not_covs <- t(covs * phi) * phi
  
  # Generate responses' pmfs
  r_means <- not_means[(p + 1) : (p + m)]
  responses <- lapply(r_means, new_predictor_binary)
  
  approx_not_and <- approx_mult_prod(
    means = not_means, covs = not_covs, n = n,
    response = "binary", responses = responses,
    predictors = predictors,
    verbose = verbose
  )
  
  tau <- c(rep(1, p), -1)
  out_covs <- t(approx_not_and$covs * tau) * tau 
  out_means <- approx_not_and$means * tau + c(rep(0, p), 1)
  
  calculate_lm(means = out_means, covs = out_covs, n = n, add_intercept = add_intercept)
}

# Below functions are defunct/outdated -----------------------------------------

#' Approximate a linear model for a series of logical OR statements
#'
#' \code{approx_or} approximates the linear model for a disjunction of m
#'   phenotypes as a function of a set of predictors.
#'
#' @param means vector of predictor and response means with the last \code{m}
#'   means being the means of \code{m} binary responses to combine in a
#'   logical or statement.
#' @param covs a matrix of the covariance of all model predictors and the
#'   responses with the order of rows/columns corresponding to the order of
#'   \code{means}.
#' @param n sample size.
#' @param predictors list of objects of class \code{predictor} corresponding
#'   to the order of the predictors in \code{means}.
#' @param add_intercept logical. Should the linear model add an intercept term?
#'
#' @examples
#' # 2 Responses ----------------------------------------------------
#' ex_data <- bin_data[c("g", "x", "y1", "y2")]
#' head(ex_data)
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)
#' predictors <- list(
#'   new_predictor_snp(maf = mean(ex_data$g) / 2),
#'   new_predictor_normal(mean = mean(ex_data$x), sd = sd(ex_data$x))
#' )
#'
#' approx_or(means = means, covs = covs, n = n, predictors = predictors,
#'   add_intercept = TRUE)
#' y1_or_y2 <- with(ex_data, y1 | y2)
#' coef(summary(lm(y1_or_y2 ~ g + x + 1, data = ex_data)))
#'
#' # 3 Responses ----------------------------------------------------
#' ex_data <- test_data_bin[c("g", "x", "y1", "y2", "y3")]
#' head(ex_data)
#' means <- colMeans(ex_data)
#' covs <- cov(ex_data)
#' n <- nrow(ex_data)
#' predictors <- list(
#'   new_predictor_snp(maf = mean(ex_data$g) / 2),
#'   new_predictor_normal(mean = mean(ex_data$x), sd = sd(ex_data$x))
#' )
#' approx_or(means = means, covs = covs, n = n, predictors = predictors,
#'  verbose = TRUE, add_intercept = TRUE)
#' yor <- with(ex_data, (y1 | y2 | y3))
#' coef(summary(lm(yor ~ 1 + g + x, data = ex_data)))
#'
approx_or_OLD <- function(means, covs, n, predictors, add_intercept = TRUE, verbose = FALSE) {
  m <- length(means) - length(predictors)
  if (m == 2) {
    approx_or_2_OLD(means = means, covs=  covs, n = n, predictors = predictors,
                add_intercept = add_intercept, verbose = verbose)
  } else if (m == 3) {
    approx_or_3_OLD(means = means, covs=  covs, n = n, predictors = predictors,
                add_intercept = add_intercept, verbose = verbose)
  } else {
    stop("approx_or currently only supports up to three responses.")
  }
}

#' Approximate a linear model for a logical OR statement with 2 disjuncts
#'
#' \code{approx_or_2} approximates the linear model for the response
#'   "y1 or y2" as a function of a set of predictors.
#'
#' @param means vector of predictor and response means with the last two means
#'   the means of two binary responses to combine in a logical or statement
#' @param covs a matrix of the covariance of all model predictors and the
#'   responses with the order of rows/columns corresponding to the order of
#'   \code{means}.
#' @param n sample size.
#' @param predictors list of objects of class \code{predictor} corresponding
#'   to the order of the predictors in \code{means}.
#' @param add_intercept logical. Should the linear model add an intercept term?
#'
approx_or_2_OLD <- function(means, covs, n, predictors, verbose = FALSE, add_intercept = TRUE) {
  # Number of responses
  m <- length(means) - length(predictors)
  p <- length(means) - m
  if (m != 2) {
    stop("approx_or_2 requires exactly 2 responses")
  }

  # response means/covariance matrix
  r_means <- means[(p + 1) : (p + m)]
  r_covs  <- covs[((p + 1) : (p + m)), ((p + 1) : (p + m))]

  # Generate responses' pmfs
  responses <- lapply(r_means, new_predictor_binary)

  # Approximate predictor/product covariances
  approx0 <- approx_mult_prod(means = means, covs = covs, n = n, response = "binary",
                           predictors = predictors, responses = responses)
  cov0 <- approx0$covs

  # Calculate response/product covariances: cov(y1, y1*y2) and cov(y2, y1*y2)
  resp_covs <- approx0$means[p + 1] * (1 - r_means) * n / (n - 1)

  cov_B <- rbind(covs[1:p, (p + 1) : (p + m)], resp_covs)

  cov_C <- cbind(t(cov_B), r_covs)

  big_cov <- rbind(cbind(cov0, cov_B), cov_C)
  # Reorder columns
  reorder_ids <- c(1:p, p+2:3, p+1)
  big_cov <- big_cov[reorder_ids, reorder_ids]
  big_means <- c(approx0$means, r_means)[reorder_ids]

  # weights for linear regression (y1 OR y2 = y1 + y2 - y1*y2)
  phi <- c(1, 1, -1)

  # linear regression of linear combination of y1, y2, y1*y2
  calculate_lm_combo(means = big_means, covs = big_cov, n = n, phi = phi, m = 3,
                     add_intercept = add_intercept)

}

#' Approximate a linear model for a logical OR statement with 3 disjuncts
#'
#' \code{approx_or_3} approximates the linear model for the response
#'   "y1 or y2 or y3" as a function of a set of predictors.
#'
#' @param means vector of predictor and response means with the last two means
#'   the means of three binary responses to combine in a logical or statement
#' @param covs a matrix of the covariance of all model predictors and the
#'   responses with the order of rows/columns corresponding to the order of
#'   \code{means}.
#' @param n sample size.
#' @param predictors list of objects of class \code{predictor} corresponding
#'   to the order of the predictors in \code{means}.
#' @param add_intercept logical. Should the linear model add an intercept term?
#'
approx_or_3_OLD <- function(means, covs, n, predictors, verbose = FALSE, add_intercept = TRUE) {
  # Number of responses
  m <- length(means) - length(predictors)
  p <- length(means) - m
  if (m != 3) {
    stop("approx_or_3 requires exactly 2 responses")
  }

  # response means/covariance matrix
  r_means <- means[(p + 1) : (p + m)]
  r_covs  <- covs[((p + 1) : (p + m)), ((p + 1) : (p + m))]

  # Generate responses' pmfs
  responses <- lapply(r_means, new_predictor_binary)

  # 2 term products ------------------------------------------------------------
  cov_B_2 <- matrix(nrow = p, ncol = m)
  means_2 <- numeric(length = m)
  vars_2  <- numeric(length = m)
  prod_names_2 <- character(length = m)
  for (i in m:1) {
    exclude_id <- p + i
    approx2 <-
      approx_mult_prod(
        means = means[-exclude_id], covs = covs[-exclude_id, -exclude_id],
        n = n, response = "binary", predictors = predictors,
        responses = responses[-i]
      )
    cov_B_2[, m - i + 1] <- approx2$covs[1:p, p + 1]
    means_2[m - i + 1] <- approx2$means[p + 1]
    vars_2[m - i + 1] <- approx2$covs[p + 1, p + 1]
    prod_names_2[m - i + 1] <- names(approx2$means[p + 1])
  }
  names(means_2) <- prod_names_2




  # 3 term product -------------------------------------------------------------
  cov_B_3 <- matrix(nrow = p, ncol = 1)
  means_3 <- numeric(length = 1)
  vars_3 <- numeric(length = 1)
  prod_names_3 <- character(length = 1)
  approx3 <- approx_mult_prod(
    means = means, covs = covs, n = n, response = "binary",
    predictors = predictors, responses = responses,
    verbose = verbose)
  cov_B_3[, 1] <- approx3$covs[1:p, p + 1]
  means_3[1] <- approx3$means[p + 1]
  vars_3[1] <- approx3$covs[p + 1, p + 1]
  prod_names_3[1] <- names(approx3$means[p + 1])
  names(means_3) <- prod_names_3

  cov_B <- cbind(cov_B_2, cov_B_3)
  colnames(cov_B) <- c(prod_names_2, prod_names_3)

  # Compare cov_B to
  # with(ex_data, cov(cbind(g, x, y1*y2, y1*y3, y2*y3, y1*y2*y3)))[1:p, (p+1):(p+m+1)]

  # cov(y1, y2*y3), cov(y2, y1*y3), and cov(y3, y1*y2)
  id_list <- list(c(1, 2, 3), c(2, 1, 3), c(3, 1, 2))
  mix_covs_2 <-
    sapply(id_list,
           function(id0) {
             approx_response_cov_recursive(
               ids = id0, r_covs = r_covs, r_means = r_means, n = n,
               responses = responses, response = "binary"
             )[["cov"]]
           })

  # cov(y1, y1*y2), cov(y1, y1*y3), cov(y2, y1*y2), cov(y2, y2*y3),
  # cov(y3, y1*y3), and cov(y3, y2*y3)
  # Using Cov(yi, yi * yj) = n/n-1 * E(yi*yj) (1 - E(yi))
  auto_covs_2 <-
    means_2[c(1, 2, 1, 3, 2, 3)] * (1 - r_means[c(1, 1, 2, 2, 3, 3)]) * n / (n - 1)

  cov_C_2 <- matrix(
    c(auto_covs_2[1:2], mix_covs_2[1],
      auto_covs_2[3], mix_covs_2[2], auto_covs_2[4],
      mix_covs_2[3], auto_covs_2[5:6]),
    nrow = 3, ncol = 3, byrow = TRUE)


  cov_C_3 <- means_3[1] * (1 - r_means) * n / (n - 1)

  cov_C <- cbind(cov_C_2, cov_C_3)
  colnames(cov_C) <- c(prod_names_2, prod_names_3)
  # Compare cov_C
  # with(ex_data, cov(cbind(y1, y2, y3, y1*y2, y1*y3, y2*y3, y1*y2*y3)))[1:3, 4:7]\

  # Covariance of with(ex_data, cov(cbind(y1*y2, y1*y3, y2*y3)))
  # Approximate whole thing and then correct diagonal with exact information
  cov_D_1 <- means_3[1] - (means_2 %*% t(means_2))
  diag(cov_D_1) <- vars_2

  # Covariance of y1y2y3 with y1y2, y1y3, and y2y3
  cov_D_2 <- means_3[1] * (1 - means_2) * n / (n-1)

  cov_D <- rbind(cbind(cov_D_1, cov_D_2),
                 c(cov_D_2, vars_3[1]))

  big_cov <- rbind(
    cbind(covs, rbind(cov_B, cov_C)),
    cbind(t(cov_B), t(cov_C), cov_D)
  )

  big_means <- c(means, means_2, means_3)

  # Regression -----------------------------------------------------------------

  # weights for linear regression:
  # (y1 OR y2 = y1 + y2 + y3 - y1*y2 - y1*y3 - y2*y3 + y1*y2*y3)
  phi <- c(rep(1, 3), rep(-1, 3), 1)

  # linear regression of linear combination of y1, y2, y1*y2
  calculate_lm_combo(means = big_means, covs = big_cov, n = n, phi = phi, m = 7,
                     add_intercept = add_intercept)

}
