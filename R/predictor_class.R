#' Create an object of class "predictor"
#'
#' @param f a function that gives the probability mass/distribution function of
#'   a random variable.
#' @param predictor_type a character describing the random variable. Either "discrete"
#'   or "continuous".
#' @param lb,ub if \code{predictor_type == "continuous"} double giving the
#'   lower/upper bound of the pdf \code{f}.
#' @param support if \code{predictor_type == "discrete"} vector of the support of
#'   the pmf for \code{f}.
#'
#' @return an object of class predictor.
#' @examples
#' new_predictor(
#'   f = function(x0) dnorm(x0, mean = 0, sd = 1),
#'   predictor_type = "continuous", lb = -Inf, ub = Inf
#' )
#' @seealso \code{\link{new_predictor_normal}},
#'   \code{\link{new_predictor_snp}} and \code{\link{new_predictor_binary}}.
#'
#'
#' @export
new_predictor <- function(f = function() {}, predictor_type = character(), lb, ub, support) {
  x <- list(f = f, predictor_type = predictor_type)

  if (!missing(lb) & !missing(ub)) {
    x <- c(x, lb = lb, ub = ub)
  } else if (!missing(support)) {
    x <- c(x, list(support = support))
  }

  validate_predictor(x)
  class(x) <- "predictor"
  return(x)
}

validate_predictor <- function(x) {
  stopifnot(is.function(x$f))
  stopifnot(is.character(x$predictor_type))

  if (x$predictor_type == "discrete") {
    stopifnot(is.vector(x$support))
  } else if (x$predictor_type == "continuous") {
    stopifnot(is.double(x$lb), is.double(x$ub))
  } else {
    stop('Invalid argument to predictor_type (Currently supports "discrete" and "continuous"')
  }
}


#' Shortcut to create a predictor object for a continuous variable
#'
#' \code{new_predictor_normal} calls \code{new_predictor}
#'
#' @param mean predictor mean (double).
#' @param sd predictor standard deviation (double)
#'
#' @return an object of class predictor.
#'
#' @importFrom stats dnorm
#'
#' @examples
#' new_predictor_normal(mean = 10, sd = 1)
#' @export
new_predictor_normal <- function(mean, sd) {
  new_predictor(
    f = function(x0) dnorm(x0, mean = mean, sd = sd),
    predictor_type = "continuous", lb = -Inf, ub = Inf
  )
}

#' Shortcut to create a predictor object for a SNP's minor allele counts
#'
#' \code{new_predictor_snp} calls \code{new_predictor}
#'
#' @param maf minor allele frequency
#'
#' @return an object of class predictor.
#'
#' @examples
#' new_predictor_snp(maf = 0.3)
#' @export
new_predictor_snp <- function(maf) {
  new_predictor(
    f = function(x0) dbinom(x0, size = 2, prob = maf),
    predictor_type = "discrete", support = 0:2
  )
}

#' Shortcut to create a predictor object for a binary variable
#'
#' \code{new_predictor_binary} calls \code{new_predictor}
#'
#' @param p probability of success (predictor mean)
#'
#' @return an object of class predictor.
#'
#' @importFrom stats dbinom
#'
#' @examples
#' new_predictor_binary(p = 0.75)
#' @export
new_predictor_binary <- function(p) {
  new_predictor(
    f = function(x0) dbinom(x0, size = 1, prob = p),
    predictor_type = "discrete", support = 0:1
  )
}
