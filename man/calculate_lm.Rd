% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_models.R
\name{calculate_lm}
\alias{calculate_lm}
\title{Calculate a linear model}
\usage{
calculate_lm(means, covs, n, add_intercept = FALSE)
}
\arguments{
\item{means}{a vector of means of all model predictors and the response with
the last element the response mean.}

\item{covs}{a matrix of the covariance of all model predictors and the
response with the order of rows/columns corresponding to the order of
\code{means}.}

\item{n}{sample size}

\item{add_intercept}{logical. If \code{TRUE} adds an intercept to the model.}
}
\description{
\code{calculate_lm} describes the linear model of the last listed variable
in \code{means} and \code{covs} as a function of all other variables in
\code{means} and \code{covs}.
}
\examples{
ex_data <- cont_data[c("g", "x", "y1")]
means <- colMeans(ex_data)
covs <- cov(ex_data)
n <- nrow(ex_data)
calculate_lm(means = means, covs = covs, n = n, add_intercept = TRUE)
# Compare results to...
mod <- lm(y1 ~ 1 + g + x, data = ex_data)
summary(mod)$coef

}