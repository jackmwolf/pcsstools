% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiplication_estimation.R
\name{model_product}
\alias{model_product}
\title{Approximate a linear model for a product using PCSS}
\usage{
model_product(
  formula,
  n,
  means,
  covs,
  predictors,
  responses,
  response = "continuous",
  ...
)
}
\arguments{
\item{formula}{an object of class \code{formula} whose dependent variable is
a combination of variables and  \code{*} operators. All model terms
must be accounted for in \code{means} and \code{covs}.}

\item{n}{sample size.}

\item{means}{named vector of predictor and response means.}

\item{covs}{named matrix of the covariance of all model predictors and the
responses.}

\item{predictors}{named list of objects of class \code{predictor}}

\item{responses}{named list of objects of class \code{predictor} 
corresponding to all terms being multiplied in the response. Can be 
left \code{NULL} if only multiplying two terms}

\item{response}{character. Describe distribution of all product terms.
Either \code{"continuous"} or \code{"binary"}. If \code{"binary"} 
different approximations of product means and variances are used.}

\item{...}{additional arguments}
}
\description{
\code{model_product} approximates the linear model for the product
  of m phenotypes as a function of a set of predictors.
}
\examples{
ex_data <- bin_data[c("g", "x", "y1", "y2", "y3")]
head(ex_data)
means <- colMeans(ex_data)
covs <- cov(ex_data)
n <- nrow(ex_data)
predictors <- list(
  g = new_predictor_snp(maf = mean(ex_data$g) / 2),
  x = new_predictor_normal(mean = mean(ex_data$x), sd = sd(ex_data$x))
)
responses <- lapply(means[3:length(means)], new_predictor_binary)

model_product(
  y1 * y2 * y3 ~ g + x, means = means, covs = covs, n = n, 
  predictors = predictors, responses = responses
)

summary(lm(y1 * y2 * y3 ~ g + x, data = ex_data))

}