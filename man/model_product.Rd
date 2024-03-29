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
  responses = NULL,
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
\value{
an object of class \code{"pcsslm"}.

  An object of class \code{"pcsslm"} is a list containing at least the 
  following components:
    \item{call}{the matched call}
    \item{terms}{the \code{terms} object used}
    \item{coefficients}{a \eqn{p x 4} matrix with columns for the 
      estimated coefficient, its standard error, t-statistic and
      corresponding (two-sided) p-value.}
    \item{sigma}{the square root of the estimated variance of the random
      error.}
    \item{df}{degrees of freedom, a 3-vector \eqn{p, n-p, p*}, the
      first being the number of non-aliased coefficients, the last being
      the total number of coefficients.}
    \item{fstatistic}{a 3-vector with the value of the F-statistic with its
      numerator and denominator degrees of freedom.}
    \item{r.squared}{\eqn{R^2}, the 'fraction of variance explained by the 
      model'.}
    \item{adj.r.squared}{the above \eqn{R^2} statistic \emph{'adjusted'},
      penalizing for higher \eqn{p}.}
    \item{cov.unscaled}{a \eqn{p x p} matrix of (unscaled) covariances of the
      \eqn{coef[j], j=1,...p}.}
    \item{Sum Sq}{a 3-vector with the model's Sum of Squares Regression 
      (SSR), Sum of Squares Error (SSE), and Sum of Squares Total (SST).}
}
\description{
\code{model_product} approximates the linear model for the product
  of m phenotypes as a function of a set of predictors.
}
\examples{
ex_data <- pcsstools_example[c("g1", "g2", "g3", "x1", "y4", "y5", "y6")]
head(ex_data)
means <- colMeans(ex_data)
covs <- cov(ex_data)
n <- nrow(ex_data)
predictors <- list(
  g1 = new_predictor_snp(maf = mean(ex_data$g1) / 2),
  g2 = new_predictor_snp(maf = mean(ex_data$g2) / 2),
  g3 = new_predictor_snp(maf = mean(ex_data$g3) / 2),
  x1 = new_predictor_normal(mean = mean(ex_data$x1), sd = sd(ex_data$x1))
)
responses <- lapply(means[c("y4", "y5", "y6")], new_predictor_binary)

model_product(
  y4 * y5 * y6 ~ g1 + g2 + g3 + x1,
  means = means, covs = covs, n = n,
  predictors = predictors, responses = responses, response = "binary"
)

summary(lm(y4 * y5 * y6 ~ g1 + g2 + g3 + x1, data = ex_data))
}
\references{
{

  \insertRef{wolf_using_2021}{pcsstools}
  
}
}
