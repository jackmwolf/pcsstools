% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiplication_estimation.R
\name{get_pcor}
\alias{get_pcor}
\title{Approximate the partial correlation of Y and Z given X}
\usage{
get_pcor(covs, cors = cov2cor(covs))
}
\arguments{
\item{covs}{Covariance matrix of X, Y, and Z.}

\item{cors}{Correlation matrix of X, Y, and Z.}
}
\value{
Approximated partial correlation of the later two terms given the 
  first
}
\description{
Approximate the partial correlation of Y and Z given X
}
