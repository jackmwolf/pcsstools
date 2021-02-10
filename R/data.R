#' Simulated data with binary phenotypes
#'
#' A dataset containing simulated genetic data with one SNP, a covariate, and
#' 5 binary phenotypes.
#'
#' @format A data frame with 1000 rows and 7 columns:
#' \describe{
#'   \item{g}{Minor allele counts at one site}
#'   \item{x}{A continuous covariate}
#'   \item{y1,y2,y3,y4,y5}{binary phenotypes}
#' }
"bin_data"

#' Simulated data with continuous phenotypes
#'
#' A dataset containing simulated genetic data with one SNP, a covariate, and
#' 5 continuous phenotypes.
#'
#' @format A data frame with 1000 rows and 7 columns:
#' \describe{
#'   \item{g}{Minor allele counts at one site}
#'   \item{x}{A continuous covariate}
#'   \item{y1,y2,y3,y4,y5}{5 continuous phenotypes}
#' }
"cont_data"


#' Simulated example data
#'
#' A dataset containing simulated genetic data with 3 SNPs, 3 continuous 
#' covariates, and 6 continuous phenotypes.
#'
#' @format A data frame with 1000 rows and 12 columns:
#' \describe{
#'   \item{g1,g2,g3}{Minor allele counts at three sites}
#'   \item{x1,x2,x3}{Continuous covariates}
#'   \item{y1,y2,y3}{Continuous phenotypes}
#'   \item{y4,y5,y6}{Binary phenotypes}
#' }
"pcsstools_example"