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
