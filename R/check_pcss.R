#' Check that independent and dependent variables are accounted for through PCSS
#' @param xterms,yterms character vector of model's independent variables or
#'   variables combined to the dependent variable
#' @param pcssterms character vector of variables with provided PCSS
#' @param pcsstype character describing the PCSS being checked. Either
#'   \code{"means"}, \code{"covs"}, \code{"predictors"}, or
#'   \code{"responses"}.
#'    
check_terms <- function(xterms, yterms, pcssterms, pcsstype) {
  
  if (pcsstype %in% c("means", "covs")) {
    missing.x <- setdiff(xterms, pcssterms)
    missing.y <- setdiff(yterms, pcssterms)
    
  } else if (pcsstype %in% c("predictors")) {
    missing.x <- setdiff(xterms, pcssterms)
    missing.y <- character(0)
    
  } else if (pcsstype %in% c("responses")) {
    missing.x <- character(0)
    missing.y <- setdiff(yterms, pcssterms)
    
  }
  
  if (length(missing.x) != 0) {
    stop(paste0("Independent variable(s) not listed in `", pcsstype, "`: ", 
                paste(missing.x, collapse = ", ")))
  }
  if (length(missing.y) != 0) {
    stop(paste0("Dependent variable(s) not listed in `", pcsstype, "`: ", 
                paste(missing.y, collapse = ", ")))
  }
  
  return(TRUE)
}


check_terms_combo <- function(xterms, yterms, means, covs) {
  check_terms(xterms, yterms, names(means), "means")
  check_terms(xterms, yterms, dimnames(covs)[[1]], "covs")
  return(TRUE)
}

check_terms_logical <- function(xterms, yterms, means, covs, 
                                predictors) {
  check_terms(xterms, yterms, names(means), "means")
  check_terms(xterms, yterms, dimnames(covs)[[1]], "covs")
  check_terms(xterms, yterms, names(predictors), "predictors")
  
  return(TRUE)
}


check_terms_product <- function(xterms, yterms, means, covs, 
                                predictors, responses) {
  check_terms(xterms, yterms, names(means), "means")
  check_terms(xterms, yterms, dimnames(covs)[[1]], "covs")
  check_terms(xterms, yterms, names(predictors), "predictors")
  if (length(yterms) > 2) {
    check_terms(xterms, yterms, names(responses), "responses")
  }
  
  return(TRUE)
}


