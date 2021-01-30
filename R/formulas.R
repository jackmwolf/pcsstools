#' Extract independent variables from a formula
#' @param formula an object of class \code{formula}.
#' @importFrom stats terms
extract_predictors <- function(formula = formula()) {
  terms0 <- terms(formula)

  model_terms <- attributes(terms0)$term.labels
  add_intercept <- attributes(terms0)$intercept

  re <- list(
    predictors = model_terms,
    add_intercept = add_intercept
  )
  return(re)
}

#' Extract dependent variables from a formula as a string
#' @param formula an object of class \code{formula}.
#' @importFrom stats terms
extract_response <- function(formula = formula()) {
  terms0 <- terms(formula)
  response <- as.character(attributes(terms0)$variables)[[2]]
  return(response)
}

parse_response <- function(response = character(), split = character()) {
  terms0 <- strsplit(response, split = split)[[1]]
  # Trim white space at start and end
  terms0 <- sapply(terms0, trimws, simplify = T)

  return(terms0)
}

parse_product <- function(response = character(), ...) {
  terms0 <- unname(parse_response(response = response, split = "\\*", ...))
  return(terms0)
}

parse_or <- function(response = character(), ...) {
  terms0 <- unname(parse_response(response = response, split = "\\|", ...))
  return(terms0)
}

parse_and <- function(response = character(), ...) {
  terms0 <- unname(parse_response(response = response, split = "\\&", ...))
  return(terms0)
}

parse_sum <- function(response = character(), ...) {
  terms0 <- unname(parse_response(response = response, split = "\\+", ...))
  return(terms0)
}

# parse_linear_combo <- function(response = character(), ...) {
#   # Replace "-" with "+ -1 *"
#
#   gsub("^-", "-1 * ", response)
#
#   gsub("[{- }{y%}]", "-1 *", response)
#   terms0 <- parse_sum(response)
# }

#' Check that independent and dependent variables are accounted for through PCSS
#' @param xterms, yterms character vector of model's independent variables or
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


