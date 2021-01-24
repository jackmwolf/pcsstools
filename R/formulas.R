#' Extract independent variables from a formula
#' @param formula an object of class \code{formula}.
#' @param all_vars character vector of possible terms
#' @importFrom stats terms
extract_predictors <- function(formula = formula(), all_vars = character()) {
  terms0 <- terms(formula)
  
  model_terms <- attributes(terms0)$term.labels
  add_intercept <- attributes(terms0)$intercept
  
  if (any(!(model_terms %in% all_vars))) {
    stop("Independent variable in formula without summary stats provided")
  }
  
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

parse_response <- function(response = character(), all_vars = character(), split = character()) {
  terms0 <- strsplit(response, split = split)[[1]]
  # Trim white space at start and end
  terms0 <- sapply(terms0, trimws, simplify = T)
  
  if (any(!(terms0 %in% all_vars))) {
    stop("Dependent variable in formula without summary stats provided")
  }
  
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

parse_linear_combo <- function(response = character(), ...) {
  # Replace "-" with "+ -1 *"
  
  gsub("^-", "-1 * ", response)
  
  gsub("[{- }{y%}]", "-1 *", response)
  terms0 <- parse_sum(response)
}