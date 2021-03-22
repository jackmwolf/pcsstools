#' Extract independent variables from a formula
#' @param formula an object of class \code{formula}.
#' @importFrom stats terms
#' @return A list with a character vector of all predictors and a logical
#'   value indicating whether the model includes an intercept term.
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
#' @return a character vector of all responses
extract_response <- function(formula = formula()) {
  terms0 <- terms(formula)
  response <- as.character(attributes(terms0)$variables)[[2]]
  return(response)
}

#' Guess the function that is applied to a set of responses
#'
#' \code{guess_response} takes a character vector of the dependent variable
#' from a \code{formula} object and identifies which function separates the
#' individual variables that make up the response. It then returns the
#' \code{model_*} function to model the appropriate response using PCSS.
#' 
#' @param response character. Output of \code{extract_response}.
#' 
#' @return A character. Either \code{"model_combo"}, \code{"model_product"}, 
#' \code{"model_or"}, \code{"model_and"}, or \code{"model_singular"}.
#'   
guess_response <- function(response = character()) {
  
  seps <- c("\\+", "\\*", "\\|", "\\&")
  names(seps) <- c("model_combo", "model_product", "model_or", "model_and")
  
  
  response_types <- sapply(seps, function(.x) grepl(.x, response))
  response_type <- names(response_types[response_types])
  
  if (length(response_type) > 1) {
    stop("Multiple response functions detected")
  } else if (length(response_type) == 0) {
    # No seperator detected. Assume response is one variable of its own
    func <- "model_singular"
  } else {
    func <- response_type
  }
  
  return(func)

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


