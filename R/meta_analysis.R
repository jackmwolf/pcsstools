#' Approximate a linear model using PCSS from multiple sources
#' 
#' \code{meta_pcsslm} pools summary statistics across sources and then 
#'   approximates a lineaer model using the pooled summary statistics as input.
#'   
#' @note \code{meta_pcsslm} currently does not support models that require
#'   response multiplication including products, or logical disjunctions or
#'   conjunctions. 
#'
#' @inheritParams pcsslm
#' @inheritParams pool_meta_stats
#' 
#' @examples
#' 
#' ## Prepare PCSS
#' data("pcsstools_example")
#' ex_data <- pcsstools_example[c("g1", "x1", "x2", "y1", "y2", "y3", "Source")]
#' pcss_list <- list(
#'   "Study_1" = list(
#'     means = colMeans(ex_data[ex_data$Source == 1, ]),
#'     covs = cov(ex_data[ex_data$Source == 1, ]),
#'     n = nrow(ex_data[ex_data$Source == 1, ])
#'   ),
#'   
#'   "Study_2" = list(
#'     means = colMeans(ex_data[ex_data$Source == 2, ]),
#'     covs = cov(ex_data[ex_data$Source == 2, ]),
#'     n = nrow(ex_data[ex_data$Source == 2, ])
#'   ),
#'   
#'   "Study_3" = list(
#'     means = colMeans(ex_data[ex_data$Source == 3, ]),
#'     covs = cov(ex_data[ex_data$Source == 3, ]),
#'     n = nrow(ex_data[ex_data$Source == 3, ])
#'   )
#' )
#' ## Covariate Adjustment
#' meta_pcsslm(y1 ~ g1 + x1 + x2, pcss_list = pcss_list)
#' 
#' ## Principal Component Analysis
#' meta_pcsslm(y1 + y2 + y3 ~ g1 + x1, pcss_list = pcss_list, comp = 1)
#' 
#' @export
#' 
meta_pcsslm <- function(formula, pcss_list = list(), ...) {
  pcss <- pool_meta_stats(pcss_list)
  
  if ("comp" %in% names(list(...))) {
    model_func <- "model_prcomp"
  } else {
    model_func <- guess_response(extract_response(formula))
  }
  
  if (!model_func %in% c("model_singular", "model_prcomp", "model_combo")) {
    stop("meta_pcsslm currently does not support models that require response multiplication.")
  }
  
  cl <- match.call()
  args <- c(append(list(formula = formula), pcss), list(...))
  
  re <- do.call(pcsslm, args)
  re$call <- cl
  return(re)
}

#' Aggregate PCSS from multiple studies into PCSS for pooled data
#'
#' @param pcss_list A named list of lists of precomputed summary statistics.
#'  Each element should correspond to a study or source and match the format
#'  described for argument \code{pcss} in \code{\link{pcsslm}}. Names,
#'  content, and orderings must be consistent across studies.
#'  
#' @return A list of precomputed summary statistics for the pooled data of
#'   the format described for argument \code{pcss} in \code{\link{pcsslm}}.
#' 
#' @author Nathan Ryder  
pool_meta_stats <- function(pcss_list) {
  
  if (
    is.null(names(pcss_list)) | 
    length(pcss_list) != length(unique(names(pcss_list)))
  ) {
    stop("Argument pcss_list must be a named list with unique names")
  }
  group_names <- names(pcss_list)
  
  # Pull PCSS from pcss_list
  mean_vectors <- lapply(pcss_list, function(.x) .x$means)
  cov_matrices <- lapply(pcss_list, function(.x) .x$covs)
  ns           <- lapply(pcss_list, function(.x) .x$n)
  
  if (inherits(ns, "list")) { ns <- unlist(ns) }
  
  ## make objects
  num_groups <- length(pcss_list)
  means_mat <- do.call(cbind, mean_vectors)
  num_vars <- nrow(means_mat)
  
  ## Check that study PCSS have the same names
  mean_names <- rownames(means_mat)
  all_cov_names <- cbind(sapply(cov_matrices, rownames), sapply(cov_matrices, colnames))
  
  if (any(apply(all_cov_names, 2, function(.x) !identical(mean_names, .x)))) {
    stop("Cannot pool PCSS across studies: Mean and covariance names and orderings must match.")
  }
    
  ## array for SS later
  means_array <- replicate(num_groups, means_mat)
  names(means_array)
  
  ## overall means and repeating matrix
  overall_means <- (means_mat %*% ns)/sum(ns)
  overall_means_mat <- overall_means[, rep(1, num_vars)]
  
  ## sum the products of differences of sample and overall means for all pairs 
  ## of variables
  ss_group <- rowSums(
    sapply(
      group_names, 
      function(group_name) {
        ns[group_name] *
          matrix(mean_vectors[[group_name]] - overall_means, ncol = 1) %*%
          matrix(mean_vectors[[group_name]] - overall_means, nrow = 1)
        }, 
      simplify = "array"), 
    dims = 2
    )
  
  ## sum the unscaled covariances of all variable pairs
  ss_within <- rowSums(
    sapply(
      group_names, function(group_name) {
        cov_matrices[[group_name]] * (ns[group_name] - 1)
        }, 
      simplify = "array"), 
    dims = 2)
  
  ## sum the SS's to get the total, then scale
  full_covs <- (ss_group + ss_within) / (sum(ns) - 1)
  colnames(full_covs) <- rownames(full_covs) <- mean_names
  
  ## return full summary stats as PCSS
  pcss_pooled <- list(
    means = drop(overall_means),
    covs = full_covs,
    n = sum(ns)
  )
  
  # If applicable, pool predictors
  if (!is.null(pcss_list[[1]]$predictors)) {
    predictors <- pool_predictors(pcss_list)
    pcss_pooled$predictors <- predictors
  }
  
  return(pcss_pooled)
}


#' Aggregate predictors from multiple studies
#' 
#' @inheritParams pool_meta_stats
#' 
pool_predictors <- function(pcss_list) {
  
  predictors <- lapply(pcss_list, function(.x) .x$predictors)
  predictor_names <- names(predictors[[1]])
  
  # Treat all predictors as mixture distributions weighted by study n
  weights <- sapply(pcss_list, function(.x) .x$n)
  weights <- weights / sum(weights)
  
  if (
    any(sapply(predictors, function(.x) names(.x) != predictor_names))
  ) {
    stop("Cannot pool PCSS predictors across studies: Predictor names must match.")
  }
  
  
  predictors <- lapply(
    predictor_names,
    function(predictor_name) {
      
      # Check that predictor types are consistent across studies
      predictor_types <- sapply(
        predictors, function(.x) .x[[predictor_name]]$predictor_type
      )
      if (any(predictor_types != predictor_types[[1]])) {
        stop("Cannot pool PCSS predictors across studies: predictor_type must be consistent across studies.")
      } else {
        predictor_type <- predictor_types[[1]]
      }
      
      
      # Evaluate f under each study and take weighted average
      fs <- lapply(predictors, function(.x) .x[[predictor_name]]$f)
      f <- function(x) sum(weights * sapply(fs, do.call, list(x)))
      
      
      
      if (predictor_type == "discrete") {
        # Take union of support across all studies
        support <- Reduce(
          union, 
          lapply(predictors, function(.x) .x[[predictor_name]]$support)
        )
        
        new_predictor(
          f = f,
          predictor_type = "discrete",
          support = support 
        )
      } else {
        lb <- min(sapply(predictors, function(.x) .x[[predictor_name]]$lb))
        ub <- max(sapply(predictors, function(.x) .x[[predictor_name]]$ub))
        
        new_predictor(
          f = f,
          predictor_type = "continuous",
          lb = lb,
          ub = ub
        )
      }
      
      }
    )
  names(predictors) <- predictor_names
  
  return(predictors)
  
}