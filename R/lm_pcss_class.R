

#' Print an object of class summary.lm_pcss
#'
#' Prints a linear model fit through pre-computed summary statistics
#'
#' @method print summary.lm_pcss
#'
#' @param x an object of class \code{"summary.lm_pcss"}
#' @param digits the number of significant digits to use when printing.
#' @param symbolic.cor logical. If TRUE, print the correlations in a symbolic
#'   form (see \link[stats]{symnum}) rather than as numbers.
#' @param signif.stars logical. If \code{TRUE}, ‘significance stars’ are printed
#'    for each coefficient.
#' @param ... further arguments passed to or from other methods.
#'
#'
#' @importFrom stats naprint pf printCoefmat symnum
#'
#' @export
#'
print.summary.lm_pcss <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  symbolic.cor = x$symbolic.cor, signif.stars = getOption("show.signif.stars"),
                                  ...) {
  cat("Model approximated using Pre-Computed Summary Statistics.\n")
  cat("\nCall:\n", paste(deparse(x$call),
    sep = "\n",
    collapse = "\n"
  ), "\n", sep = "")
  df <- x$df
  rdf <- df[2L]

  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    if (nsingular <- df[3L] - df[1L]) {
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
        sep = ""
      )
    } else {
      cat("\nCoefficients:\n")
    }
    coefs <- x$coefficients
    if (any(aliased <- x$aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4, dimnames = list(
        cn,
        colnames(coefs)
      ))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs,
      digits = digits, signif.stars = signif.stars,
      na.print = "NA", ...
    )
  }
  cat("\nResidual standard error:", format(signif(
    x$sigma,
    digits
  )), "on", rdf, "degrees of freedom")
  cat("\n")
  if (nzchar(mess <- naprint(x$na.action))) {
    cat("  (", mess, ")\n", sep = "")
  }
  if (!is.null(x$fstatistic)) {
    cat("Multiple R-squared: ", formatC(x$r.squared,
      digits = digits
    ))
    cat(
      ",\tAdjusted R-squared: ", formatC(x$adj.r.squared,
        digits = digits
      ), "\nF-statistic:", formatC(x$fstatistic[1L],
        digits = digits
      ), "on", x$fstatistic[2L], "and",
      x$fstatistic[3L], "DF,  p-value:", format.pval(pf(x$fstatistic[1L],
        x$fstatistic[2L], x$fstatistic[3L],
        lower.tail = FALSE
      ),
      digits = digits
      )
    )
    cat("\n")
  }
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1L) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2),
          nsmall = 2,
          digits = digits
        )
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}
