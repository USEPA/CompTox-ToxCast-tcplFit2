#' Calculate Area Under the Curve (AUC)
#'
#' Function that calculates area under the curve (AUC) for dose-response curves
#'
#' This function takes in the model name and the respective set of model parameters,
#' and returns the area under the curve (AUC) between the specified lower and upper
#' concentration bounds. AUC can be used to compute an efficacy/ potency
#' metric for "active" does-response curves. The model parameters should be
#' entered as a list or a vector. For models operate in log10-scale (hill and
#' gain-loss), the lower and upper concentration bounds, "ga" (gain AC50)
#' and "la" (loss AC50) will be converted to log10-scale if applicable.
#'
#'
#' @param fit_method Name of the model to calculate area under the curve (AUC) for
#' @param lower Lower concentration bound, usually is the lowest concentration in the data
#' @param upper Upper concentration bound, usually the highest concentration in the data
#' @param ps Vector (or list) of model parameters for the specified model in fit_method
#'
#' @return AUC value
#'
#' @export
#' @importFrom stats integrate
#'
#' @examples
#' conc <- c(.03,.1,.3,1,3,10,30,100)
#' fit_method <- "gnls"
#' modpars <- list(tp = 1.023, ga = 2.453, p = 1.592,
#'                 la = 4288.993, q = 5.770, er = -3.295)
#' get_AUC("exp2", min(conc), max(conc), ps = modpars)
#'
get_AUC <- function(fit_method, lower, upper, ps) {

  # special cases - models in log scale, hill and gnls
  if (fit_method %in% c("hill", "gnls")) {

    # convert concentration to log scale
    lower <- log10(lower)
    upper <- log10(upper)

    if (fit_method == "hill") {

      # convert parameter to log scale
      # paste log to the fname
      ps["ga"] <- log10(ps[["ga"]])
      fit_method <- paste0("log", fit_method)

    } else if (fit_method == "gnls") {

      # convert parameter to log scale
      # paste log to the fname
      ps["ga"] <- log10(ps[["ga"]])
      ps["la"] <- log10(ps[["la"]])
      fit_method <- paste0("log", fit_method)
    }
  }

  # calculate and return AUC
  out <- integrate(get(fit_method),
                   lower,
                   upper,
                   ps = unlist(ps))

  return(out$value)

}

