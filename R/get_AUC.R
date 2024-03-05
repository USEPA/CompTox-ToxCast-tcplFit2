#' Calculate Area Under the Curve (AUC)
#'
#' Function that calculates the area under the curve (AUC) for dose-response curves.
#'
#' This function takes in a model name and the respective set of model parameters,
#' and returns the area under the curve (AUC) between the specified lower and upper
#' concentration bounds. The AUC can be used to compute an efficacy/potency
#' metric for "active" dose-response curves. For decreasing curves, the AUC returned
#' will be negative. However, users have the option to return a positive AUC in
#' these cases. Model parameters should be entered as a numeric list or vector.
#' Models optimized on the log10-scale (hill and gain-loss), the lower and upper
#' concentration bounds, parameters "ga" (gain AC50) and "la" (loss AC50) will
#' be converted to log10-scale.
#'
#'
#' @param fit_method Name of the model to calculate the area under the curve (AUC) for.
#' @param lower Lower concentration bound, usually is the lowest concentration in the data.
#' @param upper Upper concentration bound, usually is the highest concentration in the data.
#' @param ps Numeric vector (or list) of model parameters for the specified model in `fit_method`.
#' @param return.abs Logical argument, defaults to FALSE.
#' If set to TRUE, the function will convert the negative AUC value and return a positive AUC.
#'
#' @return AUC value (numeric)
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
get_AUC <- function(fit_method, lower, upper, ps, return.abs = FALSE) {

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


  if (return.abs) {
    return(abs(out$value))
  } else {
    return(out$value)
  }

}

