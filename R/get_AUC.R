#' Calculate AUC
#'
#' Function that calculates AUC for dose-response curves
#'
#' This function takes in a type of model and a set of model parameters,
#' and return the Area under the curve (AUC) between the lower and upper
#' concentration bounds. AUC can be used to compute an efficacy/ potency
#' metric for "active" does-response curves. The model parameters should be
#' entered as a list or a vector. For models operate in log scale (hill and
#' gain-loss), the lower and upper concentration bounds, "ga" (gain AC50)
#' and "la" (loss AC50) will be converted to log-scale if applicable.
#' This version hasn't taken consideration into how to appropriately work with
#' decreasing curves. One of the future works to do.
#'
#'
#' @param fit_method name of the model
#' @param lower lower bound for x
#' @param upper upper bound for x
#' @param ps model parameter vector
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
  return(integrate(get(fit_method),
                      lower,
                      upper,
                      ps = unlist(ps))$value)

}

