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
#' @param return.abs Logical argument, if TRUE, returns the absolute value of the AUC. Defaults to FALSE.
#' @param use.log Logical argument, defaults to FALSE. By default, the function estimates AUC with
#' concentrations in normal unit. If set to TRUE, will use concentration in log10-scale for
#' estimating AUC.
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
get_AUC <- function(fit_method, lower, upper, ps, return.abs = FALSE, use.log = FALSE) {

  if (use.log) {
    message("log10-scale concentration is used for AUC estimation")
    # convert concentration to log10 scale
    lower <- log10(lower)
    upper <- log10(upper)

    # special cases - models that also need to convert parameters
    if (fit_method %in% c("hill", "gnls", "exp4", "exp5")) {
      if (fit_method == "gnls") {
        ps["ga"] <- log10(ps[["ga"]])
        ps["la"] <- log10(ps[["la"]])
        } else {
          # for loghill, logexp4 and logexp5
          ps["ga"] <- log10(ps[["ga"]])
        }
    }
    # add "log" to the function name
    fit_method <- paste0("log", fit_method)
  } else {
    # function name of hill in normal unit is "hillfn"
    if (fit_method == "hill") {fit_method <- paste0(fit_method, "fn")}
  }

  # calculate AUC value
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

# re-parameterization of functions other than hill and gnls with x on log10 scale

# poly1 in normal/raw scale
# f(x) = a*x
logpoly1 <- function(ps, x) {
  return(ps[1]*(10^x))
}

# poly2 in normal/raw scale
# f(x) = a*(x/b + x^2/b^2)
logpoly2 <- function(ps, x) {
  x0 = (10^x)/ps[2]
  return(ps[1]*(x0 + x0*x0))
}

# power in normal/raw scale
# f(x) = a*x^p
logpow = function(ps,x){
  #a = ps[1], p = ps[2]
  return(ps[1]*10^(x*ps[2])  )
}

# exp2 in normal/raw scale
# f(x) = a*(e^{(x/b)}- 1)
logexp2 = function(ps,x){
  #a = ps[1], b = ps[2]
  return(ps[1]*(exp((10^x)/ps[2]) - 1)  )
}

# exp3 in normal/raw scale
# f(x) = a*(e^{(x/b)^p} - 1)
logexp3 = function(ps,x){
  #a = ps[1], b = ps[2], p = ps[3]
  return(ps[1]*(exp((10^x/ps[2])^ps[3]) - 1)  )
}

# exp4 in normal/raw scale
# f(x) = tp*(1-2^{(-x/ga)})
logexp4 = function(ps,x){
  # for computational convenience, both x and ga are converted to log10-scale
  #tp = ps[1], ga = ps[2]
  return(ps[1]*(1-2^(-10^(x-ps[2])))  )
}

# exp5 in normal/raw scale
# f(x) = tp*(1-2^{(-(x/ga)^p)})
logexp5 = function(ps,x){
  # for computational convenience, both x and ga are converted to log10-scale
  #tp = ps[1], ga = ps[2], p = ps[3]
  return(ps[1]*(1-2^(-10^((x-ps[2])*ps[3])))  )
}
