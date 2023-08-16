#' Calculate AUC
#'
#' Function that calculates AUC for dose-response curves
#'
#' Assume for now we want to calculate AUC for the winning model.
#' The first half of tcplhit2_core selects a winning model and extract the model
#' name, model parameters from model fitting result (result from
#' tcplfit2_core) and populate them into the environment. We will have
#' fit_method - name of the winning model (string) and modpars - a list of
#' model parameters.
#'
#'
#' @param fit_method name of the winning model
#' @param lower lower bound for x
#' @param upper upper bound for x
#' @param ps model parameter vector
#' @param ... placeholder, should be fine without it
#'
#' @return AUC value
#'
#' @export
#' @importFrom stats integrate
#'
#' @examples
#' conc <- c(.03,.1,.3,1,3,10,30,100)
#' modpars <- c(1, 2, er = -3.295307) #a, b, er
#' get_AUC("exp2", min(conc), max(conc), ps = modpars)
#'
get_AUC <- function(fit_method, lower, upper, ps) {

  if (fit_method %in% c("hill", "gnls")) {
    fit_method <- paste0("log", fit_method)
    lower <- log10(lower)
    upper <- log10(upper)
  }

  return(integrate(get(fit_method),
                      lower,
                      upper,
                      ps = unlist(ps))$value)

}
