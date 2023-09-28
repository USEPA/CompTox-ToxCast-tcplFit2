#' Calculate AUC After Hit-calling
#'
#' Function that calculates AUC for dose-response curves after hit-calling
#'
#' This function calculates AUC value for the winning model selected after
#' hit-calling. This is a wrapper function for get_AUC. Designated to take
#' one-row output from tcplhit2_core and parse the model details, then
#' pass to get_AUC to obtain AUC for the winning model.
#'
#' @param hit_results output from tcplhit2_core
#'
#' @return AUC value of the winning model
#'
#' @export
#' @importFrom stringr str_split
#'
#' @examples
#'
#' conc <- c(.03, .1, .3, 1, 3, 10, 30, 100)
#' resp <- c(0, .2, .1, .4, .7, .9, .6, 1.2)
#' params <- tcplfit2_core(conc, resp, .8)
#' output <- tcplhit2_core(params, conc, resp, 0.8, 0.5)
#' post_hit_AUC(output)
#'
#'
post_hit_AUC <- function(hit_results, ...) {

  # parameter list
  param <- c("a","tp","b","ga","p", "la", "q", "er")

  # get concentrations
  conc <- as.numeric(str_split(hit_results[1,"conc"],"\\|")[[1]])
  # get fitted values of the winning model
  modpars <- unlist(hit_results[1, param])
  modpars <- modpars[!is.na(modpars)]
  # get the winning model name
  fit_method <- hit_results[["fit_method"]]
  out <- get_AUC(fit_method, min(conc), max(conc), ps = modpars, ...)

  # return AUC
  return(out)

}
