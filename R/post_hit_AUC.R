#' Calculate AUC After Hit-calling
#'
#' Function that calculates AUC for dose-response curves after hit-calling
#'
#' This function calculates AUC value for the winning model selected after
#' hit-calling.
#'
#'
#' @param output output from concRespCore or tcplhit2_core
#'
#' @return AUC value of the winning model
#'
#' @export
#' @importFrom stringr str_split
#'
#' @examples
#'
#' conc <- list(.03, .1, .3, 1, 3, 10, 30, 100)
#' resp <- list(0, .2, .1, .4, .7, .9, .6, 1.2)
#' row <- list(conc = conc,
#'             resp = resp,
#'             bmed = 0,
#'             cutoff = 1,
#'             onesd = .5,
#'             name = "some chemical",
#'             assay = "some assay")
#' output <- concRespCore(row)
#'
#' post_hit_AUC(output)
#'
#'
post_hit_AUC <- function(output) {

  # in case user select return.detail = TRUE in concRespCore,
  # the output will be a list of 2
  if (is.null(nrow(output))) {
    output <- output[[1]]
  }

  # parameter list
  param <- c("a","tp","b","ga","p", "la", "q", "er")

  # get concentrations
  conc <- as.numeric(str_split(output[1,"conc"],"\\|")[[1]])
  # get fitted values of the winning model
  modpars <- as.numeric(output[1, param])
  modpars <- modpars[!is.na(modpars)]
  # get the winning model name
  fit_method <- output[["fit_method"]]

  # return AUC
  return(get_AUC(fit_method, min(conc), max(conc), ps = modpars))

}
