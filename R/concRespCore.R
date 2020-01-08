#' Concentration Response Core
#'
#' Core of concentration response curve fitting for pvalue based cutoff. This
#' function calls tcplfit2_core to get curve fits, and then tcplhit2_core to
#' perform the hitcalling.
#'
#' @param row A named list that must include:
#'   \itemize{
#'     \item conc - list of concentrations (not in log units)
#'     \item resp - list of corresponding responses
#'     \item bmed - median of noise estimate.
#'     \item cutoff - noise cutoff
#'     \item onesd - 1 standard deviation of the noise (for bmd calculation)
#'   }
#'   Other elements (usually identifiers, like casrn) of row will be attached to
#'   the final output.
#' @param fitmodels Vector of model names to use.
#' @param conthits conthits = T uses continuous hitcalls, otherwise they're
#'   discrete.
#' @param aicc aicc = T uses corrected AIC to choose winning method; otherwise
#'   regular AIC.
#' @param force.fit If TRUE force the fitting to proceed even if there are no points
#'   outside of the bounds (default FALSE)
#' @param bidirectional If TRUE allow fitting to happen in both directions (default TRUE)
#' @param verbose  If TRUE, write extra output from tcplfit2_core (default FALSE)
#' @param do.plot If TRUE, create a plot in the tcplfit2_core function (default FALSE)
#' @param return.details If TRUE, return the hitcalling details and the summary, if FALSE (default), just return the summary
#'
#' @return A list of two elements. The first (summary) is teh output from tcplhit2_core. The second, params is the
#' output from tcplfit2_core
#' a dataframe of one row containing
#' @export
#'
#' @examples
#' conc = list(.03,.1,.3,1,3,10,30,100)
#' resp = list(0,.2,.1,.4,.7,.9,.6, 1.2)
#' row = list(conc = conc, resp = resp, bmed = 0, cutoff = 1, onesd = .5,name="some chemical")
#' concRespCore(row, conthits = TRUE)
#' concRespCore(row, aicc = TRUE)
#'
concRespCore <- function(row,
                         fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
                                       "exp4", "exp5"),
                         conthits = T,
                         aicc = F,
                         force.fit = FALSE,
                         bidirectional = TRUE,
                         verbose = FALSE,
                         do.plot = FALSE,
                         return.details=FALSE) {
  # variable binding to pass cmd checks
  bmed <- cutoff <- onesd <- NULL
  # row needs to include cutoff and bmed
  # unpack row into the local environment, for ease: sample_id, dtxsid, casrn, name, time, pathway, size, con, resp
  list2env(row,envir=environment())
  resp = unlist(resp)
  conc = unlist(conc)

  # prepare input
  resp = resp - bmed
  conc = conc[!is.na(resp)]
  resp = resp[!is.na(resp)]
  identifiers = row[!names(row) %in% c("conc", "resp", "bmed", "onesd", "cutoff")]

  # run the fits
  params <- tcplfit2_core(conc, resp, cutoff, force.fit = conthits, bidirectional = bidirectional, fitmodels = fitmodels,
                     verbose=verbose, do.plot=do.plot)

  # calculate the hitcall
  summary <- tcplhit2_core(params,conc,resp,cutoff,onesd,bmed,conthits,aicc,identifiers)
  if(return.details) return(list(summary=summary,all.models=params))
  else return(summary)
}


