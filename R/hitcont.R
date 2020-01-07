#' Continuous Hitcalls
#'
#' Wrapper that computes continuous hitcalls for a provided concRespCore input row.
#'
#' indf parameter columns should be NA when not required by fit method. "conc"
#' and "resp" entries should be a single string with values separated by |.
#' Details on indf columns can be found in concRespCore.
#'
#' @param indf Dataframe similar to concRespCore input. Must contain "conc" and "resp"
#'   columns if xs and ys are not provided. Must contain "top", "ac50", "er",
#'   "fit_method", "caikwt", and "mll" columns as well as columns for each
#'   model parameter.
#' @param xs List of concentration vectors that can be provided for speed.
#' @param ys List of response vectors that can be provided for speed.
#' @param newcutoff Vector of new cutoff values to use. Length should be equal
#'   to rows in indf.
#' @param mc.cores Number of cores to use for large dataframes.
#'
#' @import future.apply
#' @import future
#'
#' @return Vector of hitcalls between 0 and 1 with length equal to indf row
#'   number.
#' @export
hitcont = function(indf, xs = NULL, ys = NULL, newcutoff, mc.cores = 1){

  # reformat concs and resps
  if(is.null(xs)){
    xs = strsplit(indf$conc, "\\|")
    xs = lapply(xs, as.numeric)
  }
  if(is.null(ys)){
    ys = strsplit(indf$resp, "\\|")
    ys = lapply(ys, as.numeric)
  }

  #correct parameter ordering: used to extract parameters from indf
  parnames = c("a", "tp", "b", "ga", "p", "la", "q", "er")

  #run hitcontinner for every row of indf
  if(mc.cores > 1){
    plan(multiprocess, workers = mc.cores)
    pin = future_lapply(1:nrow(indf), function(i){indf[i, parnames]})
    hitcall = future_mapply(hitcontinner, conc= xs, resp = ys, top = indf$top, cutoff = newcutoff, er = indf$er, ps = pin,
                            fit_method = indf$fit_method, caikwt = indf$caikwt, mll = indf$mll,
                            future.globals = structure(TRUE, add = c("cnst", "poly1",
                                                                     "poly2", "pow", "exp2", "exp3", "exp4", "exp5", "hillfn", "gnls") ))
    plan("default")
  } else {
    pin = lapply(1:nrow(indf), function(i){indf[i, parnames]})
    hitcall = mapply(hitcontinner, conc= xs, resp = ys, top = indf$top, cutoff = newcutoff, er = indf$er, ps = pin,
                     fit_method = indf$fit_method, caikwt = indf$caikwt, mll = indf$mll)
  }

  return(hitcall)

}
