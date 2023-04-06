#' Continuous Hitcalls
#'
#' Wrapper that computes continuous hitcalls for a provided concRespCore input row.
#'
#' indf parameter columns should be NA when not required by fit method. "conc"
#' and "resp" entries should be a single string with values separated by |.
#' Details on indf columns can be found in concRespCore.
#'
#' @param indf Dataframe similar to concRespCore output. Must contain "conc" and "resp"
#'   columns if xs and ys are not provided. Must contain "top", "ac50", "er",
#'   "fit_method", "caikwt", and "mll" columns as well as columns for each
#'   model parameter.
#' @param xs List of concentration vectors that can be provided for speed.
#' @param ys List of response vectors that can be provided for speed.
#' @param newcutoff Vector of new cutoff values to use. Length should be equal
#'   to rows in indf.
#'
#'
#' @return Vector of hitcalls between 0 and 1 with length equal to indf row
#'   number.
#' @export
#' @examples
#' conc <- list(.03, .1, .3, 1, 3, 10, 30, 100)
#' resp <- list(0, .2, .1, .4, .7, .9, .6, 1.2)
#' row <- list(
#'   conc = conc,
#'   resp = resp,
#'   bmed = 0,
#'   cutoff = 1,
#'   onesd = .5,
#'   name = "some chemical",
#'   assay = "some assay"
#' )
#' res <- concRespCore(row, conthits = TRUE)
#' hitcont(res, newcutoff = 0.2)
#'
hitcont <- function(indf, xs = NULL, ys = NULL, newcutoff) {
  # reformat concs and resps
  if (is.null(xs)) {
    xs <- strsplit(indf$conc, "\\|")
    xs <- lapply(xs, as.numeric)
  }
  if (is.null(ys)) {
    ys <- strsplit(indf$resp, "\\|")
    ys <- lapply(ys, as.numeric)
  }

  # correct parameter ordering: used to extract parameters from indf
  parnames <- c("a", "tp", "b", "ga", "p", "la", "q", "er")

  pin <- lapply(1:nrow(indf), function(i) {
    indf[i, parnames]
  })
  hitcall <- mapply(hitcontinner,
    conc = xs, resp = ys, top = indf$top, cutoff = newcutoff, er = indf$er, ps = pin,
    fit_method = indf$fit_method, caikwt = indf$caikwt, mll = indf$mll
  )


  return(hitcall)
}
