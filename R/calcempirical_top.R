#' Empirically locate top and top concentration
#'
#' Locate the top and the concentration at which top occurs empirically by sampling points within the experimental concentration range.
#'
#' @param conc Vector of concentrations
#' @param ps Vector of parameters, must be in order: a, tp, b, ga, p, la, q, er
#' @param fname Model function name (equal to model name except hill which
#'   uses "hillfn")
#' @param precision Number of points sampled on the base-10 logarithmic scaled experimental concentration range. (Defaults to 100.)
#'
#' @return A named list of the top and top concentration. The elements of the list are:
#'  \itemize{
#'    \item top - empirical top
#'    \item x_top - concentration at which empirical top occurs
#'  }
#' @export
#'
#' @examples
#' conc = c(.03,.1,.3,1,3,10,30,100)
#' ps = c(1,2,1,2,2)
#' fname = "gnls"
#' calcempirical_top(conc, ps, fname)
calcempirical_top = function(conc, ps, fname, precision = 100) {

  # replace untreated controls with pseudo-value
  if (any(conc == 0)) warning("Data contains untreated controls (conc = 0). A pseudo value replaces -Inf after log-transform.  The pseudo value is set to one log-unit below the lowest experimental `conc`.")

  logc <- log10(conc)
  logc_temp <- replace(logc, logc == -Inf, sort(unique(logc)[2]-1))
  # generate concentration sequence over entire experimental concentration range
  conc_seq <- 10^seq(from = min(logc_temp), to = max(logc_temp), length.out = precision)
  fit = do.call(fname, list(ps, conc_seq))
  # calculate the largest change from baseline (y = 0) of the modeled fit
  top = fit[which.max(abs(fit))] # empirical top
  if(all(is.na(top))) top = NA # if top returns a numeric(empty), replace with NA
  x_top = conc_seq[which.max(abs(fit))] # empirical x_top
  if(all(is.na(x_top))) x_top = NA # if x_top returns a numeric(empty), replace with NA

  return(list("top" = top, "x_top" = x_top))
}
