#' Constant Model Fit
#'
#' Function that fits a constant line and returns generic model outputs.
#'
#' success = 1 for a successful fit, 0 if optimization failed, and NA if
#' nofit = T. aic, rme, and er are set to NA in case of nofit or failure. pars
#' always equals "er".
#'
#' @param conc Vector of concentration values NOT in log units.
#' @param resp Vector of corresponding responses.
#' @param nofit If nofit = T, returns formatted output filled with missing values.
#'
#' @return List of five elements: success, aic (Aikaike Information Criteria),
#'   rme (root mean square error), er (error parameter), pars (parameter names).
#' @export
#' @importFrom methods is
#' @importFrom stats mad optim
#'
#' @examples
#' fitcnst(c(.1,1,10,100), c(1,2,0,-1))
#' fitcnst(c(.1,1,10,100), c(1,2,0,-1), nofit = TRUE)
fitcnst = function(conc, resp, nofit = F, ...){

  pars <- "er"
  myparams = c("success", "aic", "rme","er")
  #nofit generic output
  if(nofit){
    out = as.list(rep(NA_real_, length(myparams)))
    names(out) = myparams
    out[["success"]] = NA_integer_
    return(out)
  }
  er_est <- if ((rmad <- mad(resp)) > 0) log(rmad) else log(1e-32)

  ###----------------------- Fit the Constant Model -----------------------###
  fit <- optim(er_est,
                tcplObj,
                fname = "cnst",
                method = "Brent",
                lower = er_est - 2,
                upper = er_est + 2,
                control = list(fnscale = -1,
                               reltol = 1e-4,
                               maxit = 500),
                conc = conc,
                resp = resp)
  if (!is(fit, "try-error")) {
    success <- 1L
    er <- fit$par
    aic <- 2 - 2*fit$value

    ## Calculate the rmse for constant
    rme <- sqrt(mean((0 - resp)^2, na.rm = TRUE))
  } else {
    success <- 0L
    er <- NA_real_
    aic <- NA_integer_
    rme <- NA_real_
  }

  return(mget(c("success", "aic", "rme", "er", "pars")))
}
