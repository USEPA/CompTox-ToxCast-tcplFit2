#' Polynomial 2 (Quadratic) Model Fit
#'
#' Function that fits to f(x) = a*(x/b + x^2/b^2) and returns generic model outputs.
#'
#' Zero background and monotonically increasing absolute response are assumed.
#' Parameters are "a" (y scale), "b" (x scale), and error term "er".
#' success = 1 for a successful fit, 0 if optimization failed, and NA if
#' nofit = T. cov = 1 for a successful hessian inversion, 0 if it fails, and NA
#' if nofit = T. aic, rme, modl, parameters, and parameter sds are set to
#' NA in case of nofit or failure.
#'
#' @param conc Vector of concentration values NOT in log units.
#' @param resp Vector of corresponding responses.
#' @param bidirectional If TRUE, model can be positive or negative; if FALSE, it
#'   will be positive only.
#' @param verbose If TRUE, gives optimization and hessian inversion details.
#' @param nofit If nofit = T, returns formatted output filled with missing values.
#'
#' @importFrom methods is
#' @importFrom numDeriv hessian
#' @importFrom stats constrOptim median
#'
#' @return Named list containing: success, aic (Aikaike Information Criteria),
#'   cov (success of covariance calculation), rme (root mean square error),
#'   modl (vector of model values at given concentrations),
#'   parameters values, parameter sd (standard deviation) estimates, pars
#'   (vector of parameter names), sds (vector of parameter sd names).
#' @export
#'
#' @examples
#' fitpoly2(c(.03,.1,.3,1,3,10,30,100), c(0,.01,.1, .1, .2, .5, 2, 8))
fitpoly2 = function(conc, resp, bidirectional = TRUE, verbose = FALSE, nofit = F){

  fenv <- environment()
  #initialize myparams
  pars <- paste0(c("a", "b", "er"))
  sds <- paste0(c("a", "b", "er"), "_sd")
  myparams = c("success", "aic", "cov", "rme", "modl", pars, sds, "pars", "sds")

  #returns myparams with appropriate NAs
  if(nofit){
    out = as.list(rep(NA_real_, length(myparams)))
    names(out) = myparams
    out[["success"]] = out[["cov"]] = NA_integer_
    out[["pars"]] = pars
    out[["sds"]] = sds
    return(out)
  }

  #median at each conc, for multi-valued responses
  rmds <- tapply(resp, conc, median)
  #get max response and corresponding conc
  if(!bidirectional) mmed = rmds[which.max(rmds)] else mmed = rmds[which.max(abs(rmds))] #shortened this code
  mmed_conc <- as.numeric(names(mmed)) #fixed this bug

  resp_max <- max(resp)
  resp_min <- min(resp)
  conc_min <- min(conc)
  conc_max <- max(conc)

  er_est <- if ((rmad <- mad(resp)) > 0) log(rmad) else log(1e-16)

  ###--------------------- Fit the Model ----------------------###
  ## Starting parameters for the Model
  a0 = mmed #use largest response with desired directionality
  if(a0 == 0) a0 = .01  #if 0, use a smallish number
  g <- c(a0/2, # y scale (a); set to run through the max resp at the max conc
         conc_max, # x scale (b); set to max conc
         er_est) # logSigma (er)

  ## Generate the bound matrices to constrain the model.
  #                a   b    er
  Ui <- matrix(c( 1,   0,   0,
                 -1,   0,   0,
                  0,   1,   0,
                  0,  -1,   0),
                byrow = TRUE, nrow = 4, ncol = 3)

  if(!bidirectional){
    bnds <- c(0, -1e8*abs(a0), # a bounds (always positive)
              1e-8*conc_max, -1e8*conc_max) # b bounds (always increasing)
  } else {
    bnds <- c(-1e8*abs(a0), -1e8*abs(a0), # a bounds (positive or negative)
              1e-8*conc_max, -1e8*conc_max) # b bounds (always increasing or always decreasing)
  }

  Ci <- matrix(bnds, nrow = 4, ncol = 1)

  ## Optimize the model
  fit <- try(constrOptim(g,
                          tcplObj,
                          ui = Ui,
                          ci = Ci,
                          mu = 1e-6,
                          method = "Nelder-Mead",
                          control = list(fnscale = -1,
                                         reltol = 1e-10,
                                         maxit = 6000),
                          conc = conc,
                          resp = resp,
                          fname = "poly2"),
              silent = !verbose)

  ## Generate some summary statistics
  if (!is(fit, "try-error")) { # The model fit the data
    if(verbose) cat("poly2 >>>",fit$counts[1],fit$convergence,"\n")

    success <- 1L
    aic <- 2*length(fit$par) - 2*fit$value # 2*length(fit$par) - 2*fit$value
    mapply(assign,
           c(pars),
           fit$par,
           MoreArgs = list(envir = fenv))

    ## Calculate rmse for gnls
    modl <- poly2(fit$par,conc)
    rme <- sqrt(mean((modl - resp)^2, na.rm = TRUE))

    ## Calculate the sd for the gnls parameters
    fit$cov <- try(solve(-hessian(tcplObj,
                                   fit$par,
                                   conc = conc,
                                   resp = resp,
                                   fname = "poly2")),
                    silent = !verbose)

    if (!is(fit$cov, "try-error")) { # Could invert gnls Hessian

      cov <- 1L
      diag_sqrt <- suppressWarnings(sqrt(diag(fit$cov)))
      if (any(is.nan(diag_sqrt))) {
        mapply(assign,
               sds,
               NaN,
               MoreArgs = list(envir = fenv))
      } else {
        mapply(assign,
               sds,
               diag_sqrt,
               MoreArgs = list(envir = fenv))
      }

    } else { # Could not invert gnls Hessian

      cov <- 0L
      mapply(assign,
             c(sds),
             NA_real_,
             MoreArgs = list(envir = fenv))

    }

  } else { # Curve did not fit the data

    success <- 0L
    aic <- NA_real_
    cov <- NA_integer_
    rme <- NA_real_
    modl = NA_real_

    mapply(assign,
           c(pars, sds),
           NA_real_,
           MoreArgs = list(envir = fenv))

  }

  return(mget(myparams))

}
