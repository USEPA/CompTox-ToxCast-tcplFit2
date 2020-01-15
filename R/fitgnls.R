#' Gain-Loss Model Fit
#'
#' Function that fits to f(x) = tp/[(1 + (ga/x)^p )(1 + (x/la)^q )]
#' and returns generic model outputs.
#'
#' Concentrations are converted internally to log10 units and optimized with
#' f(x) = tp/[(1 + 10^(p*(ga-x)) )(1 + 10^(q*(x-la)) )], then ga, la, ga_sd,
#' and la_sd are converted back to regular units before returning.
#' Zero background and increasing initial absolute response are assumed.
#' Parameters are "tp" (top), "ga" (gain AC50), "p" (gain power), "la"
#' (loss AC50),"q" (loss power) and error term "er".
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
#' @param minwidth Minimum allowed distance between gain ac50 and loss ac50 (in
#'   log10 units).
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
#' fitgnls(c(.03,.1,.3,1,3,10,30,100), c(0,.3,1, 2, 2.1, 1.5, .8, .2))
fitgnls = function(conc, resp, bidirectional = TRUE, verbose = FALSE, nofit = F, minwidth = 1.5){

  logc = log10(conc)
  fenv <- environment()

  pars <- paste0(c("tp", "ga", "p", "la", "q", "er"))
  sds <- paste0(c("tp", "ga", "p", "la", "q", "er"), "_sd")
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

  rmds <- tapply(resp, logc, median)
  if(!bidirectional) mmed = rmds[which.max(rmds)] else mmed = rmds[which.max(abs(rmds))] #shortened this code
  mmed_conc <- as.numeric(names(mmed)) #fixed this bug

  resp_max <- max(resp)
  resp_min <- min(resp)
  logc_min <- min(logc)
  logc_max <- max(logc)

  er_est <- if ((rmad <- mad(resp)) > 0) log(rmad) else log(1e-32)

  ###--------------------- Fit the Gain-Loss Model ----------------------###
  ## Starting parameters for the Gain-Loss Model
  # cind <- (ceiling(length(meds)/2) + 1):length(meds)
  g <- c(mmed, # top (tp)
         mmed_conc - 0.5, # gain logAC50 (ga)
         1.2, # gain hill coefficient (p)
         # mmed_conc - 0.99 + minwidth + .01, # loss logAC50 (la), a little farther than min width
         mmed_conc - 0.5 + minwidth + .01, # loss logAC50 (la), start with tight gnls ranges
         5, # loss hill coefficient (q)
         er_est) # logSigma (er)
  if (g[1] == 0) g[1] <- 0.1
  ## Generate the bound matrices to constrain the model.
  #                tp   ga   p   la   q   er
  Ui <- matrix(c(  1,   0,   0,   0,   0,   0,
                  -1,   0,   0,   0,   0,   0,
                   0,   1,   0,   0,   0,   0,
                   0,  -1,   0,   0,   0,   0,
                   0,   0,   1,   0,   0,   0,
                   0,   0,  -1 ,  0,   0,   0,
                   0,   0,   0,   1,   0,   0,
                   0,   0,   0,  -1,   0,   0,
                   0,   0,   0,   0,   1,   0,
                   0,   0,   0,   0,  -1,   0,
                   0,  -1,   0,   1,   0,   0),
                byrow = TRUE, nrow = 11, ncol = 6)
  if(!bidirectional) {
    bnds <- c(0, -1.2*resp_max, # tp bounds
               logc_min - 1, -(logc_max + .5), # ga bounds
               0.3, -8, # p bounds
               logc_min - 1, -(logc_max + 2), # la bounds
               0.3, -8, # q bounds
               minwidth) # la-ga >= minwidth
  } else {
    val <- 1.2*max(abs(resp_min),abs(resp_max))
    bnds <- c(-val,-val, # tp bounds
               logc_min - 1, -(logc_max + 0.5), # ga bounds
               0.3, -8, # p bounds
               logc_min - 1, -(logc_max + 2), # la bounds
               0.3, -8, # q bounds
               minwidth) # la-ga >= minwidth
  }

  Ci <- matrix(bnds, nrow = 11, ncol = 1)

  ## Optimize the gnls model
  fit <- try(constrOptim(g,
                          tcplObj,
                          ui = Ui,
                          ci = Ci,
                          mu = 1e-6,
                          method = "Nelder-Mead",
                          control = list(fnscale = -1,
                                         reltol = 1e-10,
                                         maxit = 6000),
                          conc = logc,
                          resp = resp,
                          fname = "loggnls"),
              silent = !verbose)

  ## Generate some summary statistics
  if (!is(fit, "try-error")) { # Gain-loss fit the data
    if(verbose) cat("gnls >>>",fit$counts[1],fit$convergence,"\n")

    success <- 1L
    aic <- 2*length(fit$par) - 2*fit$value # 2*length(fit$par) - 2*fit$value
    mapply(assign,
           c(pars),
           fit$par,
           MoreArgs = list(envir = fenv))

    ## Calculate rmse for gnls
    modl = loggnls(fit$par, logc)
    rme <- sqrt(mean((modl - resp)^2, na.rm = TRUE))

    #output ga, la in regular units
    ga = 10^ga
    la = 10^la

    ## Calculate the sd for the gnls parameters
    fit$cov <- try(solve(-hessian(tcplObj,
                                   fit$par,
                                   conc = logc,
                                   resp = resp,
                                   fname = "loggnls")),
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
        #use taylor's theorem to approximate sd's in change of units
        #(only valid when sd's are much smaller than ln(10))
        ga_sd = ga*log(10)*ga_sd
        la_sd = la*log(10)*la_sd
      }

    } else { # Could not invert gnls Hessian

      cov <- 0L
      mapply(assign,
             c(sds),
             NA_real_,
             MoreArgs = list(envir = fenv))

    }

  } else { # Gain-loss did not fit the data

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
