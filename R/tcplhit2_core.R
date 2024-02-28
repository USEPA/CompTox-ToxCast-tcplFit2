#' Hitcalling Function
#'
#' Core of hitcalling function. This method chooses the winning model from tcplfit2_core,
#' extracts the top and ac50, computes the hitcall, and calculates bmd/bmdl/bmdu among other
#' statistics. Nested model selection is used to choose between poly1/poly2, then
#' the model with the lowest AIC (or AICc) is declared the winner. Continuous
#' hitcalls requires tcplfit2_core to be run with force.fit = TRUE and "cnst" never to
#' be chosen as the winner.
#'
#' @param params The output from tcplfit2_core
#' @param conc list of concentrations (not in log units)
#' @param resp list of corresponding responses
#' @param bmed median of noise estimate. Default 0
#' @param cutoff noise cutoff
#' @param onesd 1 standard deviation of the noise (for bmd calculation)
#' @param bmr_scale bmr scaling factor. Default = 1.349
#' @param conthits conthits = TRUE uses continuous hitcalls, otherwise they're
#'   discrete. Default TRUE
#' @param aicc aicc = TRUE uses corrected AIC to choose winning method; otherwise
#'   regular AIC. Default FALSE
#' @param identifiers A one-row data frame containing identifiers of the concentration-response profile,
#'   such as the chemical name or other identifiers, and any assay identifiers. The column names
#'   identify the type of value. This can be NULL. The values will be included in the output
#'   summary data frame
#' @param bmd_low_bnd Multiplier for bmd lower bound, must be between 0 and 1
#'   (excluding 0). A value of 0.1 would require the bmd is no lower than a
#'   of 1/10th of the lowest tested concentration (i.e. lower threshold).
#'   If the bmd is less than the threshold, the bmd and its confidence
#'   interval will be censored and shifted right.
#' @param bmd_up_bnd Multiplier for the bmd upper bound, must be greater than
#'   or equal to 1. A value of 10 would require the bmd is no larger than 10
#'   times the highest tested concentration (i.e. upper threshold). If the bmd
#'   is greater than the threshold, the bmd and its confidence interval will be
#'   censored and shifted left.
#'
#' @return A list of with the detailed results from all of the different model fits.
#' The elements of summary are:
#'   \itemize{
#'     \item any elements of the identifiers input
#'     \item n_gt_cutoff - number of data points above the cutoff
#'     \item cutoff - noise cutoff
#'     \item fit_method - curve fit method
#'     \item top_over_cutoff - top divided by cutoff
#'     \item rmse - RMSE of the data points around the best model curve
#'     \item a - fitting parameter methods: exp2, exp3, poly1, poly2, pow
#'     \item b - fitting parameter methods: exp2, exp3, ploy2
#'     \item p - fitting parameter methods: exp3, exp5, gnls, hill, pow
#'     \item q - fitting parameter methods: gnls,
#'     \item tp - top of the curve
#'     \item ga - ac50 for the rising curve in a gnls model or the Hill model
#'     \item la - ac50 for the falling curve in a gnls model
#'     \item er - fitted error term for plotting error bars
#'     \item bmr - benchmark response; level at which bmd is calculated = onesd*bmr_scale default bmr_scale is 1.349
#'     \item bmd - benchmark dose, curve value at bmr
#'     \item bmdl - lower limit on the bmd
#'     \item bmdu - upper limit on the bmd
#'     \item caikwt - one factor used in calculating the continuous hitcall. It is calcalated from the formula
#'       = exp(-aic(cnst)/2)/(exp(-aic(cnst)/2) + exp(-aic(fit_method)/2)) and measures how much lower the
#'       selected method AIC is than that for the constant model
#'     \item mll - another factor used in calcualting the continuous hitcall = length(modpars) - aic(fit_method)/2
#'     \item hitcall - the final hitcall, a value ranging from 0 to 1
#'     \item top - curve top
#'     \item ac50 - curve value at 50\% of top, curve value at cutoff
#'     \item lc50 - curve value at 50\% of top corresponding to the loss side of the gain-loss curve
#'     \item ac5 - curve value at 5\% of top
#'     \item ac10 - curve value at 10\% of top
#'     \item ac20 - curve value at 20\% of top
#'     \item acc - curve value at cutoff
#'     \item ac1sd - curve value at 1 standard deviation
#'     \item conc - conc string separated by |'s
#'     \item resp - response string separated by |'s
#'   }
#' @export
#'
tcplhit2_core <- function(params, conc, resp, cutoff, onesd,bmr_scale = 1.349, bmed = 0, conthits = TRUE, aicc = FALSE, identifiers = NULL, bmd_low_bnd = NULL, bmd_up_bnd = NULL) {
  # initialize parameters to NA
  a <- b <- tp <- p <- q <- ga <- la <- er <- top <- ac50 <- ac50_loss <- ac5 <- ac10 <- ac20 <- acc <- ac1sd <- bmd <- NA_real_
  bmdl <- bmdu <- caikwt <- mll <- NA_real_

  # get error distribution
  errfun = params[["errfun"]]
  if (is.null(errfun))
    warning("'errfun' is missing in the output from tcplfit2_core. 'errfun' tracks the error distribution assumed for the model fits and should be provided in 'params' list. -- see tcplfit2_core help for additional details.")

  # get aics and degrees of freedom
  aics <- sapply(params$modelnames, function(x) {
    params[[x]][["aic"]]
  })
  dfs <- sapply(params$modelnames, function(x) {
    length(params[[x]][["pars"]])
  })
  aics <- aics[!is.na(aics)]
  if (sum(!is.na(aics)) == 0) {
    # if all fits failed, use none for method
    fit_method <- "none"
    rmse <- NA_real_
  } else {
    # use nested chisq to choose between poly1 and poly2, remove poly2 if it fails.
    # pvalue hardcoded to .05
    aics <- nestselect(aics, "poly1", "poly2", dfdiff = 1, pval = .05)
    dfs <- dfs[names(dfs) %in% names(aics)]
    # it's useful to keep original aics so we can extract loglikelihoods for nested models (above) and mll (below)
    if (aicc) saics <- aics + 2 * dfs * (dfs + 1) / (length(resp) - dfs - 1) else saics <- aics
    if (conthits) {
      # if all fits, except the constant fail, use none for the fit method
      # when continuous hit calling is in use
      if(sum(!is.na(aics)) == 1 & "cnst" %in% names(aics[!is.na(aics)])){
        fit_method <- "none"
        rmse <- NA_real_
      }else{
        # get aikaike weight of winner (vs constant) for cont hitcalls
        # never choose constant as winner for cont hitcalls
        nocnstaics <- saics[names(saics) != "cnst"]
        fit_method <- names(nocnstaics)[which.min(nocnstaics)]
        caikwt <- exp(-saics["cnst"] / 2) / (exp(-saics["cnst"] / 2) + exp(-saics[fit_method] / 2))
        if (is.nan(caikwt)) {
          term <- exp(saics["cnst"] / 2 - saics[fit_method] / 2)
          if (term == Inf) {
            caikwt <- 0
          } else {
            caikwt <- 1 / (1 + term)
          }
          # caikwt <- 1
        }
      }
    } else {
      fit_method <- names(saics)[which.min(saics)]
    }
    # if the fit_method is not reported as 'none' the obtain model information
    if(fit_method!="none"){
      fitout <- params[[fit_method]]
      rmse <- fitout$rme
      modpars <- fitout[fitout$pars]
      list2env(fitout, envir = environment()) # put all parameters in environment
    }
  }
  n_gt_cutoff <- sum(abs(resp) > cutoff)

  # compute discrete or continuous hitcalls
  if (fit_method == "none") {
    hitcall <- 0
  } else if (conthits) {
    mll <- length(modpars) - aics[[fit_method]] / 2
    hitcall <- hitcontinner(conc, resp, top, cutoff, er,
      ps = unlist(modpars), fit_method,
      caikwt = caikwt, mll = mll, errfun = errfun
    )
  } else {
    hitcall <- hitloginner(conc, resp, top, cutoff, ac50)
  }

  if(is.nan(hitcall)){
    hitcall <- 0
  }


  bmr <- onesd * bmr_scale # magic bmr is default 1.349
  if (hitcall > 0) {

    # fill ac's; can put after hit logic
    ac5 <- acy(.05 * top, modpars, type = fit_method) # note: cnst model automatically returns NAs
    ac10 <- acy(.1 * top, modpars, type = fit_method)
    ac20 <- acy(.2 * top, modpars, type = fit_method)
    acc <- acy(sign(top) * cutoff, c(modpars,top = top), type = fit_method)
    ac1sd <- acy(sign(top) * onesd, modpars, type = fit_method)
    bmd <- acy(sign(top) * bmr, modpars, type = fit_method)

    # get bmdl and bmdu
    bmdl <- bmdbounds(fit_method,
      bmr = sign(top) * bmr, pars = unlist(modpars), conc, resp, onesidedp = .05,
      bmd = bmd, which.bound = "lower"
    )
    bmdu <- bmdbounds(fit_method,
      bmr = sign(top) * bmr, pars = unlist(modpars), conc, resp, onesidedp = .05,
      bmd = bmd, which.bound = "upper"
    )

    # apply bmd min
    if(!is.null(bmd_low_bnd) & !is.na(bmd)){
      # check if the argument is within its allowable range
      if (bmd_low_bnd > 0 & bmd_low_bnd <= 1) {
        # warning message for extreme values
        if (bmd_low_bnd < 1e-3){warning("The specified bmd_lower_bnd is less than 1e-3. This may result in an extremely low threshold value for BMD censoring. Suggested value is 0.1.")}
        min_conc <- min(conc[conc!=0])
        min_bmd <- min_conc*bmd_low_bnd
        if(bmd < min_bmd){
          bmd_diff <- min_bmd - bmd
          #shift all bmd to the right
          bmd <- bmd + bmd_diff
          bmdl <- bmdl + bmd_diff
          bmdu <- bmdu + bmd_diff
          }
      } else {
        warning("bmd_low_bnd must be between 0 and 1, not including 0.")
      }
    }

    # apply bmd max
    if(!is.null(bmd_up_bnd) & !is.na(bmd)){
      # check if the argument is within its allowable range
      if (bmd_up_bnd >= 1) {
        # warning message for extreme values
        if (bmd_up_bnd > 1e3) {warning("The specified bmd_up_bnd is larger than 1e3. This may result in an extremely high threshold value for BMD censoring. Suggested value is 10.")}
        max_conc <- max(conc)
        max_bmd <- max_conc*bmd_up_bnd
        if(bmd > max_bmd){
          #shift all bmd to the left
          bmd_diff <- bmd - max_bmd
          bmd <- bmd - bmd_diff
          bmdl <- bmdl - bmd_diff
          bmdu <- bmdu - bmd_diff
          }
      } else {
        warning("bmd_up_bnd must be greater than or equal to 1.")
      }
    }




  }

  top_over_cutoff <- abs(top) / cutoff
  conc <- paste(conc, collapse = "|")
  resp <- paste(resp, collapse = "|")

  # row contains the specified columns and any identifying, unused columns in the input
  name.list <- c(
    "n_gt_cutoff", "cutoff", "fit_method",
    "top_over_cutoff", "rmse", "a", "b", "tp", "p", "q", "ga", "la", "er", "bmr", "bmdl", "bmdu", "caikwt",
    "mll", "hitcall", "ac50", "ac50_loss", "top", "ac5", "ac10", "ac20", "acc", "ac1sd", "bmd", "conc", "resp", "errfun"
  )
  row <- as.data.frame(c(identifiers, mget(name.list)), stringsAsFactors = FALSE)
  return(row)
}
