#' Hitcalling Function
#'
#' Core of hitcalling function. This method chooses the winning model from tcplfit2_core,
#' extracts the top and ac50, computes the hitcall, and calculates bmd/bmdl/bmdu among other
#' statistics. Nested model selection is used to choose between poly1/poly2, then
#' the model with the lowest AIC (or AICc) is declared the winner. Continuous
#' hitcalls requires tcplfit2_core to be run with force.fit = T and "cnst" never to
#' be chosen as the winner.
#'
#' @param params The output from tcplfit2_core
#' @param conc list of concentrations (not in log units)
#' @param resp list of corresponding responses
#' @param bmed median of noise estimate. Default 0
#' @param cutoff noise cutoff
#' @param onesd 1 standard deviation of the noise (for bmd calculation)
#' @param conthits conthits = TUE uses continuous hitcalls, otherwise they're
#'   discrete. Default TRUE
#' @param aicc aicc = TRUE uses corrected AIC to choose winning method; otherwise
#'   regular AIC. Default FALSE
#' @param identifiers A one-row data frame containing identifiers of the concentration-response profile,
#'   such as the chemical name or other identifiers, and any assay identifiers. The column names
#'   identify the type of value. This can be NULL. The values will be included in the output
#'   summary data frame
#'
#' @return A list of with the detailed results from all of the different model fits.
#' The elements of summary are:
#'   \itemize{
#'     \item any elements of the identifiers input
#'     \item n_gt_cutoff - number of data points above the cutoff
#'     \item cutoff - noise cutoff
#'     \item fit_method - curve fit method
#'     \item top_over_cutoff - top divided by cutoff
#'     \item rmse - RMSE of the data points arount the best model curve
#'     \item a - fitting parameter methods: exp2, exp3, poly1, poly2, pow
#'     \item b - fitting parameter methods: exp2, exp3, ploy2
#'     \item p - fitting parameter methods: exp3, exp5, gnls, hill, pow
#'     \item q - fitting parameter methods: gnls,
#'     \item tp - top of the curve
#'     \item ga - ac50 for the rising curve in a gnls model or the Hill model
#'     \item la - ac50 for the falling curve in a gnls model
#'     \item er - fitted error term for plotting error bars
#'     \item bmr - benchmark response; level at which bmd is calculated = onesd*1.349
#'     \item bmd - benchmark dose, curve value at bmr
#'     \item bmdl - lower limit on the bmd
#'     \item bmdu - upper limit on the bmd
#'     \item caikwt - one factor used in calculating the continuous hitcall. It is calcalated from the formula
#'       = exp(-aic(cnst)/2)/(exp(-aic(cnst)/2) + exp(-aic(fit_method)/2)) and measures how much lower the
#'       selected method AIC is than that for the constant model
#'     \item mll - anoter factor used in calcualting the continuous hitcall = length(modpars) - aic(fit_method)/2
#'     \item hitcall - the final hitcall, a value ranging from 0 to 1
#'     \item top - curve top
#'     \item ac50 - curve value at 50\% of top, curve value at cutoff
#'     \item lc50 - curve value at 50\% of top corresponding to the loss side of the gain-loss curve
#'     \item ac5 - curve value at 5\% of top
#'     \item ac10 - curve value at 10\% of top
#'     \item ac20 - curve value at 20\% of top
#'     \item acc - curve value at 1 standard deviation
#'     \item conc - conc string separated by |'s
#'     \item resp - response string separated by |'s
#'   }
#' @export
#'
tcplhit2_core <- function(params, conc, resp, cutoff, onesd, bmed = 0, conthits = T, aicc = F, identifiers = NULL) {
  # initialize parameters to NA
  a <- b <- tp <- p <- q <- ga <- la <- er <- top <- ac50 <- ac50_loss <- ac5 <- ac10 <- ac20 <- acc <- ac1sd <- bmd <- NA_real_
  bmdl <- bmdu <- caikwt <- mll <- NA_real_
  # get aics and degrees of freedom
  aics <- sapply(params$modelnames, function(x) {
    params[[x]][["aic"]]
  })
  dfs <- sapply(params$modelnames, function(x) {
    if(x == "cnst") 1L
    else length(params[[x]][gsub("_sd","",names(params[[x]])[grepl("_sd",names(params[[x]]))])])
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
    } else {
      fit_method <- names(saics)[which.min(saics)]
    }
    fitout <- params[[fit_method]]
    rmse <- fitout$rme
    # hacky way to get modpars without hardcoding the names or needing the list
    # basically each model (except cnst) has an sd for each parameter
    # we can use that to select all model parameters
    modpars <- fitout[gsub("_sd","",names(fitout)[grepl("_sd",names(fitout))])]
    if(fit_method == "cnst") modpars <- fitout["er"]  #since no sd for cnst
    list2env(fitout, envir = environment()) # put all parameters in environment
  }

  n_gt_cutoff <- sum(abs(resp) > cutoff)

  # compute discrete or continuous hitcalls
  if (fit_method == "none") {
    hitcall <- 0
  } else if (conthits) {
    mll <- length(modpars) - aics[[fit_method]] / 2
    hitcall <- hitcontinner(conc, resp, top, cutoff, er,
      ps = unlist(modpars), fit_method,
      caikwt = caikwt, mll = mll
    )
  } else {
    hitcall <- hitloginner(conc, resp, top, cutoff, ac50)
  }

  bmr <- onesd * 1.349 # magic bmr is hard-coded
  if (hitcall > 0) {

    # fill ac's; can put after hit logic
    ac5 <- acy(.05 * top, modpars, type = fit_method) # note: cnst model automatically returns NAs
    ac10 <- acy(.1 * top, modpars, type = fit_method)
    ac20 <- acy(.2 * top, modpars, type = fit_method)
    acc <- acy(sign(top) * cutoff, modpars, type = fit_method)
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
  }

  top_over_cutoff <- abs(top) / cutoff
  conc <- paste(conc, collapse = "|")
  resp <- paste(resp, collapse = "|")

  # row contains the specified columns and any identifying, unused columns in the input
  name.list <- c(
    "n_gt_cutoff", "cutoff", "fit_method",
    "top_over_cutoff", "rmse", "a", "b", "tp", "p", "q", "ga", "la", "er", "bmr", "bmdl", "bmdu", "caikwt",
    "mll", "hitcall", "ac50", "ac50_loss", "top", "ac5", "ac10", "ac20", "acc", "ac1sd", "bmd", "conc", "resp"
  )
  row <- as.data.frame(c(identifiers, mget(name.list)), stringsAsFactors = F)
  return(row)
}
