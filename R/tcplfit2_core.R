#' Concentration-response curve fitting
#'
#' Concentration response curve fitting using the methods from BMDExpress
#'
#' All models are equal to 0 at 0 concentration (zero background).
#' To add more models in the future, write a fit____ function, and add the model
#' name to the fitmodels and modelnames vectors.
#'
#' @param conc Vector of concentrations (NOT in log units).
#' @param resp Vector of responses.
#' @param cutoff Desired cutoff. If no absolute responses > cutoff and
#'   force.fit = F, will only fit constant model.
#' @param force.fit If force.fit = T, will fit all models regardless of cutoff.
#' @param bidirectional If bidirectional = F, will only give positive fits.
#' @param verbose If verbose = T, will print optimization details and aics.
#' @param do.plot If do.plot = T, will generate a plot comparing model curves.
#' @param fitmodels Vector of model names to try fitting. Missing models still
#'   return a skeleton output filled with NAs.
#' @param ... Other fitting parameters (deprecated).
#'
#' @import RColorBrewer
#' @import graphics
#' @importFrom stats median
#'
#' @return List of N(models) elements, one for each of the models run (up to 10),
#' followed by a last element "modelnames", which is a  vector of model names so
#' other functions can easily cycle through the output. For a full list, see the
#' documentation for the individual fitting method functions. For each model there
#' is a sublist with elements including:
#'   \itemize{
#'     \item success - was the model successfully fit
#'     \item aic - the AIC value
#'     \item cov - success of the the covariance matrix calculation
#'     \item rme - root mean error of the data around the curve
#'     \item modl - vector of model values at the civen concentrations
#'     \item tp - the top of the curve fit
#'     \item ga - the AC50 or Hill paramters
#'     \item er - the error term
#'     \item ... other paramters specific to the model (see the documentation for the specific models)
#'     \item tp_sd, ga_sd, p_sd, etc., the values of the standard deviations of the paramters for the models
#'     \item er_sd - standard deviation of the error term
#'     \item pars - the names of the parameters
#'     \item sds - the names of the standard deviations of the paramters
#'   }
#' @export
#'
#' @examples
#' conc <- c(.03, .1, .3, 1, 3, 10, 30, 100)
#' resp <- c(0, .1, 0, .2, .6, .9, 1.1, 1)
#' output <- tcplfit2_core(conc, resp, .8,
#'   fitmodels = c("cnst", "hill"), verbose = TRUE,
#'   do.plot = TRUE
#' )
tcplfit2_core <- function(conc, resp, cutoff, force.fit = FALSE, bidirectional = TRUE, verbose = FALSE, do.plot = FALSE,
                          fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3", "exp4", "exp5"),
                          ...) {
  logc <- log10(conc)
  rmds <- tapply(resp, logc, median)
  fitmodels <- unique(c("cnst", fitmodels)) # cnst models must be present for conthits but not chosen

  # first decide which of possible models will be fit
  modelnames <- c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3", "exp4", "exp5")
  # decide whether to run each model, then use generic functions to run model by name
  for (model in modelnames) {
    # only fit when four or more concentrations, the model is in fitmodels, and
    # ( either one response is above cutoff OR force.fit == T OR it's the constant model.)
    to.fit <- (length(rmds) >= 4 && model %in% fitmodels && (length(which(abs(rmds) >= cutoff)) > 0 || force.fit ||
      model == "cnst"))
    fname <- paste0("fit", model) # requires each model function have name "fit____" where ____ is the model name
    # use do.call to call fit function; cnst has different inputs than others.
    assign(model, do.call(fname, list(
        conc = conc, resp = resp, bidirectional = bidirectional, verbose = verbose,
        nofit = !to.fit
      )))
      if (to.fit) {
        if (model %in% c("poly1", "poly2", "pow", "exp2", "exp3")) {
          # methods that grow without bound: top defined as model value at max conc
          assign(model, append(get(model), list(top = get(model)$modl[which.max(abs(get(model)$modl))]))) # top is taken to be highest model value
          assign(model, append(get(model), list(ac50 = acy(.5 * get(model)$top, get(model), type = model))))
        } else if (model %in% c("hill", "exp4", "exp5")) {
          # methods with a theoretical top/ac50
          assign(model, append(get(model), list(top = get(model)$tp)))
          assign(model, append(get(model), list(ac50 = get(model)$ga)))
        } else if (model == "gnls") {
          # gnls methods; use calculated top/ac50, etc.
          assign(model, append(get(model), list(top = acy(0, get(model), type = model, returntop = T))))
          assign(model, append(get(model), list(ac50 = acy(.5 * get(model)$top, get(model), type = model))))
          assign(model, append(get(model), list(ac50_loss = acy(.5 * get(model)$top, get(model), type = model, getloss = T))))
        }
      }

  }
  # optionally print out AICs
  if (verbose) {
    print("aic values:")
    aics <- sapply(modelnames, function(x) {
      get(x)[["aic"]]
    })
    names(aics) <- modelnames
    print(aics)
    cat("Winner: ", modelnames[which.min(aics)])
  }

  # optionallky plot all models if there's at least one model to plot
  shortnames <- modelnames[modelnames != "cnst"]
  successes <- sapply(shortnames, function(x) {
    get(x)[["success"]]
  })
  if (do.plot && sum(successes, na.rm = T) == length(shortnames)) {
    resp <- resp[order(logc)]
    par(xpd = T)
    cols <- c("black", brewer.pal(9, "Set1"))
    n <- length(logc)
    allresp <- c(resp, sapply(shortnames, function(x) {
      get(x)[["modl"]][order(logc)]
    }))
    logc <- logc[order(logc)]
    plot(rep(logc, length.out = length(allresp)), allresp, col = rep(cols, each = n), pch = 16)

    for (i in 1:length(allresp)) {
      points(logc, allresp[((i - 1) * n + 1):(i * n)], col = cols[i], type = "l")
    }

    legend("top", legend = c("resp", shortnames), col = cols, pch = 16, ncol = 10, inset = c(0, -.1))
  }

  # put all the model outputs into one list and return
  out <- c(
    mget(modelnames),
    list(modelnames = modelnames, ...)
  )

  return(out)
}
