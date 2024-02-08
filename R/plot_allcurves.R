#' Plot All Curves Fit with tcplfit2_core - ggplot2
#'
#' This function takes output from `tcplfit2_core` and generates a basic plot of
#' the observed concentration-response data with all resulting curve fits.
#' A `ggplot` object, which users may customize with additional `ggplot` layers,
#' is returned.
#'
#' @param modelfits Output from `tcplfit2_core`, contains resulting fits for all
#' models used to evaluate the observed concentration-response data.
#'
#' @param conc Vector of concentrations (NOT in log units).
#'
#' @param resp Vector of responses.
#'
#' @param log_conc Logical argument. If `TRUE`, convert the concentrations (x-axis)
#' into log-10 scale. Defaults to `FALSE`.
#'
#' @return A `ggplot` object of the observed concentration-response data and
#' all resulting curve fits from `tcplfit2_core`. (Note: The constant model is
#' not included, and only the successful fits will be displayed.)
#'
#' @export
plot_allcurves <- function(modelfits, conc, resp, log_conc = FALSE) {

  list2env(modelfits, envir = environment())

  shortnames <- modelnames[modelnames != "cnst"]
  successes <- sapply(shortnames, function(x) {
    get(x)[["success"]]
  })

  # If any of the models from the output failed,
  # display the models that do fit and throw a warning message
  # about the failed models being excluded from the plot.
  if (sum(successes, na.rm = T) != length(shortnames)) {
    failed_fits <- names(which(successes != 1))
    shortnames <- shortnames[! shortnames %in% failed_fits]
    warning(paste0("The following model(s) failed to fit: ", paste(failed_fits, collapse = ', '),
                  ", and are excluded from the plot."))
  }

  if (log_conc) {
    conc <- log10(conc)
  }

  X <- seq(min(conc), max(conc), length.out = 100)
  allresp <- NULL
  pars_lists <- sapply(shortnames, function(x) {get(x)[["pars"]]})

  for (i in 1:length(pars_lists)) {
    model <- names(pars_lists)[i]
    values <- sapply(pars_lists[[model]], function(x) {
      get(model)[[x]]})
    if (model == "hill") model <- paste0(model,"fn")
    y <- do.call(model, list(values, X))
    allresp <- cbind(allresp, y)
  }

  colnames(allresp) <- shortnames

  # Format data into a data.frame for ease of plotting.
  estDR <- cbind.data.frame(X,allresp) %>%
    reshape2::melt(data = .,measure.vars = shortnames)


  ## Plot the Model Fits ##
  p <- ggplot(data.frame(conc, resp), aes(x = conc,y = resp)) +
    geom_point(pch = 1,size = 2)+
    geom_line(data = estDR,
              aes(x = X,y = value,colour = variable,lty = variable))+
    labs(colour = "Models",lty = "Models")+
    xlab("Concentration")+
    ylab("Responses")+
    theme_bw()

  if (log_conc) p <- p + xlab(expression(paste(log[10],"(Concentration) ",mu,"M")))

  return(p)
}
