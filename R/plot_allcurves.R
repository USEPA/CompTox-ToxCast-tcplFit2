#' Plot All Fitted Curves From tcplfit2_core
#'
#' This function takes in output from `tcplfit2_core` and plots
#' the observed concentration response data with all fitted curves.
#' The function returns a ggplot object that is capable of
#' taking additional layers of ggplot add-ons.
#'
#' @param output Output from `concRespCore` or `tcplhit2_core`, containing information
#' about the winning curve fit of a compound.
#'
#' @param conc Vector of concentrations (not in log units).
#'
#' @param resp Vector of corresponding responses.
#'
#' @param log_conc Logical argument. If `TRUE`, convert the concentration (x-axis)
#' into log-10 scale. Defaults to `FALSE`.
#'
#' @return A ggplot plot of the observed concentration response data and
#' all the fitted curves.
#'
#' @export
#'

plot_allcurves <- function(output, conc, resp, log_conc = FALSE) {

  list2env(output, envir = environment())

  shortnames <- modelnames[modelnames != "cnst"]
  successes <- sapply(shortnames, function(x) {
    get(x)[["success"]]
  })

  # If any of the models from the output failed,
  # display the models that do fit and throw a warning message
  # about the failed models and they are not included in the plot.
  if (sum(successes, na.rm = T) != length(shortnames)) {
    failed_fits <- names(which(successes != 1))
    shortnames <- shortnames[! shortnames %in% failed_fits]
    warning(paste0("The following model(s) failed to fit: ", paste(failed_fits, collapse = ', '),
                  ", and will not be included in the plot."))
    }

    X <- seq(min(conc), max(conc),
             length.out = 100)
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

    if (!log_conc) {
      ## Plot the Model Fits ##
      p <- ggplot(data.frame(conc, resp), aes(x = conc,y = resp))+
        geom_point(pch = 1,size = 2)+
        geom_line(data = estDR,
                  aes(x = X,y = value,colour = variable,lty = variable))+
        labs(colour = "Models",lty = "Models")+
        xlab("Concentration")+
        ylab("Responses")+
        theme_bw()
      } else {
        p <- ggplot(data.frame(conc, resp), aes(x = log10(conc),y = resp))+
          geom_point(pch = 1,size = 2)+
          geom_line(data = estDR,
                    aes(x = log10(X),y = value,colour = variable,lty = variable))+
          labs(colour = "Models",lty = "Models")+
          xlab(expression(paste(log[10],"(Concentration) ",mu,"M")))+
          ylab("Responses")+
          theme_bw()
    }

  return(p)
}
