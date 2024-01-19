#' Plot All Fitted Curves
#'
#' This function takes in output from `tcplfit2_core` and plots
#' the observed concentration response data along with all the model fits.
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
#' @param log_conc Logical argument. If `TRUE`, convert the x-axis into log-10 scale.
#' Defaults to `FALSE`.
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
  if (sum(successes, na.rm = T) == length(shortnames)) {
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
  } else {
    stop("At least one of the model failed to fit.")
  }

  return(p)
}
