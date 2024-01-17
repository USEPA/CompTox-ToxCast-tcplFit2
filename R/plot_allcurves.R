plot_allcurves <- function() {

  shortnames <- fitmodels[fitmodels != "cnst"]
  successes <- sapply(shortnames, function(x) {
    get(x)[["success"]]
  })
  if (do.plot && sum(successes, na.rm = T) == length(shortnames)) {
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

    ## Plot the Model Fits ##
    ggplot(data.frame(conc, resp), aes(x = log10(conc),y = resp))+
      geom_point(pch = 1,size = 2)+
      geom_line(data = estDR,
                aes(x = log10(X),y = value,colour = variable,lty = variable))+
      labs(colour = "Models",lty = "Models")+
      #scale_colour_manual(values = fit_cols)+
      xlab(expression(paste(log[10],"(Concentration) ",mu,"M")))+
      ylab("Responses")+
      theme_bw()

  }
}
