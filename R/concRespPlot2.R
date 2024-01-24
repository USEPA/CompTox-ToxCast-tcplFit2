#' Concentration Response Plot - New
#'
#' This function takes in `concRespCore` or `tcplhit2_core` output and
#' generates a basic plot of the observed concentration response data with the
#' winning model curve. The function returns a ggplot object that is capable of
#' taking additional layers of ggplot add-ons.
#'
#' @param row Output from `concRespCore` or `tcplhit2_core`, containing information
#' about the winning curve fit of a compound.
#' @param log_conc Logical argument. If `TRUE`, convert the concentration (x-axis)
#' into log-10 scale. Defaults to `FALSE`.
#'
#' @return A ggplot object, a scatter plot of the concentration response data
#' overlaid with the winning model curve.
#'
#' @export
#'

concRespPlot2 <- function(row, log_conc = FALSE) {

  # get the winning model curve
  fit_method <- row[,"fit_method"]

  #reformat conc and resp as vectors
  conc <- as.numeric(str_split(row[1,"conc"],"\\|")[[1]])
  logc <- log10(conc)
  resp <- as.numeric(str_split(row[1,"resp"],"\\|")[[1]])

  #hard-code plotting points for curves
  logc_plot <- seq(from=min(logc),to=max(logc),length=100)
  conc_plot <- 10**logc_plot

  #get model parameters
  parnames = c("a", "tp", "b", "ga", "p", "la", "q")
  modpars = as.list(row[,parnames])
  modpars= modpars[!sapply(modpars, is.na)]

  #calculate and plot model curves
  if(fit_method == "hill"){
    resp_plot <- do.call("hillfn",list(ps = unlist(modpars), x = conc_plot))
  } else if(!fit_method %in% c("cnst","none") ){
    resp_plot <- do.call(fit_method,list(ps = unlist(modpars), x = conc_plot))
  }

  if (!log_conc) {
    basic <- ggplot(data.frame(conc, resp), aes(conc, resp)) +
      # Observed data points
      geom_point(pch = 1,size = 2) +
      # Winning Curve
      geom_line(data=data.frame(conc_plot, resp_plot), aes(x=conc_plot, y=resp_plot, color = fit_method)) +
      # Labels
      xlab("Concentration")+
      ylab("Response") +
      theme_bw()
  } else {
    basic <- ggplot(data.frame(logc, resp), aes(logc, resp)) +
      # Observed data points
      geom_point(pch = 1,size = 2) +
      # Winning Curve
      geom_line(data=data.frame(logc_plot, resp_plot), aes(x=logc_plot, y=resp_plot, color = fit_method)) +
      # Labels
      xlab(expression(paste(log[10],"(Concentration) ",mu,"M")))+
      ylab("Response") +
      theme_bw()

}
  return(basic)
}
