#' Concentration Response Plot - ggplot2
#'
#' This function takes output from `concRespCore` or `tcplhit2_core` to
#' generate a basic plot of the observed concentration-response data and the
#' best fit curve (winning model). A `ggplot` object, which users may customize
#' with additional `ggplot` layers, is returned.
#'
#' @param row Output from `concRespCore` or `tcplhit2_core`, containing information
#' about the best fit curve (winning model).
#' @param log_conc Logical argument. If `TRUE`, convert the concentrations (x-axis)
#' into log-10 scale. Defaults to `FALSE`.
#'
#' @return A `ggplot` object of the observed concentration-response data
#' overlaid with the best fit curve (winning model).
#' @import ggplot2
#'
#' @export
concRespPlot2 <- function(row, log_conc = FALSE) {

  # get the winning model curve
  fit_method <- row[,"fit_method"]

  #reformat conc and resp as vectors
  conc <- as.numeric(str_split(row[1,"conc"],"\\|")[[1]])
  resp <- as.numeric(str_split(row[1,"resp"],"\\|")[[1]])

  #get model parameters
  parnames = c("a", "tp", "b", "ga", "p", "la", "q")
  modpars = as.list(row[,parnames])
  modpars= modpars[!sapply(modpars, is.na)]

  # Generate a series of X for plotting the curve
  if (log_conc) {
    if (any(conc==0)) warning("Data contains untreated controls (conc = 0). A pseudo value replaces -Inf after log-transform.  The pseudo value is set to one log-unit below the lowest experimental `conc`.")
    conc <- log10(conc)
    # replace the negative infinity with a number that is one log-10 unit
    # less than the second lowest dose (in log).
    conc <- replace(conc, conc==-Inf, sort(conc)[2]-1)
    conc_plot <- seq(from=min(conc),to=max(conc),length=100)
    conc_plot <- replace(conc_plot, conc_plot==-Inf, sort(conc)[2]-1)
    # what will be passed into the object function to calculate responses for the curve
    calc.x <- 10**conc_plot
  } else {
    conc_plot <- seq(from=min(conc),to=max(conc),length=100)
    calc.x <- conc_plot
    }

  # calculate and plot model curves
  if(fit_method == "hill"){
    resp_plot <- do.call("hillfn",list(ps = unlist(modpars), x = calc.x))
  } else if(!fit_method %in% c("cnst","none") ){
    resp_plot <- do.call(fit_method,list(ps = unlist(modpars), x = calc.x))
  }

  basic <- ggplot(data.frame(conc, resp), aes(conc, resp)) +
    # Observed data points
    geom_point(pch = 1,size = 2) +
    # Winning Curve
    geom_line(data=data.frame(conc_plot, resp_plot), aes(x=conc_plot, y=resp_plot, color = fit_method)) +
    # Labels
    xlab("Concentration")+
    ylab("Response") +
    theme_bw()

  if (log_conc) basic <- basic + xlab(expression(paste(log[10],"(Concentration) ",mu,"M")))

  return(basic)
}
