## do.plot option in tcplfit2_core

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

concRespPlot2 <- function(row,ymin=-120,ymax=120,draw.error.arrows=FALSE) {

  #every variable in row goes into the environment to make it easy
  #to update this function to use new row data.
  list2env(row,envir = environment())
  name <- unlist(name)
  assay <- unlist(assay)
  cutoff <- unlist(cutoff)
  bmr <- unlist(bmr)
  er <- unlist(er)
  fit_method <- unlist(fit_method)
  ac50 <- unlist(ac50)
  top <- unlist(top)
  bmd <- unlist(bmd)
  acc <- unlist(acc)
  hitcall <- unlist(hitcall)
  bmdl <- unlist(bmdl)
  bmdu <- unlist(bmdu)

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

  #gcalculate and plot model curves
  if(fit_method == "hill"){
    resp_plot <- do.call("hillfn",list(ps = unlist(modpars), x = conc_plot))
  } else if(!fit_method %in% c("cnst","none") ){
    resp_plot <- do.call(fit_method,list(ps = unlist(modpars), x = conc_plot))
  }

  basic <- ggplot(data.frame(logc, resp), aes(logc, resp)) +
    geom_point(pch = 1,size = 2)+ +
    geom_line(data.frame(conc_plot, resp_plot), aes(conc_plot, resp_plot)) +
    geom_rect(aes(xmin = min(logc),xmax = max(logc),ymin = -cutoff,ymax = cutoff),
              alpha = 0.15,fill = "skyblue")+
    # Titles and Labels
    xlab(expression(paste(log[10],"(Concentration) ",mu,"M")))+
    ylab(expression(paste(log[2],"(Fold Induction)")))
    # ggtitle(
    #   label = paste("Level 5 Best Model Fit",
    #                 mc4[which(mc4[,spid] == "01504209"),dsstox_substance_id],
    #                 sep = "\n")) +


    # Next, add the various potency layers.
    # BMD
  final_plot <- basic +
    geom_hline(
      aes(yintercept = bmr),
    ) +
    geom_segment(
      aes(x = log10(bmd), xend = log10(bmd), y = -0.5, yend = bmr)
    )

  return(final_plot)
}
