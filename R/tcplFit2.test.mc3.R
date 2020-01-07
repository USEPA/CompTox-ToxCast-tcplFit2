#' Run a test of the tcploFit2 code
#'
#' This function shows how to run multiple chemicals using data from invitrodb
#'
#' @importFrom stats sd
#' @importFrom stats mad
tcplFit2.test.mc3 <- function() {

  # read in the data
  file <- "data/mc3.RData"
  load(file=file)
  print(dim(mc3))

  # set up a 3 x 2 grid for the plots
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  # determine the background variation
  temp <- mc3[mc3$logc<= -2,"resp"]
  bmad <- mad(temp)
  onesd <- sd(temp)
  cutoff <- 3*bmad

  # select six samples. Note that there may be more than one sample processed for a given chemical
  spid.list <- unique(mc3$spid)
  spid.list <- spid.list[1:6]

  for(spid in spid.list) {
    # select the data for just this sample
    temp <- mc3[is.element(mc3$spid,spid),]

    # The data file has stored concentration in log10 form, so fix that
    conc <- 10**temp$logc
    resp <- temp$resp

    # pull out all of the chemical identifiers and the name of the assay
    dtxsid <- temp[1,"dtxsid"]
    casrn <- temp[1,"casrn"]
    name <- temp[1,"name"]
    assay <- temp[1,"assay"]

    # create the row object
    row <- list(conc = conc, resp = resp, bmed = 0, cutoff = cutoff, onesd = onesd,assay=assay,dtxsid=dtxsid,casrn=casrn,name=name)

    # run the concentration-response modeling for a single sample
    res <- concRespCore(row,fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
                                          "exp4", "exp5"),conthits = T, aicc = F,bidirectional=F)

    # plot the results
    concRespPlot(res$summary,ymin=-10,ymax=100)
  }
  browser()
}

