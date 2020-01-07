#' Run a test of the tcploFit2 code
#'
#' This function shows how to run multiple chemicals using data from invitrodb
#'
tcplFit2.test.mc3 <- function() {

  file <- "data/mc3_TOX21_ERa_BLA_Agonist_ratio.RData"
  load(file=file)
  print(dim(mc3))

  par(mfrow=c(3,2),mar=c(4,4,2,2))

  # determine the background variation
  temp <- mc3[mc3$logc<= -2,"resp"]
  bmad <- mad(temp)
  onesd <- sd(temp)
  cutoff <- 3*bmad
  spid.list <- unique(mc3$spid)
  spid.list <- spid.list[1:6]
  for(spid in spid.list) {
    temp <- mc3[is.element(mc3$spid,spid),]
    conc <- 10**temp$logc
    resp <- temp$resp
    dtxsid <- temp[1,"dtxsid"]
    casrn <- temp[1,"casrn"]
    name <- temp[1,"name"]
    assay <- temp[1,"assay"]
    row <- list(conc = conc, resp = resp, bmed = 0, cutoff = cutoff, onesd = onesd,assay=assay,dtxsid=dtxsid,casrn=casrn,name=name)
    res <- concRespCore(row,fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
                                          "exp4", "exp5"),conthits = T, aicc = F,bidirectional=F)
    concRespPlot(res,ymin=-10,ymax=100)
  }
  browser()
}

