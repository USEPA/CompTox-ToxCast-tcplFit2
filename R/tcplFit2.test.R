#' Run a test of the tcploFit2 code
#'
#' This is jsut a test routine to show how tcplFit2 runs
#'
tcplFit2.test <- function() {

  conc <- list(.03,.1,.3,1,3,10,30,100)
  resp <- list(0,.2,.1,.4,.7,.9,.6, 1.2)
  row <- list(conc = conc, resp = resp, bmed = 0, cutoff = 1, onesd = .5)
  res <- concRespCore(row,fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
                                        "exp4", "exp5"), conthits = T, aicc = F, do.plot=T,verbose=T)
  return(res)
}

