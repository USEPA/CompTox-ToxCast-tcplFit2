test_that("concRespCore works", {
  data("signatures")
  row = list(conc=as.numeric(str_split(signatures[1,"conc"],"\\|")[[1]]),
             resp=as.numeric(str_split(signatures[1,"resp"],"\\|")[[1]]),
             bmed=0,
             cutoff=signatures[1,"cutoff"],
             onesd=signatures[1,"onesd"],
             name=signatures[1,"name"],
             assay=signatures[1,"signature"])
  out = concRespCore(row)

  expect_equal(out$fit_method, "exp4")
  expect_equal(out$hitcall, 0.99, tolerance = 1e-2)
  expect_equal(out$tp, 0.749, tolerance = 1e-2)
  expect_equal(out$ga, 9.59, tolerance = 1e-2)

})

test_that("HTPP global data internal check", {
  load("~/CompTox-ToxCast-tcplFit2/R/sysdata.rda")

  my_global <- NULL
  for (this.chem in unique(htpp_global_subset$chem_id)){
    ## for each chemical, prepare data for concRespCore
    temp <- htpp_global_input |> filter(chem_id == this.chem)
    ## prepare "row" list object per Description
    row <- list(
      pg_id = unique(temp$pg_id),
      stype = unique(temp$stype),
      chem_id = this.chem,
      min_conc = min(temp$conc),
      max_conc = max(temp$conc),
      n_conc = length(unique(temp$conc)),
      ctr_mean = CONTROL_GMAH$BMED,
      ctr_sd = CONTROL_GMAH$CUTOFF,
      conc = temp$conc,
      resp = temp$d,
      bmed = CONTROL_GMAH$BMED,
      cutoff = CONTROL_GMAH$CUTOFF,
      onesd = CONTROL_GMAH$ONESD,
      approach = "global",
      endpoint = "global"
    )
    ## run concRespCore on the ‘row’ list object
    newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2", "exp3","exp4", "exp5"), conthits = T, aicc = F, force.fit = FALSE, bidirectional = FALSE, AUC = FALSE))
    if(is.null(newLine) | class(newLine)=="try-error"){
      newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,  bidirectional = TRUE, fitmodels=c("cnst"), AUC = FALSE))
    }else if(class(newLine) == "try-error"){
      newLine <- row
    }else{
      rownames(newLine) <- ""
    }
    my_global <- rbind(my_global, newLine)
  }

  # compare results

})






