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

  skip_on_cran()

  # load necessary data
  load(here::here("R", "sysdata.rda"))

  ## Confirmed with Derik Haggard, for chemical C02,
  ## concentrations above 2.98 uM are removed from curve-fitting
  ## because higher concentrations cause cytotoxicity above the standard threshold of 50%.
  ## (i.e. concentration that is below the fitted EC50 value for a chemical is removed.)

  ## dropping higher concentrations for C02
  htpp_global_input = subset(htpp_global_input, htpp_global_input$trt != "C02" |
                               htpp_global_input$conc <= 2.98 )

  ## One can use the following code to verify that for curve fitting, C02 uses
  ## 5 unique concentrations, with the highest concentration being 2.98.
  ## htpp_global_subset[htpp_global_subset$trt == "C02",c("trt", "min_conc", "max_conc","n_conc", "conc")]
  ## unique(htpp_global_input[htpp_global_input$trt == "C02","conc"])


  my_global <- NULL
  for (this.chem in unique(htpp_global_subset$trt)){
    ## for each chemical, prepare data for concRespCore
    temp <- htpp_global_input[htpp_global_input$trt == this.chem, ]
    ## prepare "row" list object per Description
    row <- list(
      pg_id = unique(temp$pg_id),
      stype = unique(temp$stype),
      trt = this.chem,
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

  ## Compare results
  ## Compare BMD, BMDU, BMDL, top_over_cutoff, hit-call, top, and AC50

  ## vector operation, Order by chemical id to make sure we are comparing the same
  ## chemical
  my_global<- my_global[order(my_global$trt),]
  htpp_global_subset<- htpp_global_subset[order(htpp_global_subset$trt),]

  ## Differences in decimals places are normal rounding errors
  ## check if the differences in the BMD estimates exceed a threshold value
  ## BMD could be NA, replace NA with -1 so it wouldn't cause trouble with all().
  ## The maximum difference between BMD's is 0.0048828.
  my_global$bmd[is.na(my_global$bmd)] <- (-1)
  htpp_global_subset$bmd[is.na(htpp_global_subset$bmd)] <- (-1)
  ## Adjust the BMD values to 3 significant digits to be consistent with how they are being applied
  my_global[, c("bmd", "bmdu", "bmdl")] <- signif(my_global[, c("bmd", "bmdu", "bmdl")], 3)
  htpp_global_subset[, c("bmd", "bmdu", "bmdl")] <- signif(htpp_global_subset[, c("bmd", "bmdu", "bmdl")], 3)
  bmd_check <- all(abs(my_global$bmd - htpp_global_subset$bmd) < 1e-5)
  expect_true(bmd_check)

  ## Compare BMDU and BMDL with 3 significant digits
  ## BMDU and BMDL could be NA. Replacing NA with -1.
  my_global$bmdu[is.na(my_global$bmdu)] <- (-1)
  htpp_global_subset$bmdu[is.na(htpp_global_subset$bmdu)] <- (-1)
  bmdu_check <- all(abs(my_global$bmdu - htpp_global_subset$bmdu) < 1e-5)
  expect_true(bmdu_check)

  my_global$bmdl[is.na(my_global$bmdl)] <- (-1)
  htpp_global_subset$bmdl[is.na(htpp_global_subset$bmdl)] <- (-1)
  bmdl_check <- all(abs(my_global$bmdl - htpp_global_subset$bmdl) < 1e-5)
  expect_true(bmdl_check)

  ## check if the differences in hit-calls exceed a threshold value
  hitcall_check <- all(abs(my_global$hitcall - htpp_global_subset$hitcall) < 1e-5)
  expect_true(hitcall_check)

  ## check if the differences in top_over_cutoff exceed a threshold value
  ## top_over_cutoff is NA when the model is none - replace NA with -1
  my_global$top_over_cutoff[is.na(my_global$top_over_cutoff)] <- (-1)
  htpp_global_subset$top_over_cutoff[is.na(htpp_global_subset$top_over_cutoff)] <- (-1)
  top_cutoff_check <- all(abs(my_global$top_over_cutoff - htpp_global_subset$top_over_cutoff) < 1e-5)
  expect_true(hitcall_check)

  ## check if the differences in top exceed a threshold value
  my_global$top[is.na(my_global$top)] <- (-1)
  htpp_global_subset$top[is.na(htpp_global_subset$top)] <- (-1)
  top_check <- all(abs(my_global$top - htpp_global_subset$top) < 1e-5)
  expect_true(top_check)

  ## check if the differences in AC50 exceed a threshold value
  my_global$ac50[is.na(my_global$ac50)] <- (-1)
  htpp_global_subset$ac50[is.na(htpp_global_subset$ac50)] <- (-1)
  ac50_check <- all(abs(my_global$ac50 - htpp_global_subset$ac50) < 1e-5)
  expect_true(ac50_check)
})
