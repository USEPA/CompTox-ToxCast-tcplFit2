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

  ## Confirmed with Derik, concentrations above 2.98 uM can be remove for C02
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

  ## vector operation, order by trt to make sure we are comparing the appropriate output
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


test_that("HTPP category data internal check", {

  skip_on_cran()

  ## load necessary data
  ## Code below is commented out, not needed when running all tests in the package at once (such as with testthat::test_local().)
  ## Do need to un-comment and run this code to load the data if one is running this test interactively in the console.

  #load(here::here("R", "sysdata.rda"))

  ## Confirmed with Derik, concentrations above 2.98 uM can be remove for C02
  ## because higher concentrations cause cytotoxicity above the standard threshold of 50%.
  ## (i.e. concentration that is below the fitted EC50 value for a chemical is removed.)

  ## dropping higher concentrations for C02
  htpp_cat_input = subset(htpp_cat_input, htpp_cat_input$trt !=  "C02" |
                            htpp_cat_input$conc <= 2.98 )

  ## One can use the following code to verify that for curve fitting, C02 uses
  ## 5 unique concentrations, with the highest concentration being 2.98.
  ## htpp_cat_subset[htpp_cat_subset$trt == "C02",c("trt", "min_conc", "max_conc","n_conc", "conc")]
  ## unique(htpp_cat_input[htpp_cat_input$trt == "C02","conc"])

  my_category <- NULL
  for (this.chem in unique(htpp_cat_subset$trt)){
    temp <- htpp_cat_input[htpp_cat_input$trt == this.chem, ]
    for (this.cat in unique(temp$category_name_r)) {
      this.sub <- temp[temp$category_name_r == this.cat, ]
      ## prepare "row" list object per Description
      Metadata <- CONTROL_CMAH[CONTROL_CMAH$category_name_r == this.cat, ]
      row <- list(
        pg_id = unique(this.sub$pg_id),
        stype = unique(this.sub$stype),
        trt = this.chem,
        min_conc = min(this.sub$conc),
        max_conc = max(this.sub$conc),
        n_conc = length(unique(this.sub$conc)),
        ctr_mean = Metadata$BMED,
        ctr_sd = Metadata$CUTOFF,
        conc = this.sub$conc,
        resp = this.sub$d,
        bmed = Metadata$BMED,
        cutoff = Metadata$CUTOFF,
        onesd = Metadata$ONESD,
        approach = "category",
        endpoint = this.cat
      )
      ## run concRespCore on the ‘row’ list object
      newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2",
                                                     "exp3","exp4", "exp5"), conthits = T, aicc = F,
                                  force.fit = FALSE, bidirectional = FALSE, AUC = FALSE))
      if(is.null(newLine) | class(newLine)=="try-error"){
        newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,
                                    bidirectional = TRUE, fitmodels=c("cnst"), AUC = FALSE))
      }else if(class(newLine) == "try-error"){
        newLine <- row
      }else{
        rownames(newLine) <- ""
      }

      my_category <- rbind(my_category, newLine)
    }
  }

  ## Compare results
  ## Compare BMD, BMDU, BMDL, top_over_cutoff, hit-call, top, and AC50

  ## Compare by vector operation, order by trt to make sure we are comparing the appropriate output
  my_category<- my_category[order(my_category$trt),]
  htpp_cat_subset<- htpp_cat_subset[order(htpp_cat_subset$trt),]

  ## Differences in decimals places are normal rounding errors
  ## check if the differences in the BMD estimates exceed a threshold value
  ## BMD could be NA, replace NA with -1 so it wouldn't cause trouble with all()
  my_category$bmd[is.na(my_category$bmd)] <- (-1)
  htpp_cat_subset$bmd[is.na(htpp_cat_subset$bmd)] <- (-1)
  ## Adjust the BMD values to 3 significant digits to be consistent with how they are being applied
  my_category[, c("bmd", "bmdu", "bmdl")] <- signif(my_category[, c("bmd", "bmdu", "bmdl")], 3)
  htpp_cat_subset[, c("bmd", "bmdu", "bmdl")] <- signif(htpp_cat_subset[, c("bmd", "bmdu", "bmdl")], 3)
  bmd_check <- all(abs(my_category$bmd - htpp_cat_subset$bmd) < 1e-5)
  expect_true(bmd_check)

  ## Compare BMDU and BMDL with 3 significant digits
  ## BMDU and BMDL could be NA. Replacing NA with -1.
  my_category$bmdu[is.na(my_category$bmdu)] <- (-1)
  htpp_cat_subset$bmdu[is.na(htpp_cat_subset$bmdu)] <- (-1)
  bmdu_check <- all(abs(my_category$bmdu - htpp_cat_subset$bmdu) < 1e-5)
  expect_true(bmdu_check)

  my_category$bmdl[is.na(my_category$bmdl)] <- (-1)
  htpp_cat_subset$bmdl[is.na(htpp_cat_subset$bmdl)] <- (-1)
  bmdl_check <- all(abs(my_category$bmdl - htpp_cat_subset$bmdl) < 1e-5)
  expect_true(bmdl_check)

  # check if the differences in hit-calls exceed a threshold value
  hitcall_check <- all(abs(my_category$hitcall - htpp_cat_subset$hitcall) < 1e-5)
  expect_true(hitcall_check)

  # check if the differences in top_over_cutoff exceed a threshold value
  ## top_over_cutoff is NA when the model is none - replace NA with -1
  my_category$top_over_cutoff[is.na(my_category$top_over_cutoff)] <- (-1)
  htpp_cat_subset$top_over_cutoff[is.na(htpp_cat_subset$top_over_cutoff)] <- (-1)
  top_cutoff_check <- all(abs(my_category$top_over_cutoff - htpp_cat_subset$top_over_cutoff) < 1e-5)
  expect_true(hitcall_check)

  ## check if the differences in top exceed a threshold value
  my_category$top[is.na(my_category$top)] <- (-1)
  htpp_cat_subset$top[is.na(htpp_cat_subset$top)] <- (-1)
  top_check <- all(abs(my_category$top - htpp_cat_subset$top) < 1e-5)
  expect_true(top_check)

  ## check if the differences in AC50 exceed a threshold value
  my_category$ac50[is.na(my_category$ac50)] <- (-1)
  htpp_cat_subset$ac50[is.na(htpp_cat_subset$ac50)] <- (-1)
  ac50_check <- all(abs(my_category$ac50 - htpp_cat_subset$ac50) < 1e-5)
  expect_true(ac50_check)

})

test_that("HTPP feature data internal check", {

  skip_on_cran()

  ## load necessary data
  ## Code below is commented out, not needed when running all tests in the package at once (such as with testthat::test_local().)
  ## Do need to un-comment and run this code to load the data if one is running this test interactively in the console.

  #load(here::here("R", "sysdata.rda"))

  ## Confirmed with Derik, concentrations above 2.98 uM can be remove for C02
  ## and concentrations above 10 uM can be removed for A05 because with these two samples
  ## higher concentrations cause cytotoxicity above the standard threshold of 50%.
  ## (i.e. concentration that is below the fitted EC50 value for a chemical is removed.)

  chem1 <- subset(htpp_feature_input, htpp_feature_input$trt == "C02" &
                    htpp_feature_input$conc <= 2.98)
  chem2 <- subset(htpp_feature_input, htpp_feature_input$trt == "A05" &
                    htpp_feature_input$conc <= 10)
  htpp_feature_input = subset(htpp_feature_input, htpp_feature_input$trt !=  "C02" &
                                htpp_feature_input$trt !=  "A05" )
  htpp_feature_input <- rbind(htpp_feature_input, chem1, chem2)

  my_feature <- NULL
  for (this.chem in unique(htpp_feature_subset$trt)){
    temp <- htpp_feature_input[htpp_feature_input$trt == this.chem, ]
    for (this.fname in unique(temp$Feature)) {
      this.sub <- temp[temp$Feature == this.fname, ]
      ## prepare "row" list object per Description
      Metadata <- CONTROL_FMAH[CONTROL_FMAH$fname == this.fname, ]
      row <- list(
        pg_id = unique(this.sub$pg_id),
        stype = unique(this.sub$stype),
        trt = this.chem,
        min_conc = min(this.sub$conc),
        max_conc = max(this.sub$conc),
        n_conc = length(unique(this.sub$conc)),
        ctr_mean = Metadata$BMED,
        ctr_sd = Metadata$CUTOFF,
        conc = this.sub$conc,
        resp = this.sub$d,
        bmed = Metadata$BMED,
        cutoff = Metadata$CUTOFF,
        onesd = Metadata$ONESD,
        approach = "feature",
        endpoint = this.fname
      )
      ## run concRespCore on the ‘row’ list object
      newLine <- try(concRespCore(row, fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2",
                                                     "exp3","exp4", "exp5"), conthits = T, aicc = F,
                                  force.fit = FALSE, bidirectional = TRUE, AUC = FALSE))
      if(is.null(newLine) | class(newLine)=="try-error"){
        newLine <- try(concRespCore(row, conthits = F, aicc = FALSE, force.fit = FALSE,
                                    bidirectional = TRUE, fitmodels=c("cnst"), AUC = FALSE))
      }else if(class(newLine) == "try-error"){
        newLine <- row
      }else{
        rownames(newLine) <- ""
      }

      my_feature <- rbind(my_feature, newLine)
    }
  }

  ## Compare results
  ## Compare BMD, BMDU, BMDL, top_over_cutoff, hit-call, top, and AC50

  ## Compare by vector operation, order by trt to make sure we are comparing the appropriate output
  my_feature <- my_feature[order(my_feature$trt),]
  htpp_feature_subset<- htpp_feature_subset[order(htpp_feature_subset$trt),]

  ## Differences in decimals places are normal rounding errors
  ## check if the differences in the BMD estimates exceed a threshold value
  ## BMD could be NA, replace NA with -1 so it wouldn't cause trouble with all()
  my_feature$bmd[is.na(my_feature$bmd)] <- (-1)
  htpp_feature_subset$bmd[is.na(htpp_feature_subset$bmd)] <- (-1)

  ## Adjust the BMD values to 3 significant digits to be consistent with how they are being applied
  my_feature[, c("bmd", "bmdu", "bmdl")] <- signif(my_feature[, c("bmd", "bmdu", "bmdl")], 3)
  htpp_feature_subset[, c("bmd", "bmdu", "bmdl")] <- signif(htpp_feature_subset[, c("bmd", "bmdu", "bmdl")], 3)
  bmd_check <- all(abs(my_feature$bmd - htpp_feature_subset$bmd) < 1e-5)
  expect_true(bmd_check)

  ## Compare BMDU and BMDL with 3 significant digits
  ## BMDU and BMDL could be NA. Replacing NA with -1.
  my_feature$bmdl[is.na(my_feature$bmdl)] <- (-1)
  htpp_feature_subset$bmdl[is.na(htpp_feature_subset$bmdl)] <- (-1)
  bmdl_check <- all(abs(my_feature$bmdl - htpp_feature_subset$bmdl) < 1e-5)
  expect_true(bmdl_check)

  my_feature$bmdu[is.na(my_feature$bmdu)] <- (-1)
  htpp_feature_subset$bmdu[is.na(htpp_feature_subset$bmdu)] <- (-1)
  bmdu_check <- all(abs(my_feature$bmdu - htpp_feature_subset$bmdu) < 1e-5)
  expect_true(bmdu_check)

  # check if the differences in hit-calls exceed a threshold value
  hitcall_check <- all(abs(my_feature$hitcall - htpp_feature_subset$hitcall) < 1e-5)
  expect_true(hitcall_check)

  # check if the differences in top_over_cutoff exceed a threshold value
  my_feature$top_over_cutoff[is.na(my_feature$top_over_cutoff)] <- (-1)
  htpp_feature_subset$top_over_cutoff[is.na(htpp_feature_subset$top_over_cutoff)] <- -1
  top_cutoff_check <- all(abs(my_feature$top_over_cutoff - htpp_feature_subset$top_over_cutoff) < 1e-5)
  expect_true(hitcall_check)

  ## check if the differences in top exceed a threshold value
  my_feature$top[is.na(my_feature$top)] <- (-1)
  htpp_feature_subset$top[is.na(htpp_feature_subset$top)] <- (-1)
  top_check <- all(abs(my_feature$top - htpp_feature_subset$top) < 1e-5)
  expect_true(top_check)

  ## check if the differences in AC50 exceed a threshold value
  my_feature$ac50[is.na(my_feature$ac50)] <- (-1)
  htpp_feature_subset$ac50[is.na(htpp_feature_subset$ac50)] <- (-1)
  ac50_check <- all(abs(my_feature$ac50 - htpp_feature_subset$ac50) < 1e-5)
  expect_true(ac50_check)
})


test_that("HTTr signature data internal check", {

  skip_on_cran()

  ## load necessary data
  ## Code below is commented out, not needed when running all tests in the package at once (such as with testthat::test_local().)
  ## Do need to un-comment and run this code to load the data if one is running this test interactively in the console.

  #load(here::here("R", "sysdata.rda"))

  # turn signature_input into a list of rows for lapply to use
  signature_input = as.list(as.data.frame(t(signature_input), stringsAsFactors = F))

  # run curve-fitting - this is taking a little more than 1 minute to run
  my_signature = lapply(X=signature_input, FUN =concRespCore,
                        fitmodels = c("cnst", "hill", "poly1", "poly2", "pow", "exp2",
                                      "exp3", "exp4", "exp5"),
                        bmr_scale=1.349, aicc=F, conthits=T, bmd_low_bnd=0.1,
                        bmd_up_bnd=10, verbose=F)
  # turn the result lists into a data frame
  my_signature = as.data.frame(data.table::rbindlist(my_signature))

  ## Compare BMD, BMDU, BMDL, top_over_cutoff, hit-call, top, and AC50

  ## Compare by vector operation, order by trt to make sure we are comparing the appropriate output
  my_signature <- my_signature[order(my_signature$trt, my_signature$signature),]
  signature_sub <- signature_sub[order(signature_sub$trt, signature_sub$signature),]

  ## Differences in decimals places are normal rounding errors
  ## check if the differences in the BMD estimates exceed a threshold value
  ## BMD could be NA, replace NA with -1 so it wouldn't cause trouble with all()
  my_signature$bmd[is.na(my_signature$bmd)] <- (-1)
  signature_sub$bmd[is.na(signature_sub$bmd)] <- (-1)

  ## Adjust the BMD values to 3 significant digits to be consistent with how they are being applied
  my_signature[, c("bmd", "bmdu", "bmdl")] <- signif(my_signature[, c("bmd", "bmdu", "bmdl")], 3)
  signature_sub[, c("bmd", "bmdu", "bmdl")] <- signif(signature_sub[, c("bmd", "bmdu", "bmdl")], 3)
  bmd_check <- all(abs(my_signature$bmd - signature_sub$bmd) < 1e-5)
  expect_true(bmd_check)

  ## Compare BMDU and BMDL with 3 significant digits
  ## BMDU and BMDL could be NA. Replacing NA with -1.
  my_signature$bmdl[is.na(my_signature$bmdl)] <- (-1)
  signature_sub$bmdl[is.na(signature_sub$bmdl)] <- (-1)
  bmdl_check <- all(abs(my_signature$bmdl - signature_sub$bmdl) < 1e-5)
  expect_true(bmdl_check)

  my_signature$bmdu[is.na(my_signature$bmdu)] <- (-1)
  signature_sub$bmdu[is.na(signature_sub$bmdu)] <- (-1)
  bmdu_check <- all(abs(my_signature$bmdu - signature_sub$bmdu) < 1e-5)
  expect_true(bmdu_check)

  # check if the differences in hit-calls exceed a threshold value
  hitcall_check <- all(abs(my_signature$hitcall - signature_sub$hitcall) < 1e-5)
  expect_true(hitcall_check)

  # check if the differences in top_over_cutoff exceed a threshold value
  my_signature$top_over_cutoff[is.na(my_signature$top_over_cutoff)] <- (-1)
  signature_sub$top_over_cutoff[is.na(signature_sub$top_over_cutoff)] <- (-1)
  top_cutoff_check <- all(abs(my_signature$top_over_cutoff - signature_sub$top_over_cutoff) < 1e-5)
  expect_true(hitcall_check)

  ## check if the differences in top exceed a threshold value
  my_signature$top[is.na(my_signature$top)] <- (-1)
  signature_sub$top[is.na(signature_sub$top)] <- (-1)
  top_check <- all(abs(my_signature$top - signature_sub$top) < 1e-5)
  expect_true(top_check)

  ## check if the differences in AC50 exceed a threshold value
  my_signature$ac50[is.na(my_signature$ac50)] <- (-1)
  signature_sub$ac50[is.na(signature_sub$ac50)] <- (-1)
  ac50_check <- all(abs(my_signature$ac50 - signature_sub$ac50) < 1e-5)
  expect_true(ac50_check)

})

test_that("HTTr gene data internal check", {

  skip_on_cran()

  ## load necessary data
  ## Code below is commented out, not needed when running all tests in the package at once (such as with testthat::test_local().)
  ## Do need to un-comment and run this code to load the data if one is running this test interactively in the console.

  #load(here::here("R", "sysdata.rda"))

  # turn gene input into a list of rows for lapply to use
  gene_input = as.list(as.data.frame(t(gene_input), stringsAsFactors = F))

  # run curve-fitting
  my_gene = lapply(X=gene_input, FUN = concRespCore,
                   fitmodels = c("cnst", "hill", "poly1", "poly2", "pow",
                                 "exp2", "exp3", "exp4", "exp5"), aicc = F)

  # turn the result lists into a data frame
  my_gene = as.data.frame(data.table::rbindlist(my_gene))

  ## Compare BMD, BMDU, BMDL, top_over_cutoff, hit-call, top, and AC50

  ## Compare by vector operation, order by trt to make sure we are comparing the appropriate output
  my_gene <- my_gene[order(my_gene$trt, my_gene$gene),]
  gene_sub <- gene_sub[order(gene_sub$trt, gene_sub$gene),]

  ## Differences in decimals places are normal rounding errors
  ## check if the differences in the BMD estimates exceed a threshold value
  ## BMD could be NA, replace NA with -1 so it wouldn't cause trouble with all()
  my_gene$bmd[is.na(my_gene$bmd)] <- (-1)
  gene_sub$bmd[is.na(gene_sub$bmd)] <- (-1)

  ## Adjust the BMD values to 3 significant digits to be consistent with how they are being applied
  my_gene[, c("bmd", "bmdu", "bmdl")] <- signif(my_gene[, c("bmd", "bmdu", "bmdl")], 3)
  gene_sub[, c("bmd", "bmdu", "bmdl")] <- signif(gene_sub[, c("bmd", "bmdu", "bmdl")], 3)
  bmd_check <- all(abs(my_gene$bmd - gene_sub$bmd) < 1e-5)
  expect_true(bmd_check)

  ## Compare BMDU and BMDL with 3 significant digits
  ## BMDU and BMDL could be NA. Replacing NA with -1.
  my_gene$bmdl[is.na(my_gene$bmdl)] <- (-1)
  gene_sub$bmdl[is.na(gene_sub$bmdl)] <- (-1)
  bmdl_check <- all(abs(my_gene$bmdl - gene_sub$bmdl) < 1e-5)
  expect_true(bmdl_check)

  my_gene$bmdu[is.na(my_gene$bmdu)] <- (-1)
  gene_sub$bmdu[is.na(gene_sub$bmdu)] <- (-1)
  bmdu_check <- all(abs(my_gene$bmdu - gene_sub$bmdu) < 1e-5)
  expect_true(bmdu_check)

  # check if the differences in hit-calls exceed a threshold value
  hitcall_check <- all(abs(my_gene$hitcall - gene_sub$hitcall) < 1e-5)
  expect_true(hitcall_check)

  # check if the differences in top_over_cutoff exceed a threshold value
  my_gene$top_over_cutoff[is.na(my_gene$top_over_cutoff)] <- (-1)
  gene_sub$top_over_cutoff[is.na(gene_sub$top_over_cutoff)] <- (-1)
  top_cutoff_check <- all(abs(my_gene$top_over_cutoff - gene_sub$top_over_cutoff) < 1e-5)
  expect_true(hitcall_check)

  ## check if the differences in top exceed a threshold value
  my_gene$top[is.na(my_gene$top)] <- (-1)
  gene_sub$top[is.na(gene_sub$top)] <- (-1)
  top_check <- all(abs(my_gene$top - gene_sub$top) < 1e-5)
  expect_true(top_check)

  ## check if the differences in AC50 exceed a threshold value
  my_gene$ac50[is.na(my_gene$ac50)] <- (-1)
  gene_sub$ac50[is.na(gene_sub$ac50)] <- (-1)
  ac50_check <- all(abs(my_gene$ac50 - gene_sub$ac50) < 1e-5)
  expect_true(ac50_check)

})


