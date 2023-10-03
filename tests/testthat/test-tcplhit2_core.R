test_that("tcplhit2 works", {
  data("signatures")
  conc=as.numeric(str_split(signatures[1,"conc"],"\\|")[[1]])
  resp=as.numeric(str_split(signatures[1,"resp"],"\\|")[[1]])
  cutoff=signatures[1,"cutoff"]
  onesd=signatures[1,"onesd"]

  params = tcplfit2_core(conc, resp, cutoff = cutoff)
  output = tcplhit2_core(params, conc, resp, cutoff, onesd)

  expect_equal(output$fit_method, "exp4")
  expect_equal(output$hitcall, 0.99, tolerance = 1e-2)
  expect_equal(output$tp, 0.749, tolerance = 1e-2)
  expect_equal(output$ga, 9.59, tolerance = 1e-2)
})

test_that("tcplhit2 BMD upper boundary censor checks", {
  data("mc3")

  # use example data from mc3
  spid <- unique(mc3$spid)[26]
  ex_df <- mc3[is.element(mc3$spid,spid),]
  conc <- 10**ex_df$logc # back-transforming concentrations on log10 scale
  resp <- ex_df$resp

  temp <- mc3[mc3$logc<= -2,"resp"]
  bmad <- mad(temp)
  onesd <- sd(temp)
  cutoff <- 3*bmad

  params <- tcplfit2_core(conc, resp, cutoff,
                          force.fit = T, bidirectional = F)
  # hit-calling with and without using the BMD boundary argument
  output_before <- tcplhit2_core(params, conc, resp, cutoff, onesd,
                                 bmed=0, conthits=T, aicc=F)
  # Estimated BMD is greater than the upper 'threshold dose'.
  output_after <- tcplhit2_core(params, conc, resp, cutoff, onesd,
                                bmed=0, conthits=T, aicc=F, bmd_up_bnd = 2)
  # Estimated BMD is below the upper 'threshold dose' and
  # above the 'reference dose' (highest expt dose).
  output_after_ten <- tcplhit2_core(params, conc, resp, cutoff, onesd,
                                bmed=0, conthits=T, aicc=F, bmd_up_bnd = 10)

  # checks for bmd estimates
  expect_equal(output_before$bmd, 299.9274, tolerance = 1e-3)
  expect_equal(output_after$bmd, 160, tolerance = 1e-3)
  expect_equal(output_after_ten$bmd, 299.9274, tolerance = 1e-3)

  # checks for bmd confident lower bound
  expect_equal(output_before$bmdl, 217.343, tolerance = 1e-3)
  expect_equal(output_after$bmdl, 77.41564, tolerance = 1e-3)
  expect_equal(output_after_ten$bmdl, 217.343, tolerance = 1e-3)

  # checks for bmd confident upper bound
  expect_equal(output_before$bmdu, 475.0638, tolerance = 1e-3)
  expect_equal(output_after$bmdu, 335.1364, tolerance = 1e-3)
  expect_equal(output_after_ten$bmdu, 475.0638, tolerance = 1e-3)
})

test_that("tcplhit2 BMD lower boundary censor checks", {
  data("mc3")

  temp <- mc3[mc3$logc<= -2,"resp"]
  bmad <- mad(temp)
  onesd <- sd(temp)
  cutoff <- 3*bmad

  # load example data from mc3
  spid <- unique(mc3$spid)[94]
  ex_df <- mc3[is.element(mc3$spid,spid),]
  conc <- 10**ex_df$logc # back-transforming concentrations on log10 scale
  resp <- ex_df$resp
  conc2 <- conc[conc>0.41]
  resp2 <- resp[which(conc>0.41)]

  params <- tcplfit2_core(conc2, resp2, cutoff,
                          force.fit = T, bidirectional = F)

  # hit-calling with and without using the BMD boundary argument
  output_before <- tcplhit2_core(params, conc2, resp2, cutoff, onesd,
                                 bmed=0, conthits=T, aicc=F)
  # Estimated BMD is less than the lower 'threshold dose'.
  output_after <- tcplhit2_core(params, conc2, resp2, cutoff, onesd,
                                bmed=0, conthits=T, aicc=F, bmd_low_bnd = 0.8)
  # Estimated BMD is above the lower 'threshold dose' and
  # below the 'reference dose' (lowest expt dose).
  output_after_ten <- tcplhit2_core(params, conc2, resp2, cutoff, onesd,
                                bmed=0, conthits=T, aicc=F, bmd_low_bnd = 0.1)

  # checks for bmd estimates
  expect_equal(output_before$bmd, 0.302, tolerance = 1e-3)
  expect_equal(output_after$bmd, 0.48, tolerance = 1e-3)
  expect_equal(output_after_ten$bmd, 0.302, tolerance = 1e-3)

  # checks for bmd confident lower bound
  expect_equal(output_before$bmdl, 0.1048, tolerance = 1e-3)
  expect_equal(output_after$bmdl, 0.2826, tolerance = 1e-3)
  expect_equal(output_after_ten$bmdl, 0.1048, tolerance = 1e-3)

  # checks for bmd confident upper bound
  expect_equal(output_before$bmdu, 0.815, tolerance = 1e-3)
  expect_equal(output_after$bmdu, 0.9928, tolerance = 1e-3)
  expect_equal(output_after_ten$bmdu, 0.815, tolerance = 1e-3)
})

test_that("tcplhit2 BMD boundary censor checks", {
  data("mc3")

  temp <- mc3[mc3$logc<= -2,"resp"]
  bmad <- mad(temp)
  onesd <- sd(temp)
  cutoff <- 3*bmad

  # find some example data from mc3
  spid <- unique(mc3$spid)[78]
  ex_df <- mc3[is.element(mc3$spid,spid),]
  conc <- 10**ex_df$logc # back-transforming concentrations on log10 scale
  resp <- ex_df$resp

  params <- tcplfit2_core(conc, resp, cutoff,
                          force.fit = T, bidirectional = F)

  # Estimated BMD is in the experimental concentration range
  output_before <- tcplhit2_core(params, conc, resp, cutoff, onesd,
                                 bmed=0, conthits=T, aicc=F)

  # Estimated BMD is below the upper 'reference dose' (highest expt dose)
  # and is above the lower 'reference dose' (lowest expt. dose).
  output_after <- tcplhit2_core(params, conc, resp, cutoff, onesd,
                                bmed=0, conthits=T, aicc=F,
                                bmd_low_bnd = 0.8, bmd_up_bnd = 2)

  # checks for bmd estimates, expect them to be the same
  expect_equal(output_before$bmd, 21.63718, tolerance = 1e-3)
  expect_equal(output_after$bmd, 21.63718, tolerance = 1e-3)

  # checks for bmd confident lower bound, expect them to be the same
  expect_equal(output_before$bmdl, 21.14863, tolerance = 1e-3)
  expect_equal(output_after$bmdl, 21.14863, tolerance = 1e-3)

  # checks for bmd confident upper bound, expect them to be the same
  expect_equal(output_before$bmdu, 22.21793, tolerance = 1e-3)
  expect_equal(output_after$bmdu, 22.21793, tolerance = 1e-3)
})

