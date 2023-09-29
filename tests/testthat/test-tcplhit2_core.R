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

test_that("tcplhit2 BMD boundary checks", {
  data("mc3")

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
  output_before <- tcplhit2_core(params, conc, resp, cutoff, onesd,
                                 bmed=0, conthits=T, aicc=F)
  output_after <- tcplhit2_core(params, conc, resp, cutoff, onesd,
                                bmed=0, conthits=T, aicc=F, bmd_up_bnd = 2)

  expect_equal(output_before$bmd, 299.9274, tolerance = 1e-3)
  expect_equal(output_after$bmd, 160, tolerance = 1e-3)

  expect_equal(output_before$bmdl, 217.343, tolerance = 1e-3)
  expect_equal(output_after$bmdl, 77.41564, tolerance = 1e-3)

  expect_equal(output_before$bmdu, 475.0638, tolerance = 1e-3)
  expect_equal(output_after$bmdu, 335.1364, tolerance = 1e-3)
})
