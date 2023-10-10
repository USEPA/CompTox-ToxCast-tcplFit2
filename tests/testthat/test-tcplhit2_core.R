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
