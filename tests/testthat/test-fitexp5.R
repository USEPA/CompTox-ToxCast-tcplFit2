test_that("fitexp5 works", {
  ex_conc <- c(0.01, 0.1, 0.5, 1, 3, 10, 30)
  resp <- exp5(ps = c(tp = 793,ga = 6.25,p = 1.25,er = 0.1),x = ex_conc)
  expect_equal(fitexp5(ex_conc, resp)$tp, 793)
  expect_equal(fitexp5(ex_conc, resp)$ga, 6.25)
  expect_equal(fitexp5(ex_conc, resp)$p, 1.25)
})
