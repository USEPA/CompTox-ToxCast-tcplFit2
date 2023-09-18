test_that("fitpow works", {
  ex_conc <- c(0.01, 0.1, 0.5, 1, 3, 10, 30)
  resp <- pow(ps = c(a = 1.5,p = 1.34,er = 0.1),x = ex_conc)
  expect_equal(fitpow(ex_conc, resp)$a, 1.5)
  expect_equal(fitpow(ex_conc, resp)$p, 1.34)
})
