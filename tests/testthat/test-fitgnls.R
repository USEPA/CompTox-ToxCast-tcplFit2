test_that("fitgnls works", {
  ex_conc <- c(0.01, 0.1, 0.5, 1, 3, 10, 30)
  resp <- gnls(ps = c(tp = 750,ga = 15,p = 1.45,la = 50,q = 1.34, er = 0.1),x = ex_conc)
  expect_equal(fitgnls(ex_conc, resp)$tp, 435, tolerance=1e-2)
  expect_equal(fitgnls(ex_conc, resp)$ga, 8.66, tolerance=1e-2)
  expect_equal(fitgnls(ex_conc, resp)$p, 1.63, tolerance=1e-2)
  expect_equal(fitgnls(ex_conc, resp)$la, 273.86, tolerance=1e-2)
  expect_equal(fitgnls(ex_conc, resp)$q, 1.32, tolerance=1e-2)
})
