test_that("fitpoly2 works", {
  ex_conc <- c(0.01, 0.1, 0.5, 1, 3, 10, 30)
  resp <- poly2(ps = c(a = 0.13,b = 2,er = 0.1), x = ex_conc)
  expect_equal(fitpoly2(ex_conc, resp)$a, 0.13)
  expect_equal(fitpoly2(ex_conc, resp)$b, 2)
})
