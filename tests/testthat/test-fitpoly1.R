test_that("fitpoly1 works", {
  ex_conc <- c(0.01, 0.1, 0.5, 1, 3, 10, 30)
  resp <- poly1(ps = c(a = 3.5,er = 0.1),x = ex_conc)
  expect_equal(fitpoly1(ex_conc, resp)$a, 3.5)
})
