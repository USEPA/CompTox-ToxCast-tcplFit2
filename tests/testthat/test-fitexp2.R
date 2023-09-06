test_that("fitexp2 works", {
  ex_conc <- c(0.01, 0.1, 0.5, 1, 3, 10, 30)
  resp <- exp2(ps = c(a = 0.45,b = 13.5,er = 0.1), x = ex_conc)
  expect_equal(fitexp2(ex_conc, resp)$a, 0.45)
  expect_equal(fitexp2(ex_conc, resp)$b, 13.5)
})
