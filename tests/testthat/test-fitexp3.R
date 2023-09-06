test_that("fitexp3 works", {
  ex_conc <- c(0.01, 0.1, 0.5, 1, 3, 10, 30)
  resp <- exp3(ps = c(a = 1.67,b = 12.5,p = 0.87,er = 0.1), x = ex_conc)
  expect_equal(fitexp3(ex_conc, resp)$a, 1.67)
  expect_equal(fitexp3(ex_conc, resp)$b, 12.5)
  expect_equal(fitexp3(ex_conc, resp)$p, 0.87)
})

