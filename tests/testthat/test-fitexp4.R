test_that("fitexp4 works", {
  ex_conc <- c(0.01, 0.1, 0.5, 1, 3, 10, 30)
  resp <- exp4(ps = c(tp = 895,ga = 15,er = 0.1),x = ex_conc)
  expect_equal(fitexp4(ex_conc, resp)$tp, 805.5)
  expect_equal(fitexp4(ex_conc, resp)$ga, 13.11, tolerance=1e-2)
})
