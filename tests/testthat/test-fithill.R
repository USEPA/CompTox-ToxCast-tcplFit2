test_that("fithill works", {
  ex_conc <- c(0.01, 0.1, 0.5, 1, 3, 10, 30)
  resp <- hillfn(ps = c(tp = 750,ga = 5,p = 1.76,er = 0.1), x = ex_conc)
  expect_equal(fithill(ex_conc, resp)$tp, 750)
  expect_equal(fithill(ex_conc, resp)$ga, 5)
  expect_equal(fithill(ex_conc, resp)$p, 1.76)
})

