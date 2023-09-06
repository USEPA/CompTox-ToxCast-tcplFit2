test_that("concRespCore works", {
  data("signatures")
  row = list(conc=as.numeric(str_split(signatures[1,"conc"],"\\|")[[1]]),
             resp=as.numeric(str_split(signatures[1,"resp"],"\\|")[[1]]),
             bmed=0,
             cutoff=signatures[1,"cutoff"],
             onesd=signatures[1,"onesd"],
             name=signatures[1,"name"],
             assay=signatures[1,"signature"])
  out = concRespCore(row)

  expect_equal(out$fit_method, "exp4")
  expect_equal(out$hitcall, 0.99, tolerance = 1e-2)
  expect_equal(out$tp, 0.749, tolerance = 1e-2)
  expect_equal(out$ga, 9.59, tolerance = 1e-2)

})
