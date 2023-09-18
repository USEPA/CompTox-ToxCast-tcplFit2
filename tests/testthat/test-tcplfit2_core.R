test_that("tcplfit2_core works", {
  data("signatures")
  conc=as.numeric(str_split(signatures[1,"conc"],"\\|")[[1]])
  resp=as.numeric(str_split(signatures[1,"resp"],"\\|")[[1]])
  cutoff=signatures[1,"cutoff"]
  output <- tcplfit2_core(conc, resp, cutoff = cutoff)

  #check hill result
  expect_equal(output[["hill"]]$tp, 0.6795, tolerance = 1e-3)
  expect_equal(output[["hill"]]$ga, 8.3549, tolerance = 1e-3)
  expect_equal(output[["hill"]]$p, 2.6533, tolerance = 1e-3)

  #check gain-loss result
  expect_equal(output[["gnls"]]$tp, 0.7883, tolerance = 1e-3)
  expect_equal(output[["gnls"]]$ga, 9.0696, tolerance = 1e-3)
  expect_equal(output[["gnls"]]$p, 2.5650, tolerance = 1e-3)
  expect_equal(output[["gnls"]]$la, 286.8063, tolerance = 1e-3)
  expect_equal(output[["gnls"]]$q, 0.8523, tolerance = 1e-3)

  #check poly1 result
  expect_equal(output[["poly1"]]$a, 0.02295, tolerance = 1e-3)

  # check poly2 result
  expect_equal(output[["poly2"]]$a, 274.2645, tolerance = 1e-3)
  expect_equal(output[["poly2"]]$b, 11976.98, tolerance = 1e-3)

  # check power result
  expect_equal(output[["pow"]]$a, 0.0687, tolerance = 1e-3)
  expect_equal(output[["pow"]]$p, 0.6722, tolerance = 1e-3)

  # check exp2 result
  expect_equal(output[["exp2"]]$a, 201.6733, tolerance = 1e-3)
  expect_equal(output[["exp2"]]$b, 8800.5, tolerance = 1e-3)

  # check exp3 result
  expect_equal(output[["exp3"]]$a, 2.7861, tolerance = 1e-3)
  expect_equal(output[["exp3"]]$b, 269.5222, tolerance = 1e-3)
  expect_equal(output[["exp3"]]$p, 0.6895, tolerance = 1e-3)

  # check exp4 result
  expect_equal(output[["exp4"]]$tp, 0.7491, tolerance = 1e-3)
  expect_equal(output[["exp4"]]$ga, 9.5911, tolerance = 1e-3)

  # check exp5 result
  expect_equal(output[["exp5"]]$tp, 0.6572, tolerance = 1e-3)
  expect_equal(output[["exp5"]]$ga, 8.4072, tolerance = 1e-3)
  expect_equal(output[["exp5"]]$p, 2.2023, tolerance = 1e-3)

})
