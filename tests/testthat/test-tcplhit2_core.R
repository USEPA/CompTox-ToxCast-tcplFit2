test_that("tcplhit2 works", {
  data("signatures")
  conc=as.numeric(str_split(signatures[1,"conc"],"\\|")[[1]])
  resp=as.numeric(str_split(signatures[1,"resp"],"\\|")[[1]])
  cutoff=signatures[1,"cutoff"]
  onesd=signatures[1,"onesd"]

  params = tcplfit2_core(conc, resp, cutoff = cutoff)
  output = tcplhit2_core(params, conc, resp, cutoff, onesd)

  expect_equal(output$fit_method, "exp4")
  expect_equal(output$hitcall, 0.99, tolerance = 1e-2)
  expect_equal(output$tp, 0.749, tolerance = 1e-2)
  expect_equal(output$ga, 9.59, tolerance = 1e-2)
})

test_that("tcplhit2 BMD boundary check", {
  # Preparation for tests on BMD boundary arguments
  # Simulate some data

  # The choice of X's is meant to imitate concentrations in real data
  # The use of normal noise departs from the assumption of t-distribution error,
  # mainly for simulation convenience. Normal distribution allows easier set
  # up for desired sd and bmr to meet different data cases.
  X <- rep(c(0.002, 0.003, 0.01, 0.03, 1, 2, 3, 5, 7, 10), each = 5)

  set.seed(167)
  Case_1 <- hillfn(ps = c(7.5,0.032,0.65,0.1),x = X) + rnorm(n = length(X), sd = 0.3)

  set.seed(135)
  Case_2 <- hillfn(ps = c(7.536,0.02958,0.7349,0.1),x = X) + rnorm(n = length(X), sd = 0.6)

  set.seed(918)
  Case_3 <- pow(ps = c(0.2,3.5627,0.1),x = X) + rnorm(n = length(X), sd = 0.5)

  set.seed(400)
  Case_4 <- poly1(ps = c(0.2,0.1),x = X) + rnorm(n = length(X), sd = 4.5)

  set.seed(66)
  Case_5 <- poly1(ps = c(0.25),x = X) + rnorm(n = length(X), sd = 4.5)

  Y <- rbind(Case_1, Case_2, Case_3, Case_4, Case_5)

  # Y contains 5 rows of simulated responses, each row corresponds to
  # one data cases in order:

  # 1. BMD < lower threshold
  # 2. lower threshold < BMD < lowest expt dose
  # 3. lowest expt dose < BMD < upper expt dose
  # 4. upper expt dose < BMD < upper threshold
  # 5. BMD > upper threshold

  df <- matrix(nrow = 5, ncol = 14)
  colnames(df) <- c("onesd", "cutoff", "bmd", "bmdu", "bmdl",
                    "combo1", "combo1u", "combo1l",
                    "combo2", "combo2u", "combo2l",
                    "combo3", "combo3u", "combo3l")
  bmad <- mad(Y[, 1:10])
  cutoff <- 3*bmad
  df[, "cutoff"] <- cutoff

  # fit each data case and hit-call with different argument combinations
  for (i in 1:5) {
    temp <- Y[i, 1:10]
    onesd <- sd(temp)

    df[i, "onesd"] <- onesd

    params <- tcplfit2_core(X, Y[i, ], cutoff,
                            force.fit = T, bidirectional = F)
    # No thresholds specified, no censoring
    output <- tcplhit2_core(params, X, Y[i,], cutoff, onesd,
                            bmed=0, conthits=T, aicc=F)

    df[i, "bmd"] <- output$bmd
    df[i, "bmdu"] <- output$bmdu
    df[i, "bmdl"] <- output$bmdl

    # For these tests we use a lower and upper BMD threshold multiple of 0.3 and
    # 5, respectively, to demonstrate use cases work for the simulated
    # concentration-response data. Threshold multiples are not necessarily meant
    # to align with the typical suggested values of 0.1 and 10 for lower and
    # upper threshold multiples, respectively.

    # Combo 1: Only the upper boundary threshold is specified
    combo1 <- tcplhit2_core(params, X, Y[i,], cutoff, onesd,
                            bmed=0, conthits=T, aicc=F, bmd_up_bnd = 5)

    # Combo 2: Only the lower boundary threshold is specified
    combo2 <- tcplhit2_core(params, X, Y[i,], cutoff, onesd,
                            bmed=0, conthits=T, aicc=F, bmd_low_bnd = 0.3)

    # Combo 3: Both lower and upper threshold is specified
    combo3 <- tcplhit2_core(params, X, Y[i,], cutoff, onesd,
                            bmed=0, conthits=T, aicc=F, bmd_up_bnd = 5, bmd_low_bnd = 0.3)

    # record BMDs
    df[i, "combo1"] <- combo1$bmd; df[i, "combo1u"] <- combo1$bmdu; df[i, "combo1l"] <- combo1$bmdl
    df[i, "combo2"] <- combo2$bmd; df[i, "combo2u"] <- combo2$bmdu; df[i, "combo2l"] <- combo2$bmdl
    df[i, "combo3"] <- combo3$bmd; df[i, "combo3u"] <- combo3$bmdu; df[i, "combo3l"] <- combo3$bmdl

  }
  df <- as.data.frame(cbind(Cases = seq(1,5,1), df))



  # Argument combination: Only the upper boundary threshold is specified
  # only case 5 where BMD > upper threshold should be affected

  # Case 1 unaffected
  expect_equal(df[1, "combo1"], df[1, "bmd"], tolerance = 1e-3)
  expect_equal(df[1, "combo1u"], df[1, "bmdu"], tolerance = 1e-3)
  expect_equal(df[1, "combo1l"], df[1, "bmdl"], tolerance = 1e-3)
  # Case 2 unaffected
  expect_equal(df[2, "combo1"], df[2, "bmd"], tolerance = 1e-3)
  expect_equal(df[2, "combo1u"], df[2, "bmdu"], tolerance = 1e-3)
  expect_equal(df[2, "combo1l"], df[2, "bmdl"], tolerance = 1e-3)
  # Case 3 unaffected
  expect_equal(df[3, "combo1"], df[3, "bmd"], tolerance = 1e-3)
  expect_equal(df[3, "combo1u"], df[3, "bmdu"], tolerance = 1e-3)
  expect_equal(df[3, "combo1l"], df[3, "bmdl"], tolerance = 1e-3)
  # Case 4 unaffected
  expect_equal(df[4, "combo1"], df[4, "bmd"], tolerance = 1e-3)
  expect_true(is.na(df[4, "combo1u"]))
  expect_equal(df[4, "combo1l"], df[4, "bmdl"], tolerance = 1e-3)
  # Case 5 - BMD will be censored to the upper threshold
  expect_equal(df[5, "combo1"], max(X)*5, tolerance = 1e-3)
  expect_true(is.na(df[5, "combo1u"]))
  expect_equal(df[5, "combo1l"], df[5, "bmdl"]-(df[5, "bmd"] - max(X)*5), tolerance = 1e-3)

  # Argument combination: Only the lower boundary threshold is specified
  # use a multiplier of 0.3
  # only case 1 where BMD < lower threshold should be affected

  # Case 1 - BMD will be censored to the lower threshold
  expect_equal(df[1, "combo2"], min(X)*0.3, tolerance = 1e-3)
  expect_equal(df[1, "combo2u"], df[1, "bmdu"] + (min(X)*0.3-df[1, "bmd"]) , tolerance = 1e-3)
  expect_equal(df[1, "combo2l"], df[1, "bmdl"] + (min(X)*0.3-df[1, "bmd"]), tolerance = 1e-3)
  # Case 2 unaffected
  expect_equal(df[2, "combo2"], df[2, "bmd"], tolerance = 1e-3)
  expect_equal(df[2, "combo2u"], df[2, "bmdu"], tolerance = 1e-3)
  expect_equal(df[2, "combo2u"], df[2, "bmdu"], tolerance = 1e-3)
  # Case 3 unaffected
  expect_equal(df[3, "combo2"], df[3, "bmd"], tolerance = 1e-3)
  expect_equal(df[3, "combo2u"], df[3, "bmdu"], tolerance = 1e-3)
  expect_equal(df[3, "combo2l"], df[3, "bmdl"], tolerance = 1e-3)
  # Case 4 unaffected
  expect_equal(df[4, "combo2"], df[4, "bmd"], tolerance = 1e-3)
  expect_true(is.na(df[4, "combo2u"]))
  expect_equal(df[4, "combo2l"], df[4, "bmdl"], tolerance = 1e-3)
  # Case 5 unaffected
  expect_equal(df[5, "combo2"], df[5, "bmd"], tolerance = 1e-3)
  expect_true(is.na(df[5, "combo2u"]))
  expect_equal(df[5, "combo2l"], df[5, "bmdl"], tolerance = 1e-3)

  # Argument combination: Both lower and upper threshold is specified
  # only case 1 and case 5 should be affected

  # Case 1 - BMD will be censored to the lower threshold
  expect_equal(df[1, "combo3"], min(X)*0.3, tolerance = 1e-3)
  expect_equal(df[1, "combo3u"], df[1, "bmdu"] + (min(X)*0.3-df[1, "bmd"]) , tolerance = 1e-3)
  expect_equal(df[1, "combo3u"], df[1, "bmdu"] + (min(X)*0.3-df[1, "bmd"]) , tolerance = 1e-3)
  # Case 2 unaffected
  expect_equal(df[2, "combo3"], df[2, "bmd"], tolerance = 1e-3)
  expect_equal(df[2, "combo3u"], df[2, "bmdu"], tolerance = 1e-3)
  expect_equal(df[2, "combo3l"], df[2, "bmdl"], tolerance = 1e-3)
  # Case 3 unaffected
  expect_equal(df[3, "combo3"], df[3, "bmd"], tolerance = 1e-3)
  expect_equal(df[3, "combo3u"], df[3, "bmdu"], tolerance = 1e-3)
  expect_equal(df[3, "combo3l"], df[3, "bmdl"], tolerance = 1e-3)
  # Case 4 unaffected
  expect_equal(df[4, "combo3"], df[4, "bmd"], tolerance = 1e-3)
  expect_true(is.na(df[4, "combo3u"]))
  expect_equal(df[4, "combo3l"], df[4, "bmdl"], tolerance = 1e-3)
  # Case 5  - BMD will be censored to the upper threshold
  expect_equal(df[5, "combo3"], max(X)*5, tolerance = 1e-3)
  expect_true(is.na(df[5, "combo3u"]))
  expect_equal(df[5, "combo3l"], df[5, "bmdl"]-(df[5, "bmd"] - max(X)*5), tolerance = 1e-3)
})

