## R script to pull in internal R unit tests data and export sysdata.rda

## Pull in HTPP unit test data subsets
load("~/CompTox-ToxCast-tcplFit2/data-raw/HTPP/htpp_unittest.RData")

## Pull in HTTr unit test data subsets
load("~/CompTox-ToxCast-tcplFit2/data-raw/HTTr/httr_unittest.RData")

## Export internal sysdata.rda
usethis::use_data(htpp_global_subset, htpp_cat_subset, htpp_feature_subset,
                  htpp_global_input, htpp_cat_input, htpp_feature_input,
                  CONTROL_GMAH, CONTROL_CMAH, CONTROL_FMAH,
                  signature_sub, gene_sub,
                  signature_input, gene_input,
                  internal = TRUE, overwrite = T)
## Include session info
utils::sessionInfo()
