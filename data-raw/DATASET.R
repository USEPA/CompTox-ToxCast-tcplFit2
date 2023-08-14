## code to prepare `DATASET` dataset goes here
signatures <- read.csv(file = "./data-raw/signatures.csv", stringsAsFactors=FALSE)
usethis::use_data(signatures, overwrite = TRUE)

## For example tcpl data
#connect to db first
mc0 <- tcplLoadData(lvl = 0, fld = "acid", val = 1829)
usethis::use_data(mc0)
