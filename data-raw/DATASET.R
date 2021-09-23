## code to prepare `DATASET` dataset goes here
signatures <- read.csv(file = "./data-raw/signatures.csv", stringsAsFactors=FALSE)
usethis::use_data(signatures, overwrite = TRUE)
