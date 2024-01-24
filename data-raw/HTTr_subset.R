## This R script should be ran from the command line using
## R CMD BATCH data-raw/HTTr_subset.R

## Script used to create the HTTr subsets at signature and gene levels for unit tests.

## The original masked HTPP data sets and the masking steps can be found under the following folder:
## "<...>\NCCT_ToxCast\Derik Haggard\HTTr\tcplfit2_HTTr_unitTestData"

## load HTTr output files
load("~/CompTox-ToxCast-tcplFit2/data-raw/HTTr/httr_inputData_MASKED.RData")
load("~/CompTox-ToxCast-tcplFit2/data-raw/HTTr/httr_tcplOutput_MASKED.RData")

## load necessary package
library(dplyr)

## Define the helper function for random sampling
init_sampler <- function(data, size) {
  # confident actives and inactives
  inactives <- data %>% filter(hitcall < 0.1 & top_over_cutoff < 0.9)
  actives <- data %>% filter(hitcall > 0.9 & top_over_cutoff > 1.5)

  # cases would not be considered active, but might have some with top_over_cutoff > 1
  noise_inactives <- data %>% filter(hitcall >= 0.1 & hitcall <= 0.7)

  # general borderline cases
  borderlines <- data %>% filter(top_over_cutoff > 0.9 & top_over_cutoff < 1.2) %>%
    filter(hitcall > 0.7 & hitcall < 0.95)
  # break down the borderline cases into borderline actives and inactives
  borderline_actives <- borderlines %>% filter(hitcall > 0.9 & hitcall < 0.95) %>%
    filter(top_over_cutoff > 1 & top_over_cutoff < 1.5)
  borderline_inactives <- borderlines %>% filter(hitcall < 0.9 | top_over_cutoff < 1)

  sub1 <- inactives[sample(x = nrow(inactives), size  = size), ]
  sub2 <- actives[sample(x = nrow(actives), size  = size), ]
  sub3 <- borderline_actives[sample(x = nrow(borderline_actives), size  = size), ]
  sub4 <- borderline_inactives[sample(x = nrow(borderline_inactives), size  = size), ]
  sub5 <- noise_inactives[sample(x = nrow(noise_inactives), size = size),]

  return(rbind(sub1, sub2, sub3, sub4, sub5))
}

## Create signature-level subset
set.seed(2233)
signature_sub <- init_sampler(SIGNATURE_CR, 25)

## Added some more samples for under-represented curve fit groups
all_active <- SIGNATURE_CR[SIGNATURE_CR$hitcall > 0.9 & SIGNATURE_CR$top_over_cutoff > 1.5,]
exp2_addins <- all_active[all_active$fit_method == "exp2", ]
pow_addins <- all_active[all_active$fit_method == "pow", ]

## Final signature subset
set.seed(63)
signature_sub <- rbind(signature_sub,
                       all_active[all_active$fit_method == "exp3",],
                       exp2_addins[sample(x=nrow(exp2_addins),size = 3),],
                       pow_addins[sample(x=nrow(pow_addins),size = 3),])


## Create gene-level subset
set.seed(1233)
gene_sub <- init_sampler(GENE_CR, 25)

## Added some more samples for under-represented curve fit groups
all_active <- GENE_CR[GENE_CR$hitcall > 0.9 & GENE_CR$top_over_cutoff > 1.5,]
poly2_acitve <- all_active[all_active$fit_method == "poly2", ]
exp3_active <- all_active[all_active$fit_method == "exp3", ]

## Final gene subset
set.seed(1132)
gene_sub <- rbind(gene_sub,
                  poly2_acitve[sample(x=nrow(poly2_acitve),size = 3),],
                  exp3_active[sample(x=nrow(exp3_active),size = 3),])

## Select the corresponding input data
signature_input <- NULL

for (i in 1:nrow(signature_sub)) {
  this.trt <- signature_sub[i, "trt"]
  this.sig <- signature_sub[i, "signature"]
  rows <- signaturescoremat[signaturescoremat$trt == this.trt &
                              signaturescoremat$signature == this.sig,]
  signature_input <- rbind(signature_input, rows)
}

gene_input <- NULL

for (i in 1:nrow(gene_sub)) {
  this.trt <- gene_sub[i, "trt"]
  this.gene <- gene_sub[i, "gene"]
  rows <- genemat[genemat$trt == this.trt &
                    genemat$gene == this.gene,]
  gene_input <- rbind(gene_input, rows)
}

# Data sets created are intermediate subset files that do not need to be tracked.
# Save them to the "HTTr' sub-directory.
# A separate R script will pull in internal R unit tests data and export sysdata.rda.
save(signature_sub, gene_sub,
     signature_input, gene_input,
     file = "~/CompTox-ToxCast-tcplFit2/data-raw/HTTr/httr_unittest.RData")
utils::sessionInfo()

