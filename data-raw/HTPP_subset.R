## This R script should be ran from the command line using
## R CMD BATCH data-raw/HTPP_subset.R

## Script used to create the HTPP subsets at global, category, and feature levels for unit tests

## The original HTPP data sets can be found under the following folder:
## "L:\Lab\NCCT_ToxCast\Derik Haggard\HTPP\tcplfit2_HTPP_unitTestData\res_htpp_mcf7_apcra_unitTestData"
## load HTPP data
load("~/CompTox-ToxCast-tcplFit2/data-raw/HTPP/htpp_inputData.RData", verbose = T)
load("~/CompTox-ToxCast-tcplFit2/data-raw/HTPP/htpp_tcplOutput.RData", verbose = T)

## load necessary package
library(dplyr)

## chemicals selected for global subset
## Try to include at least one response from each curve, and balance active and inactive responses
## poly2 is not available in the corresponding output data
global_chem_id <- c("EPAPLT0617B02", # poly1 - borderline, low hit-call 0.2, top/cutoff 1.14
                    "EPAPLT0617C08", # poly1 - borderline, low hit-call 0.11, top/cutoff 0.767
                    "EPAPLT0617E10", # the only exp3, robust active
                    "EPAPLT0617D04", # exp2, active
                    "EPAPLT0617B01", # exp4, active
                    "EPAPLT0617C02", # exp5, active
                    "EPAPLT0617A11", # pow, active
                    "EPAPLT0617A09", # hill, active
                    "EPAPLT0617B09", # exp4, inactive, NA bmd
                    "EPAPLT0617D01", # poly1 inactive
                    "EPAPLT0617B08", # none, inactive, NA bmd
                    "EPAPLT0617C03"  # poly1, inactive
)

## subset the input and output files, n = 12
htpp_global_subset <- tcpl_global[tcpl_global$chem_id %in% global_chem_id,]
htpp_global_input <- htpp_global_mah[htpp_global_mah$chem_id %in% global_chem_id,]


## chemicals selected for category subset
## There are 49 unique categories, most of them start with prefixes:
## AGP, DNA, ER, Mito, and RNA, and there are two called Shape and Position.
## Strategy:
## Group categories into buckets by their prefixes,
## and try to sample an equal number of responses from each of the buckets.
## want to include robust active and inactive cases, and borderline cases.

## Steps:
## Find borderline cases. Look for chemicals that have `top_over_cutoff` close to 1.

## cat_borderline <- tcpl_cat %>% filter(!is.na(top_over_cutoff)) %>%
##   filter(top_over_cutoff > 0.9 & top_over_cutoff < 1.5) %>%
##   filter(stype == "test sample")

## For these borderline cases, the hit-call ranges between 0.137 and 0.967.
## Selected 5 samples: selected the sample with the lowest hit-call among these borderline
## cases, selected one sample which has a higher hit-call (but below 0.9),
## selected one sample with a hit-call close 0.5 (the "random guess" threshold),
## and selected two samples with hit-call between 0.15 and 0.5.

## Then moved on to fill the rest with some robust cases from buckets that have not been sampled on,
## prioritize models that have not been included yet.
cat_chem_id <- c( "EPAPLT0617C09-AGP_Intensity_Ring", # exp4 borderline, hitcall 0.5
                  "EPAPLT0617A11-ER_Texture_Ring", # exp5 borderline, hitcall 0.826
                  "EPAPLT0617B09-ER_Symmetry_Cells", #poly1 borderline, hitcall 0.137
                  "EPAPLT0617B04-DNA_Profile_Cytoplasm", # exp2 borderline, hitcall 0.1829
                  "EPAPLT0617B01-AGP_Profile_Nuclei", #exp5 borderline, hitcall 0.4043162

                  "EPAPLT0617A10-DNA_Radial_Nuclei", #poly1, inactive
                  "EPAPLT0617D04-Mito_Profile_Cytoplasm", # poly2 active
                  "EPAPLT0617B02-Mito_Profile_Cytoplasm", # exp5 inactive
                  "EPAPLT0617C08-RNA_Symmetry_Nuclei", # "none" model, inactive
                  "EPAPLT0617C02-RNA_Texture_Nuclei", # hill active
                  "EPAPLT0617E10-Shape", # exp3 active
                  "EPAPLT0617B08-Shape", # poly1 inactive
                  "EPAPLT0617A03-Position", # poly1 inactive
                  "EPAPLT0617E10-Position" # pow active
)

## subset the input and output files, n = 14
htpp_cat_subset <- NULL
htpp_cat_input <- NULL

for (each in cat_chem_id) {
  chem <- strsplit(each, split="-")[[1]][1]
  cat <- strsplit(each, split="-")[[1]][2]
  row <- tcpl_cat[tcpl_cat$chem_id == chem & tcpl_cat$endpoint == cat,]
  sub <- htpp_cat_mah[htpp_cat_mah$chem_id == chem & htpp_cat_mah$category_name_r == cat,]
  htpp_cat_subset <- rbind(htpp_cat_subset, row)
  htpp_cat_input <- rbind(htpp_cat_input, sub)
}


## chemicals selected for feature subset
## There are 1300 unique feature. The strategy is just to consider covering
## all types of curves and have a good mix of robust and borderline cases.
## Create three subsets from the output data: the very active ones and the
## inactive ones, and the ones on borderline.
## When selecting samples from subset 3, keep a balance between number of borderline with
## high hit-call and low hit=call.

## temp <- tcpl_feature[tcpl_feature$stype == "test sample",]
## sub1 <- temp %>% filter(hitcall < 0.1)
## sub2 <- temp %>% filter(hitcall > 0.9)
## sub3 <- temp %>% filter(top_over_cutoff > 0.9 & temp$top_over_cutoff < 1.5) %>%
##   filter(hitcall > 0.1 & hitcall < 0.9)

feature_chem_id <- c(
  ## from sub1
  "EPAPLT0617B08-f_266", # exp2 - inactive
  "EPAPLT0617A01-f_1028", # exp4 - inactive
  "EPAPLT0617C10-f_1141",  # exp5 - inactive
  "EPAPLT0617A04-f_1105",  # hill - inactive
  "EPAPLT0617A03-f_1164",  # pow - inactive
  "EPAPLT0617A01-f_124", # poly1 - inactive
  ## from sub2
  "EPAPLT0617D04-f_406", #  exp2 active
  "EPAPLT0617E10-f_75",  #  exp3 active
  "EPAPLT0617A01-f_1137", #  exp4 active
  "EPAPLT0617D04-f_1180", #  exp5 active
  "EPAPLT0617B01-f_517", # hill active
  "EPAPLT0617C02-f_330", #  poly1 active
  "EPAPLT0617D04-f_549", #  poly2 active
  "EPAPLT0617A05-f_502", #  pow active
  ## from sub3
  "EPAPLT0617A03-f_1142", # hill, borderline, hit-call around 0.69
  "EPAPLT0617A03-f_1224", # exp2, borderline, hit-call around 0.5
  "EPAPLT0617C08-f_1185", # poly1 borderline, hit-call around 0.209
  "EPAPLT0617B08-f_1297", # exp5, borderline, hit-call around 0.309
  "EPAPLT0617D06-f_37", # exp4, borderline, hit-call around 0.722
  "EPAPLT0617B06-f_393" # pow, borderline, hit-call around 0.787
)

## subset the input and output files, n = 20
htpp_feature_subset <- NULL
htpp_feature_input <- NULL

for (each in feature_chem_id) {
  chem <- strsplit(each, split="-")[[1]][1]
  fname <- strsplit(each, split="-")[[1]][2]
  row <- tcpl_feature[tcpl_feature$chem_id == chem & tcpl_feature$endpoint == fname,]
  sub <- htpp_well_norm %>% filter(chem_id == chem) %>% select(c(1:23, which(colnames(htpp_well_norm) == fname)))
  colnames(sub)[ncol(sub)] <- "d"
  sub$Feature <- fname
  htpp_feature_subset <- rbind(htpp_feature_subset, row)
  htpp_feature_input <- rbind(htpp_feature_input, sub)
}

## Calculate Vehicle Control Information for Three Levels
## Global - uses the same cutoff and BMED for the whole data set
## Category - each category uses specific cutoff and BMED values
## Feature - each feature uses specific cutoff and BMED values

CONTROL_GMAH <- htpp_global_mah %>%
  dplyr::filter(stype == "test sample") %>%
  dplyr::filter(dose_level == 1|dose_level == 2) %>%
  dplyr::summarise(BMED = mean(d),CUTOFF = sd(d),ONESD = sd(d)/1.349)

CONTROL_CMAH <- htpp_cat_mah %>%
  dplyr::filter(stype == "test sample") %>%
  dplyr::group_by(category_name_r) %>%
  dplyr::filter(dose_level == 1|dose_level == 2) %>%
  dplyr::summarise(BMED = mean(d),CUTOFF = sd(d),ONESD = sd(d)/1.349)


df <- htpp_well_norm %>%
  dplyr::filter(stype == "test sample")  %>%
  dplyr::filter(dose_level == 1|dose_level == 2)

CONTROL_FMAH <- NULL
for (i in 24:ncol(htpp_well_norm)){
  temp <- df %>% select(i)
  CONTROL_FMAH <- rbind(CONTROL_FMAH,
  data.frame(BMED = mean(temp[[1]]),
             CUTOFF = sd(temp[[1]]),
             ONESD = sd(temp[[1]])/1.349,
             fname = colnames(htpp_well_norm)[i]))
  }

# save all the data frames created as internal data files
usethis::use_data(htpp_global_subset, htpp_cat_subset, htpp_feature_subset,
                  htpp_global_input, htpp_cat_input, htpp_feature_input,
                  CONTROL_GMAH, CONTROL_CMAH, CONTROL_FMAH, internal = TRUE, overwrite = T)
utils::sessionInfo()

