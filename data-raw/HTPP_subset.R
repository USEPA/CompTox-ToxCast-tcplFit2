## Scripts used to create HTPP subsets at global, category, and feature levels for unit test

## load HTTP data
load(system.file("data-raw/HTPP", "htpp_inputData.RData", package = "tcplfit2"))

## chemicals selected for global subset
## Try to include at least one response from each curve, and balance active and inactive responses
## poly2 and exp3 are not available in the corresponding output data
global_chem_id <- c("EPAPLT0202E05", # the only exp2 - active
                    "EPAPLT0202F18", # the only exp5 - active
                    "EPAPLT0202G15", # power - active
                    "ETOP", # hill - active
                    "EPAPLT0202A11", # exp4 - inactive
                    "DEX", # exp4 - active
                    "EPAPLT0202D18", # poly1 - inactive
                    "EPAPLT0202F22", # poly1 - inactive
                    "EPAPLT0202C01", # csnt - inactive
                    "EPAPLT0202F11" # poly1 - hit-call around 0.5, top/cutoff 1.2
)

## create the subset, n = 10
htpp_global_subset <- tcpl_global[tcpl_global$chem_id %in% global_chem_id]


## chemicals selected for category subset
## There are 49 unique categories, most of them start with prefixes:
## AGP, DNA, ER, Mito, and RNA, and there's one called Shape.
## Group categories into six buckets by the prefixes.
## Strategy:
## Trying to sample from all six buckets and cover all models,
## want to include robust active and inactive cases, and borderline cases.
## exp3 is not available in the corresponding output data

## Ran table(tcpl_cat$fit_method), there's only only poly2 in the output, include it in the subset.
## Then move on to find borderline cases because they are rarer.
## Look for chemicals that have top_over_cutoff close to 1, meaning the top is just below or just above the cutoff

##cat_borderline <- tcpl_cat %>% filter(top_over_cutoff > 0.9 & top_over_cutoff < 1.5) %>%
  ##filter(stype == "test sample")

## For these cases, the hit-call is between 0.2 and 0.8.
## Select EPAPLT0202E05 from DNA_Axial_Nuclei with a hit-call of 0.89, borderline high hit-call
## Select EPAPLT0202E12 from Mito_Intensity_Ring with a hit-call of 0.5, right at the "random guess" threshold
## Select EPAPLT0202A06 from  ER_Profile_Cytoplasm with a hit-call of 0.203, borderline low hit-call

## Move on to select some robust cases from Shape, RNA and AGP buckets, prioritize
## models that have not been included yet. Then select one more case from Mito,
## and DNA categories.
## Example codes I used to filter and select:
## temp <- tcpl_cat[grepl("RNA", tcpl_cat$endpoint), ]
## table(temp$fit_method) check if there's unique cases in this
## temp <- temp %>% filter(stype == "test sample" & fit_method == "exp4") # select the first or the second row

cat_chem_id <- c(   "EPAPLT0202A05-ER_Radial_Cells", # only poly2, robust high hit-call
                    "EPAPLT0202E05-DNA_Axial_Nuclei", # pow borderline high hitcall
                    "EPAPLT0202A06-ER_Profile_Cytoplasm", # poly1 borderline low hitcall
                    "EPAPLT0202E12-Mito_Intensity_Ring", # poly1 borderline 0.5 hitcall
                    "EPAPLT0202A05-AGP_Compactness_Cells", # hill high hit-cal
                    "EPAPLT0202A11-AGP_Texture_Ring", # exp4 low hit-call
                    "EPAPLT0202F18-RNA_Intensity_Nuclei", # exp5 high hit-call
                    "EPAPLT0202D10-RNA_Intensity_Nuclei", # csnt low hit-call
                    "EPAPLT0202A12-Shape", # poly1 high hit-call
                    "DEX-Shape", # exp4 low hit-call
                    "EPAPLT0202D18-DNA_Symmetry_Nuclei", # exp5 low hit-call
                    "EPAPLT0202E05-Mito_Compactness_Cells" # exp2 high hit-call
                    )

## create the subset, n = 12
htpp_cat_subset <- data.frame(matrix(ncol = ncol(tcpl_cat), nrow = 0))

for (each in cat_chem_id) {
  chem <- strsplit(each, split="-")[[1]][1]
  cat <- strsplit(each, split="-")[[1]][2]
  row <- tcpl_cat[tcpl_cat$chem_id == chem & tcpl_cat$endpoint == cat,]
  htpp_cat_subset <- rbind(htpp_cat_subset, row)
}


## chemicals selected for feature subset
## There are 1300 unique feature. The strategy is just to consider covering
## all types of curves and have a good mix of robust and borderline cases.
## Create three subsets from the output data: the very active ones and the
## inactive ones, and the ones on borderline.
## Select one sample from each type of curves available in each subset
## and then check if that covers all the curves used to fit.
## When selecting samples from subset 3, keep a balance between number of borderline with
## high hit-call and low hit=call.

# temp <- tcpl_feature[tcpl_feature$stype == "test sample"]
# sub1 <- temp %>% filter(hitcall < 0.1)
# table(sub1$fit_method) - sub1 has cnst, exp2, exp4, exp5, hill, poly1 and pow
# sub2 <- temp %>% filter(hitcall > 0.9)
# table(sub2$fit_method) - sub2 has exp2, exp3, exp4, exp5, hill, poly1, poly2 and pow
# sub3 <- temp %>% filter(top_over_cutoff > 0.9 & temp$top_over_cutoff < 1.2)
# table(sub3$fit_method) - sub3 has exp2, exp4, exp5, hill, poly1 and pow

feature_chem_id <- c(
  ## from sub1
  "EPAPLT0202H02-f_1073", # f_1073 cnst inactive
  "EPAPLT0202H02-f_984", # f_984 exp4 inactive
  "EPAPLT0202A14-f_703",  # f_703 exp5 inactive
  "EPAPLT0202E12-f_441",  # f_441 hill inactive
  "EPAPLT0202D10-f_902",  # f_902 pow inactive
  "EPAPLT0202A01-f_841", # f_841 poly1 inactive
  ## from sub2
  "EPAPLT0202A05-f_1105", # f_1105 exp2 active
  "EPAPLT0202E12-f_329",  # f_329 exp3 active
  "EPAPLT0202A05-f_993", # f_993 exp4 active
  "EPAPLT0202F18-f_436", # f_436 exp5 active
  "EPAPLT0202E12-f_556", # f_556 poly1 active
  "EPAPLT0202A05-f_1160", # f_1160 poly2 active
  "EPAPLT0202A05-f_519", # f_519 pow active
  ## from sub3
  "EPAPLT0202E14-f_192", # only hill in sub3, high hit-call
  "EPAPLT0202A01-f_250", # exp2, low hit-call (exp2 all have very low hit-call
  "EPAPLT0202A06-f_1178", # poly1, hit-call between 0.25-0.3
  "EPAPLT0202A01-f_958", # pow, hit-call around 0.5
  "EPAPLT0202A11-f_722" # poly1, hit-call around 0.79
)

## create the subset, n = 18
htpp_feature_subset <- data.frame(matrix(ncol =  ncol(tcpl_feature), nrow = 0))

for (each in feature_chem_id) {
  chem <- strsplit(each, split="-")[[1]][1]
  fname <- strsplit(each, split="-")[[1]][2]
  row <- tcpl_feature[tcpl_feature$chem_id == chem & tcpl_feature$endpoint == fname,]
  htpp_feature_subset <- rbind(htpp_feature_subset, row)
}

# usethis::use_data(htpp_global_subset, htpp_cat_subset, htpp_feature_subset, internal = TRUE)
# utils::sessionInfo()

