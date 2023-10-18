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
                    "STAURO", # none - inactive
                    "EPAPLT0202D18", # poly1 - inactive
                    "EPAPLT0202F22", # poly1 - inactive
                    "EPAPLT0202C01" # csnt - inactive
)

## create the subset, n = 10
htpp_global_subset <- tcpl_global[tcpl_global$chem_id %in% global_chem_id]


## chemicals selected for category subset
## There are 49 unique categories, most of them start with prefixes:
## AGP, DNA, ER, Mito, and RNA, and there's one called Shape.
## Group categories by the prefixes and select two responses,
## one active and one inactive, from each group, and select two from Shape.
## exp3 is not available in the corresponding output data
cat_chem_id <- c(   "EPAPLT0202A05-ER_Radial_Cells", # only poly2 from ER_Radial_Cells, active
                    "EPAPLT0202G01-ER_Intensity_Cytoplasm", # exp4 from ER_Intensity_Cytoplasm inactive
                    "EPAPLT0202A05-AGP_Intensity_Ring", # pow from AGP_Intensity_Ring active
                    "EPAPLT0202A11-AGP_Texture_Ring", # exp4 from AGP_Texture_Ring inactive
                    "EPAPLT0202D18-DNA_Symmetry_Nuclei", # exp5 from DNA_Symmetry_Nuclei inactive
                    "EPAPLT0202A05-DNA_Profile_Nuclei", # hill from DNA_Profile_Nuclei active
                    "EPAPLT0202A06-Mito_Texture_Ring", # poly1 from Mito_Texture_Ring inactive
                    "EPAPLT0202E05-Mito_Compactness_Cells", # exp2 from Mito_Compactness_Cells active
                    "EPAPLT0202A12-Shape", #from Shape active poly1
                    "EPAPLT0202E15-Shape", # inactive poly1 from shape
                    "EPAPLT0202D10-RNA_Intensity_Nuclei", # csnt from RNA_Intensity_Nuclei inactive
                    "EPAPLT0202F18-RNA_Intensity_Nuclei" # exp5 from RNA_Axial_Nuclei active
                    )

## create the subset, n = 12
htpp_cat_subset <- data.frame(matrix(ncol = ncol(tcpl_cat), nrow = 0))

for (each in cat_chem_id) {
  chem <- strsplit(each, split="-")[[1]][1]
  cat <- strsplit(each, split="-")[[1]][2]
  row <- tcpl_cat[tcpl_cat$chem_id == chem & tcpl_cat$endpoint == cat]
  htpp_cat_subset <- rbind(htpp_cat_subset, row)
}


## chemicals selected for feature subset
## There are 1300 unique feature, just consider to include at least one response for each
## curve. Create two subsets from the output data: the very active ones and the
## inactive ones. Select one response from each type of curves available in each subset
## and that covers all the models used to fit.

#temp <- tcpl_feature[tcpl_feature$stype == "test sample"]
#sub1 <- temp[temp$hitcall < 0.1]
#sub2 <- temp[temp$hitcall > 0.9]

feature_chem_id <- c(
  ## from sub1
  "EPAPLT0202H02-f_1073", # f_1073 cnst inactive
  "EPAPLT0202D20-f_380", # f_380 exp2 inactive
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
  "EPAPLT0202F18-f_608", # f_608 hill active
  "EPAPLT0202E12-f_556", # f_556 poly1 active
  "EPAPLT0202A05-f_1160", # f_1160 poly2 active
  "EPAPLT0202A05-f_519" # f_519 pow active
)

## create the subset, n = 15
htpp_feature_subset <- data.frame(matrix(ncol =  ncol(tcpl_feature), nrow = 0))

for (each in feature_chem_id) {
  chem <- strsplit(each, split="-")[[1]][1]
  fname <- strsplit(each, split="-")[[1]][2]
  row <- tcpl_feature[tcpl_feature$chem_id == chem & tcpl_feature$endpoint == fname]
  htpp_feature_subset <- rbind(htpp_feature_subset, row)
}


# usethis::use_data(htpp_global_subset, htpp_cat_subset, htpp_feature_subset, internal = TRUE)
# utils::sessionInfo()

