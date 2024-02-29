
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CREATE AGGREGATED CATEGORIES OF INTEREST BASED ON MAPBIOMAS CATEGORIES AT THE MUNI LEVEL
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# GROUNDHOG (REPRODUCIBILITY SOLUTION TO HANDLING DIFFERENT VERSIONS OF R AND ITS PACKAGES)

# check if groundhog is installed and load it
if ("groundhog" %in% installed.packages()) {
  library("groundhog")
} else {
  install.packages("groundhog")
  library("groundhog")
}

# define date of reference to load all packages
groundhog.date <- "2022-04-01"

# guarantee version 1.5 of groundhog is being used
groundhog::meta.groundhog(date = "2022-04-01")


# HERE
groundhog::groundhog.library("here", groundhog.date) # load package here


# TICTOC
groundhog::groundhog.library("tictoc", groundhog.date) # load package tictoc


# DECLARE LOCATION OF CURRENT SCRIPT TO SET UP PROJECT ROOT CORRECTLY
here::i_am("code/projectSpecific/prepData/muniCategories_projectSpecific_prepData.R", uuid = "feeeefe6-67dd-4e8a-8b8e-73d8f33e4a4e")


# START TIME
tictoc::tic(msg = "muniCategories_projectSpecific_prepData script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sf", groundhog.date)  # manipulate spatial data (vector format)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# SPATIAL MUNI-LEVEL SAMPLE
load(here::here("data/projectSpecific/prepData/sampleMuniSpatial_prepData.Rdata"))



# LAND COVER AND USE (MAPBIOMAS - MUNI LEVEL)
load(here::here("data/raw2clean/landUseCoverMuni_mapbiomas/output", "clean_landUseCoverMuni.Rdata"))





# DATA PREP ------------------------------------------------------------------------------------------------------------------------------------------

# add aggregated land use/cover categories
clean.landUseCoverMuni <-
  clean.landUseCoverMuni %>%
  dplyr::mutate(forest          = mapbiomasLandCoverId_3,
                pasture         = mapbiomasLandCoverId_15,
                soybean         = mapbiomasLandCoverId_39,
                otherCrops      = mapbiomasLandCoverId_41 + mapbiomasLandCoverId_20 + mapbiomasLandCoverId_21) %>%
  dplyr::select(-mapbiomasLandCoverId_3, -mapbiomasLandCoverId_15, -mapbiomasLandCoverId_39, -mapbiomasLandCoverId_41, -mapbiomasLandCoverId_20,
                -mapbiomasLandCoverId_21) %>%
  dplyr::rowwise(muni_code) %>%
  dplyr::mutate(otherCategories = sum(c_across(starts_with("mapbiomas")))) %>% # sum the value of all the rest of mapbiomas categories not selected above
  dplyr::select(-starts_with("mapbiomasLandCover")) %>%
  dplyr::ungroup()



# SPATIAL SAMPLE
# column selection
sampleMuniSpatial.prepData <-
  sampleMuniSpatial.prepData %>%
  dplyr::select(muni_code, muni_area, biomeAmazon_share)





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# MERGE ALL DATASETS WITH THE SPATIAL SAMPLE
muniCategories.prepData <-
  sampleMuniSpatial.prepData %>%
  dplyr::left_join(clean.landUseCoverMuni, by = c("muni_code"))


# clear environment
rm(sampleMuniSpatial.prepData, clean.landUseCoverMuni)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
# all variables already labelled


# POST-TREATMENT OVERVIEW
# summary(spatial.mapbiomasCategories)
# View(spatial.mapbiomasCategories)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(muniCategories.prepData,
     file = here::here("data/projectSpecific/prepData",
                      paste0("muniCategories_prepData", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/prepData")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
