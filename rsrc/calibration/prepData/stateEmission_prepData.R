
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: PREPATE DATA TO ESTIMATE PARAMETER K (EMISSION FACTOR OF AGRICULTURAL SECTOR)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "stateEmission_prepData.R script", log = T)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# EMISSIONS DATA
load(here::here("data/raw2clean/emission_seeg/output/clean_emission.Rdata"))


# LAND COVER AND USE (MAPBIOMAS - MUNI LEVEL)
load(here::here("data/raw2clean/landUseCoverMuni_mapbiomas/output", "clean_landUseCoverMuni.Rdata"))


# CROSS-SECTION MUNI-LEVEL SAMPLE
load(here::here("data/calibration/prepData/sampleMuniCrossSection_prepData.Rdata"))





# DATA PREP ------------------------------------------------------------------------------------------------------------------------------------------

# MAPBIOMAS LAND USE AND COVER

# column and year selection + columns aggregation
clean.landUseCoverMuni <-
  clean.landUseCoverMuni %>%
  dplyr::mutate(agriculturalUse_area = mapbiomasLandCoverId_15 + mapbiomasLandCoverId_39 +
                  mapbiomasLandCoverId_41 + mapbiomasLandCoverId_20 + mapbiomasLandCoverId_21) %>% # generate variable of interest
  dplyr::select(muni_code, year, agriculturalUse_area)





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# MERGE ALL DATASETS
stateEmission.prepData <-
  clean.landUseCoverMuni %>%
  dplyr::left_join(sampleMuniCrossSection.prepData, by = c("muni_code")) %>%
  dplyr::left_join(clean.emission, by = c("state_uf", "year")) %>%
  dplyr::filter(year >= 1990) %>% # select period of interest
  dplyr::filter(!is.na(biomeAmazon_share)) %>% # remove municipalities outside amazon biome
  dplyr::group_by(state_uf, year, emission_co2e, netEmission_co2e) %>%
  dplyr::summarise(agriculturalUse_area = sum(agriculturalUse_area),
                   biomeAmazon_share = sum(biomeAmazon_share*muni_area)/sum(muni_area)) %>%
  dplyr::mutate(emission_co2e = emission_co2e*biomeAmazon_share,
                netEmission_co2e = netEmission_co2e*biomeAmazon_share,
                agriculturalUse_area = agriculturalUse_area*biomeAmazon_share) %>%
  dplyr::select(state_uf, year, emission_co2e, netEmission_co2e, agriculturalUse_area)


# clear environment
rm(clean.emission, clean.landUseCoverMuni, sampleMuniCrossSection.prepData)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(stateEmission.prepData$emission_co2e) <- "agricultural use emission factor adjusted by fraction of state area inside amazon biome (CO2e)"
sjlabelled::set_label(stateEmission.prepData$netEmission_co2e) <- "agricultural use net emission factor adjusted by fraction of state area inside amazon biome (CO2e"
sjlabelled::set_label(stateEmission.prepData$agriculturalUse_area) <- "agricultural area adjusted by fraction of state area inside amazon biome (ha)"



# POST-TREATMENT OVERVIEW
# summary(stateEmission.prepData)
# View(stateEmission.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(stateEmission.prepData,
     file = file.path("data/calibration/prepData",
                      paste0("stateEmission_prepData", ".Rdata")))



# # END TIMER
# tictoc::toc(log = T)

# # export time to csv table
# ExportTimeProcessing("code/calibration")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
