
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




# START TIMER
tictoc::tic(msg = "stateEmission_prepData.R script", log = TRUE)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# EMISSIONS DATA
load("data/clean/emission.Rdata")


# LAND COVER AND USE (MAPBIOMAS - MUNI LEVEL)
load("data/clean/land_use_cover_muni.Rdata")


# CROSS-SECTION MUNI-LEVEL SAMPLE
load("data/prepData/sampleMuniCrossSection_prepData.Rdata")



# DATA PREP ------------------------------------------------------------------------------------------------------------------------------------------

# MAPBIOMAS LAND USE AND COVER

# column and year selection + columns aggregation
land_use_cover_muni <-
  land_use_cover_muni %>%
  dplyr::mutate(agriculturalUse_area = mapbiomasLandCoverId_15 + mapbiomasLandCoverId_39 +
                  mapbiomasLandCoverId_41 + mapbiomasLandCoverId_20 + mapbiomasLandCoverId_21) %>% # generate variable of interest
  dplyr::select(muni_code, year, agriculturalUse_area)





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# MERGE ALL DATASETS
stateEmission_prepData <-
  land_use_cover_muni %>%
  dplyr::left_join(sampleMuniCrossSection_prepData, by = c("muni_code")) %>%
  dplyr::left_join(emission, by = c("state_uf", "year")) %>%
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
rm(emission, land_use_cover_muni, sampleMuniCrossSection_prepData)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(stateEmission_prepData$emission_co2e) <- "agricultural use emission factor adjusted by fraction of state area inside amazon biome (CO2e)"
sjlabelled::set_label(stateEmission_prepData$netEmission_co2e) <- "agricultural use net emission factor adjusted by fraction of state area inside amazon biome (CO2e"
sjlabelled::set_label(stateEmission_prepData$agriculturalUse_area) <- "agricultural area adjusted by fraction of state area inside amazon biome (ha)"




# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(stateEmission_prepData,
     file = "data/prepData/stateEmission_prepData.Rdata")



# END TIMER
tictoc::toc(log = TRUE)





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
