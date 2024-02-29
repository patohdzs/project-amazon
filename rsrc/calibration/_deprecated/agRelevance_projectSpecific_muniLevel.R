
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: COMBINE AGRICULTURAL RELEVANCE VARIABLES (AG CENSUS + PAM + PPM + MAPBIOMAS) WITH SAMPLE OF INTEREST
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "agRelevance_projectSpecific_muniLevel script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(sf) # manipulate spatial data
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns




# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# PANEL SAMPLE
load("data/projectSpecific/muniLevel/sampleSpatial_muniLevel.Rdata")



# AGRICULTURAL CENSUS - SELECTED CROPS AND CATTLE
load(file.path("data/raw2clean/agCensus5crops_ibge/output", "clean_agCensus5crops.Rdata"))
load(file.path("data/raw2clean/agCensusCattle_ibge/output", "clean_agCensusCattle.Rdata"))



# MUNICIPAL AGRICULTURAL PRODUCTION (PAM)
load(file.path("data/raw2clean/pam5crops_ibge/output", "clean_pam5crops.Rdata"))
load(file.path("data/raw2clean/pamTempCrops_ibge/output", "clean_pamTempCrops.Rdata"))

# MUNICIPAL LIVESTOCK SURVEY (PPM)
load(file.path("data/raw2clean/ppmCattle_ibge/output", "clean_ppmCattle.Rdata"))



# LAND COVER AND USE (MAPBIOMAS - MUNI LEVEL)
load(file.path("data/raw2clean/landUseCoverMuni_mapbiomas/output", "clean_landUseCoverMuni.Rdata"))





# DATA PREP ------------------------------------------------------------------------------------------------------------------------------------------

# AGRICULTURAL CENSUS
clean.agCensus5crops <-
  clean.agCensus5crops %>%
  dplyr::select(-muni_name) %>% # remove useless column
  dplyr::group_by(muni_code) %>% # set muni as reference group
  dplyr::mutate(harvested_area_3crops = sum(harvested_area_rice, harvested_area_corn, harvested_area_cassava),
                harvested_area_otherTempCrops = sum(harvested_area_tempCrops, -harvested_area_soybean, -harvested_area_sugarcane)) %>%
  dplyr::rename_at(.vars = vars(dplyr::starts_with("harvested")), .funs = ~ paste0(., "_agCensus")) # rename columns by pattern

clean.agCensusCattle <-
  clean.agCensusCattle %>%
  dplyr::select(-muni_name, -livestock_all) %>% # remove useless columns
  dplyr::rename(livestock_cattle_agCensus = livestock_cattle) # rename column



# PAM 5 CROPS
# column and year selection
clean.pam5crops <-
  clean.pam5crops %>%
  dplyr::filter(year == 2017) %>% # select year of interest
  dplyr::left_join(clean.pamTempCrops, by = c("muni_code", "year")) %>% # tempCrops variables has a lot of NAs because its sample was downloaded only for Legal Amazon municipalities
  dplyr::select(muni_code, dplyr::starts_with("harvested"), dplyr::starts_with("planted")) %>% # select columns of interest
  dplyr::group_by(muni_code) %>% # set muni as reference group
  dplyr::mutate(harvested_area_3crops = sum(harvested_area_rice, harvested_area_corn, harvested_area_cassava),
                harvested_area_otherTempCrops = sum(harvested_area_tempCrops, -harvested_area_soybean, -harvested_area_sugarcane),
                planted_area_3crops = sum(planted_area_rice, planted_area_corn, planted_area_cassava),
                planted_area_otherTempCrops = sum(planted_area_tempCrops, -planted_area_soybean, -planted_area_sugarcane)) %>%
  dplyr::rename_at(.vars = vars(dplyr::starts_with("harvested")), .funs = ~ paste0(., "_pam")) %>% # rename columns by pattern
  dplyr::rename_at(.vars = vars(dplyr::starts_with("planted")), .funs = ~ paste0(., "_pam")) # rename columns by pattern


# PPM CATTLE
# column and year selection
clean.ppmCattle <-
  clean.ppmCattle %>%
  dplyr::filter(year == 2017) %>%
  dplyr::rename(livestock_cattle_ppm = livestock_cattle) %>%
  dplyr::select(-muni_name, -year)



# MAPBIOMAS LAND USE AND COVER

# column and year selection + columns aggregation
clean.landUseCoverMuni <-
  clean.landUseCoverMuni %>%
  dplyr::filter(year == 2017) %>%
  dplyr::mutate(pasture_area_mapbiomas                  = mapbiomasLandCoverId_15,
                soybean_area_mapbiomas                  = mapbiomasLandCoverId_39,
                sugarcane_area_mapbiomas                = mapbiomasLandCoverId_20,
                otherTempCrops_area_mapbiomas           = mapbiomasLandCoverId_41,
                nonObserved_area_mapbiomas              = mapbiomasLandCoverId_27) %>%
  dplyr::select(-mapbiomasLandCoverId_15, -mapbiomasLandCoverId_39, -mapbiomasLandCoverId_20, -mapbiomasLandCoverId_41, -mapbiomasLandCoverId_36,
                -mapbiomasLandCoverId_21, -mapbiomasLandCoverId_27, -mapbiomasLandCoverId_36, -mapbiomasLandCoverId_21) %>%
  dplyr::rowwise(muni_code) %>%
  dplyr::mutate(nonFarming_area_mapbiomas = sum(c_across(starts_with("mapbiomas")))) %>% # sum the value of all the rest of mapbiomas categories not selected above
  dplyr::select(-year, -starts_with("mapbiomasLandCover")) %>%
  dplyr::ungroup()



# SPATIAL SAMPLE
# column selection
sampleSpatial.muniLevel <-
  sampleSpatial.muniLevel %>%
  dplyr::select(muni_code, muni_area, biomeAmazon_share)





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# MERGE ALL DATASETS WITH THE SPATIAL SAMPLE
spatial.agRelevance.muniLevel <-
  sampleSpatial.muniLevel %>%
  dplyr::left_join(clean.agCensus5crops,   by = c("muni_code")) %>%
  dplyr::left_join(clean.pam5crops,        by = c("muni_code")) %>%
  dplyr::left_join(clean.agCensusCattle,   by = c("muni_code")) %>%
  dplyr::left_join(clean.ppmCattle,        by = c("muni_code")) %>%
  dplyr::left_join(clean.landUseCoverMuni, by = c("muni_code"))


# clear environment
rm(clean.agCensus5crops, clean.agCensusCattle, clean.pam5crops, clean.pamTempCrops, clean.ppmCattle, sampleSpatial.muniLevel, clean.landUseCoverMuni)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(spatial.agRelevance.muniLevel$muni_code)                              <- "municipality code (7-digit, IBGE - 2015)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$muni_area)                              <- "municipality area (ha, calculated from shapefile under SIRGAS2000 Polyconic projection)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$biomeAmazon_share)                      <- "share of the municipality area in the Amazon biome"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_rice_agCensus)           <- "harvested area - rice (ha - Agricultural Census 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_sugarcane_agCensus)      <- "harvested area - sugarcane (ha - Agricultural Census 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_cassava_agCensus)        <- "harvested area - cassava (ha - Agricultural Census 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_corn_agCensus)           <- "harvested area - corn (ha - Agricultural Census 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_soybean_agCensus)        <- "harvested area - soybean (ha - Agricultural Census 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_tempCrops_agCensus)      <- "harvested area - total of all temporary crops (ha - Agricultural Census 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_3crops_agCensus)         <- "harvested area - total of corn, cassava, and rice (ha - Agricultural Census 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_otherTempCrops_agCensus) <- "harvested area - total of temporary crops, except soybean and sugarcane (ha - Agricultural Census 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_rice_pam)                <- "harvested area - rice (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_sugarcane_pam)           <- "harvested area - sugarcane (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_cassava_pam)             <- "harvested area - cassava (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_corn_pam)                <- "harvested area - corn (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_soybean_pam)             <- "harvested area - soybean (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_tempCrops_pam)           <- "harvested area - total of all temporary crops (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_3crops_pam)              <- "harvested area - total of corn, cassava, and rice (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$harvested_area_otherTempCrops_pam)      <- "harvested area - total of temporary crops, except soybean and sugarcane (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$planted_area_rice_pam)                  <- "planted area - rice (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$planted_area_sugarcane_pam)             <- "planted area - sugarcane (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$planted_area_cassava_pam)               <- "planted area - cassava (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$planted_area_corn_pam)                  <- "planted area - corn (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$planted_area_soybean_pam)               <- "planted area - soybean (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$planted_area_tempCrops_pam)             <- "planted area - total of all temporary crops (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$planted_area_3crops_pam)                <- "planted area - total of corn, cassava, and rice (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$planted_area_otherTempCrops_pam)        <- "planted area - total of temporary crops, except soybean and sugarcane (ha - PAM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$livestock_cattle_agCensus)              <- "number of cattle (head count - Agricultural Census)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$livestock_cattle_ppm)                   <- "number of cattle (head count - PPM 2017)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$pasture_area_mapbiomas)                 <- "land use area - pasture (ha - MapBiomas Collection 5)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$soybean_area_mapbiomas)                 <- "land use area - soybean (ha - MapBiomas Collection 5)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$sugarcane_area_mapbiomas)               <- "land use area - sugarcane (ha - MapBiomas Collection 5)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$otherTempCrops_area_mapbiomas)          <- "land use area - total of temporary crops, except soybean and sugarcane (ha - MapBiomas Collection 5)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$nonObserved_area_mapbiomas)             <- "land use area - non observed (ha - MapBiomas Collection 5)"
sjlabelled::set_label(spatial.agRelevance.muniLevel$nonFarming_area_mapbiomas)              <- "land use area - non farming (ha - MapBiomas Collection 5)"






# POST-TREATMENT OVERVIEW
# summary(spatial.agRelevance.muniLevel)
# View(spatial.agRelevance.muniLevel)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(spatial.agRelevance.muniLevel,
     file = file.path("data/projectSpecific/muniLevel",
                      paste0("spatial_agRelevance_muniLevel", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/muniLevel")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
