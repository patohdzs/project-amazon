
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: DEFINE MUNI-LEVEL SAMPLE AND SAVE IT IN SPATIAL, CROSS-SECTION AND PANEL FORMATS
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "sampleConstructionMuni_prepData.R script", log = T)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# BRAZILIAN MUNICIPALITIES DIVISION 2015 SHAPEFILE
load(here::here("data/raw2clean/muniDivision2015_ibge/output/clean_muniDivision2015.Rdata"))



# AMAZON BIOME BOUNDARY
load(here::here("data/raw2clean/amazonBiome_ibge/output/clean_amazonBiome.Rdata"))





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# calculate municipality area in square kilometers
clean.muniDivision2015$muni_area <-
  sf::st_area(clean.muniDivision2015) %>%
  units::set_units(ha) %>%
  unclass()

# select columns of interest
clean.muniDivision2015 <-
  clean.muniDivision2015 %>%
  dplyr::select(muni_code, muni_name, state_uf, muni_area, geometry)

clean.amazonBiome <-
  clean.amazonBiome %>%
  dplyr::select(biome_name, geometry)


# combine municipalities with biomes
aux.biomeMuni <- sf::st_intersection(clean.muniDivision2015, clean.amazonBiome)

# clear environment
rm(clean.amazonBiome)

# calculate biome areas inside each municipality
aux.biomeMuni$biome_area <-
  sf::st_area(aux.biomeMuni) %>%
  units::set_units(ha) %>%
  unclass()

# select municipalities of interest - in Amazon or Cerrado
aux.biomeMuni <-
  aux.biomeMuni %>%
  sf::st_drop_geometry() %>% # drop spatial dimension
  tidyr::pivot_wider(names_from = biome_name, values_from = biome_area, values_fill = 0) %>% # convert from long to wide - biome areas in columns
  dplyr::filter(AMAZON > 0) %>% # select municipalities of interest
  dplyr::mutate(biomeAmazon_share  = AMAZON/muni_area) %>%
  dplyr::select(-AMAZON)

# merge original municipalities shapefile with biomes information
sampleMuniSpatial.prepData <-
  clean.muniDivision2015 %>%
  dplyr::left_join(aux.biomeMuni) %>%
  dplyr::filter(biomeAmazon_share > 0)

# clear environment
rm(aux.biomeMuni)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(sampleMuniSpatial.prepData$muni_code)          <- "municipality code (7-digit, IBGE - 2015)"
sjlabelled::set_label(sampleMuniSpatial.prepData$muni_name)          <- "municipality name"
sjlabelled::set_label(sampleMuniSpatial.prepData$state_uf)           <- "state name (abbreviation)"
sjlabelled::set_label(sampleMuniSpatial.prepData$muni_area)          <- "municipality area (ha, calculated from shapefile under SIRGAS2000 Polyconic projection)"
sjlabelled::set_label(sampleMuniSpatial.prepData$biomeAmazon_share)  <- "share of the municipality area in the Amazon biome"




# OTHER EXPORT FORMATS
# extract data.frame from sampleSpatial data
sampleMuniCrossSection.prepData <-
  sampleMuniSpatial.prepData %>%
  sf::st_drop_geometry()

# create panel data from sampleCrossSection and add missing label
sampleMuniPanel.prepData <- tidyr::expand_grid(sampleMuniCrossSection.prepData, year = 2000:2019)
sjlabelled::set_label(sampleMuniPanel.prepData$year)     <- "year of reference (calendar or PRODES year)"



# POST-TREATMENT OVERVIEW
# summary(sampleMuniSpatial.prepData)
# View(sampleMuniSpatial.prepData)
# plot(sampleMuniSpatial.prepData$geometry)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(sampleMuniSpatial.prepData,
     file = file.path("data/calibration/prepData",
                      "sampleMuniSpatial_prepData.Rdata"))

save(sampleMuniCrossSection.prepData,
     file = file.path("data/calibration/prepData",
                      "sampleMuniCrossSection_prepData.Rdata"))

save(sampleMuniPanel.prepData,
     file = file.path("data/calibration/prepData",
                      "sampleMuniPanel_prepData.Rdata"))



# # END TIMER
# tictoc::toc(log = T)

# # export time to csv table
# ExportTimeProcessing("code/calibration")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------