
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




# START TIMER
tictoc::tic(msg = "sampleConstructionMuni_prepData.R script", log = TRUE)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# BRAZILIAN MUNICIPALITIES DIVISION 2015 SHAPEFILE
load("data/clean/muni_division_2015.Rdata")



# AMAZON BIOME BOUNDARY
load("data/clean/amazon_biome.Rdata")





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# calculate municipality area in square kilometers
muni_division_2015$muni_area <-
  sf::st_area(muni_division_2015) %>%
  units::set_units(ha) %>%
  unclass()

# select columns of interest
muni_division_2015 <-
  muni_division_2015 %>%
  dplyr::select(muni_code, muni_name, state_uf, muni_area, geometry)

amazon_biome <-
  amazon_biome %>%
  dplyr::select(biome_name, geometry)


# combine municipalities with biomes
aux_biomeMuni <- sf::st_intersection(muni_division_2015, amazon_biome)

# clear environment
rm(amazon_biome)

# calculate biome areas inside each municipality
aux_biomeMuni$biome_area <-
  sf::st_area(aux_biomeMuni) %>%
  units::set_units(ha) %>%
  unclass()

# select municipalities of interest - in Amazon or Cerrado
aux_biomeMuni <-
  aux_biomeMuni %>%
  sf::st_drop_geometry() %>% # drop spatial dimension
  tidyr::pivot_wider(names_from = biome_name, values_from = biome_area, values_fill = 0) %>% # convert from long to wide - biome areas in columns
  dplyr::filter(AMAZON > 0) %>% # select municipalities of interest
  dplyr::mutate(biomeAmazon_share  = AMAZON/muni_area) %>%
  dplyr::select(-AMAZON)

# merge original municipalities shapefile with biomes information
sampleMuniSpatial_prepData <-
  muni_division_2015 %>%
  dplyr::left_join(aux_biomeMuni) %>%
  dplyr::filter(biomeAmazon_share > 0)

# clear environment
rm(aux_biomeMuni)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(sampleMuniSpatial_prepData$muni_code)          <- "municipality code (7-digit, IBGE - 2015)"
sjlabelled::set_label(sampleMuniSpatial_prepData$muni_name)          <- "municipality name"
sjlabelled::set_label(sampleMuniSpatial_prepData$state_uf)           <- "state name (abbreviation)"
sjlabelled::set_label(sampleMuniSpatial_prepData$muni_area)          <- "municipality area (ha, calculated from shapefile under SIRGAS2000 Polyconic projection)"
sjlabelled::set_label(sampleMuniSpatial_prepData$biomeAmazon_share)  <- "share of the municipality area in the Amazon biome"




# OTHER EXPORT FORMATS
# extract data.frame from sampleSpatial data
sampleMuniCrossSection_prepData <-
  sampleMuniSpatial_prepData %>%
  sf::st_drop_geometry()

# create panel data from sampleCrossSection and add missing label
sampleMuniPanel_prepData <- tidyr::expand_grid(sampleMuniCrossSection_prepData, year = 2000:2019)
sjlabelled::set_label(sampleMuniPanel_prepData$year)     <- "year of reference (calendar or PRODES year)"



# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(sampleMuniSpatial_prepData,
     file = file.path("data/prepData",
                      "sampleMuniSpatial_prepData.Rdata"))

save(sampleMuniCrossSection_prepData,
     file = file.path("data/prepData",
                      "sampleMuniCrossSection_prepData.Rdata"))

save(sampleMuniPanel_prepData,
     file = file.path("data/prepData",
                      "sampleMuniPanel_prepData.Rdata"))



# END TIMER
tictoc::toc(log = TRUE)



# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------