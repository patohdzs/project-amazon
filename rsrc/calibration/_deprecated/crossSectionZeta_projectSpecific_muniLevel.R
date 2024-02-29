
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: COMBINE VARIABLES RELEVANT FOR ZETA CALIBRATION
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "crossSectionZeta_projectSpecific_muniLevel script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns
library(sf) # manipulate spatial data
library(raster) # manipulate spatial data




# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# SAMPLE SPATIAL MUNI LEVEL
load(file.path("data/projectSpecific/muniLevel/sampleSpatial_muniLevel.Rdata"))

# AG CENSUS AGRICULTURAL USE AREA
load(file.path("data/raw2clean/agCensusAgUseArea_ibge/output", "clean_agCensusAgUseArea.Rdata"))

# AG CENSUS WAGE EXPENSE
load(file.path("data/raw2clean/agCensusWage_ibge/output", "clean_agCensusWage.Rdata"))

# HISTORICAL TEMPERATURE
raster.temp <- raster::stack("data/raw2clean/temperature_worldClim/output/clean_temperature.tif")

# HISTORICAL PRECIPITATION
raster.precip <- raster::stack("data/raw2clean/precipitation_worldClim/output/clean_precipitation.tif")





# DATA PREP ------------------------------------------------------------------------------------------------------------------------------------------

# HISTORICAL CLIMATE

# match spatial sample crs with raster
sampleSpatial.muniLevel <- sf::st_transform(sampleSpatial.muniLevel, sf::st_crs(raster.precip))

# crop rasters
raster.precip <- raster::crop(raster.precip, sampleSpatial.muniLevel)
raster.temp <- raster::crop(raster.temp, sampleSpatial.muniLevel)

# calculate total yearly precipitation
raster.precip <- sum(raster.precip)

# calculate average yearly temperature
raster.temp <- mean(raster.temp)

# extract total precipitation data by muni
sampleSpatial.muniLevel$historical_precip <- raster::extract(raster.precip, sampleSpatial.muniLevel, fun = sum, na.rm = T)[,1]

# extract average temperature data by muni
sampleSpatial.muniLevel$historical_temp <- raster::extract(raster.temp, sampleSpatial.muniLevel, fun = mean, na.rm = T)[,1]

# reproject spatial sample
sampleSpatial.muniLevel <- sf::st_transform(sampleSpatial.muniLevel, sf::st_crs(5880))

# clean environment
rm(raster.precip, raster.temp)




# ADD LON LAT (MUNI CENTROIDS)
aux.centroids <- sampleSpatial.muniLevel %>% sf::st_centroid() %>% sf::st_coordinates()
sampleSpatial.muniLevel$lon <- aux.centroids[,"X"]
sampleSpatial.muniLevel$lat <- aux.centroids[,"Y"]

# clear environment
rm(aux.centroids)





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# MERGE ALL DATASETS
crossSection.zeta <-
  sampleSpatial.muniLevel %>%
  dplyr::left_join(clean.agCensusWage, by = c("muni_code")) %>%
  dplyr::left_join(clean.agCensusAgUseArea, by = c("muni_code")) %>%
  dplyr::filter(!is.na(biomeAmazon_share)) %>% # remove municipalities outside amazon biome
  dplyr::mutate(wage_expense = 1000*wage_expense/3.192) %>%  # change from thousand BRL to BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
  dplyr::mutate(wage_ha = dplyr::if_else(agUseArea_value == 0 & !is.na(wage_expense), 0, wage_expense/agUseArea_value)) %>%
  dplyr::select(muni_code, muni_area, biomeAmazon_share, wage_expense, wage_ha, agUseArea_value, pastureArea_value,
                lon, lat, historical_precip, historical_temp, geometry)


# clear environment
rm(clean.agCensusWage, clean.agCensusAgUseArea, sampleSpatial.muniLevel)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(crossSection.zeta$lon) <- "longitude of municipality centroid (calculate under EPSG:5880)"
sjlabelled::set_label(crossSection.zeta$lat) <- "latitude  of municipality centroid (calculate under EPSG:5880)"
sjlabelled::set_label(crossSection.zeta$historical_precip) <- "historical (1970-2000) total annual precipitation (mm)"
sjlabelled::set_label(crossSection.zeta$historical_temp) <- "historical (1970-2000) mean annual average temperature (celsius degrees) "
sjlabelled::set_label(crossSection.zeta$wage_expense) <- "wage expense (USD)"
sjlabelled::set_label(crossSection.zeta$wage_ha) <- "wage expense divided by agricultural use area (USD/ha)"



# POST-TREATMENT OVERVIEW
# summary(crossSection.zeta)
# View(crossSection.zeta)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(crossSection.zeta,
     file = file.path("data/projectSpecific/muniLevel",
                      paste0("crossSection_zeta_muniLevel", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/muniLevel")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
