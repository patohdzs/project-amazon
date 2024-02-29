
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: COMBINE VARIABLES RELEVANT FOR THETA CALIBRATION
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "crossSectionTheta2006_projectSpecific_muniLevel script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns
library(sf) # manipulate spatial data
library(raster) # manipulate spatial data




# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# LAND COVER AND USE (MAPBIOMAS - MUNI LEVEL)
load(file.path("data/raw2clean/landUseCoverMuni_mapbiomas/output", "clean_landUseCoverMuni.Rdata"))


# SAMPLE SPATIAL MUNI LEVEL
load(file.path("data/projectSpecific/muniLevel/sampleSpatial_muniLevel.Rdata"))


# AG CENSUS 2006 CATTLE FOR SLAUGHTER
load(file.path("data/raw2clean/agCensus2006CattleSlaughter_ibge/output", "clean_agCensus2006CattleSlaughter.Rdata"))


# AG CENSUS 2006 PASTURE AREA
load(file.path("data/raw2clean/agCensus2006PastureArea_ibge/output", "clean_agCensus2006PastureArea.Rdata"))


# AGRICUTLURAL COMMODITY PRICES
load("data/raw2clean/commodityPrices_seabpr/output/clean_commodityPrices.Rdata")


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
crossSection2006.theta <-
  sampleSpatial.muniLevel %>%
  dplyr::left_join(clean.agCensus2006CattleSlaughter, by = c("muni_code")) %>%
  dplyr::left_join(clean.agCensus2006PastureArea, by = c("muni_code")) %>%
  dplyr::left_join(clean.landUseCoverMuni %>% dplyr::filter(year == 2006) %>% dplyr::select(muni_code, mapbiomas2006_forest_ha = mapbiomasLandCoverId_3)) %>%
  dplyr::filter(!is.na(biomeAmazon_share)) %>% # remove municipalities outside amazon biome
  dplyr::mutate(cattleSlaughter2006_value = 1000*cattleSlaughter_value_2006/2.1761) %>%  # change from thousand BRL to BRL to USD (commercial exchange rate - selling - average - annual - 2006 - ipeadata))
  dplyr::mutate(cattleSlaughter2006_value_ha = dplyr::if_else(pastureArea_value_2006 == 0 & !is.na(cattleSlaughter2006_value), 0, cattleSlaughter2006_value/pastureArea_value_2006)) %>%
  dplyr::select(muni_code, muni_area, biomeAmazon_share, starts_with("cattle"), starts_with("pasture"), mapbiomas2006_forest_ha,
                lon, lat, historical_precip, historical_temp, geometry)



# clear environment
rm(clean.landUseCoverMuni, clean.agCensus2006CattleSlaughter, clean.agCensus2006PastureArea, sampleSpatial.muniLevel)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(crossSection2006.theta$lon) <- "longitude of municipality centroid (calculate under EPSG:5880)"
sjlabelled::set_label(crossSection2006.theta$lat) <- "latitude  of municipality centroid (calculate under EPSG:5880)"
sjlabelled::set_label(crossSection2006.theta$historical_precip) <- "historical (1970-2000) total annual precipitation (mm)"
sjlabelled::set_label(crossSection2006.theta$historical_temp) <- "historical (1970-2000) mean annual average temperature (celsius degrees) "
sjlabelled::set_label(crossSection2006.theta$cattleSlaughter2006_value_ha) <- "value of cattle sold for slaughter per pasture area 2006 Agricultural Census (USD/ha)"
sjlabelled::set_label(crossSection2006.theta$mapbiomas2006_forest_ha) <- "forest area in ha in 2006"


# POST-TREATMENT OVERVIEW
# summary(crossSection2006.theta)
# View(crossSection2006.theta)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(crossSection2006.theta,
     file = file.path("data/projectSpecific/muniLevel",
                      paste0("crossSection2006_theta_muniLevel", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/muniLevel")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
