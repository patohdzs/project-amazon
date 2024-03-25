
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: COMBINE VARIABLES RELEVANT FOR THETA CALIBRATION AT THE MUNI LEVEL
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -


# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "muniTheta_prepData.R script", log = T)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------


# Gamma MUNI LEVEL
load(here::here("data/calibration/prepData/muni_Biomass2017_prepData.Rdata"))

gamma_merge<-gamma_muni_2017 %>%
  group_by(muni_code)%>%
  summarise(
    agb_2017=mean(agb_2017,na.rm = TRUE)
  )

gamma_merge_df <- st_set_geometry(gamma_merge, NULL)


# SAMPLE SPATIAL MUNI LEVEL
load("data/calibration/prepData/sampleMuniSpatial_prepData.Rdata")


# 2017 AG CENSUS CATTLE SOLD
load("data/raw2clean/agCensus2017CattleSold_ibge/output/clean_agCensus2017CattleSold.Rdata")


# 2017 AG CENSUS AGRICULTURAL USE AREA
load("data/raw2clean/agCensus2017AgUseArea_ibge/output/clean_agCensus2017AgUseArea.Rdata")


# LAND COVER AND USE (MAPBIOMAS - MUNI LEVEL)
load("data/raw2clean/landUseCoverMuni_mapbiomas/output/clean_landUseCoverMuni.Rdata")


# 2006 AG CENSUS CATTLE FOR SLAUGHTER
load("data/raw2clean/agCensus2006CattleSlaughter_ibge/output/clean_agCensus2006CattleSlaughter.Rdata")


# 2006 AG CENSUS AGRICULTURAL USE AREA
load("data/raw2clean/agCensus2006AgUseArea_ibge/output/clean_agCensus2006AgUseArea.Rdata")


# DEFLATOR (IPA-EP-DI)
load("data/raw2clean/deflatorIPA_fgv/output/clean_deflatorIPA.Rdata")



# HISTORICAL TEMPERATURE
raster.temp <- terra::rast("data/raw2clean/temperature_worldClim/output/clean_temperature.tif")


# HISTORICAL PRECIPITATION
raster.precip <- terra::rast("data/raw2clean/precipitation_worldClim/output/clean_precipitation.tif")


# load("data/calibration/z_muni_2017.Rdata")


# z_muni_2017_merge<-z_muni_2017 %>%
#   group_by(muni_code)%>%
#   summarise(
#     agriculturaluse_2017=mean(share_agriculturalUse_2017,na.rm = TRUE),
#     forest_2017=mean(share_forest_2017,na.rm = TRUE),
#     other_2017=mean(share_other_2017,na.rm = TRUE),
#     agriculturaluse_1995=mean(share_agriculturalUse_1995,na.rm = TRUE),
#     forest_1995=mean(share_forest_1995,na.rm = TRUE),
#     other_1995=mean(share_other_1995,na.rm = TRUE),
#     agriculturaluse_2008=mean(share_agriculturalUse_2008,na.rm = TRUE),
#     forest_2008=mean(share_forest_2008,na.rm = TRUE),
#     other_2008=mean(share_other_2008,na.rm = TRUE)
#   )
# z_muni_2017_merge_df <- st_set_geometry(z_muni_2017_merge, NULL)




# DATA PREP ------------------------------------------------------------------------------------------------------------------------------------------

# AGGREGATE DEFLATOR BY YEAR
clean.deflatorIPA <-
  clean.deflatorIPA %>%
  dplyr::mutate(year = lubridate::year(date)) %>% # construct year, month and trimester variables
  dplyr::group_by(year) %>%
  dplyr::summarise(deflator_ipa = mean(deflator_ipa)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(year >= 2000) # select period of interest

# CHANGE DEFLATOR BASE TO 2017
aux.deflator2017 <- clean.deflatorIPA[clean.deflatorIPA$year == 2017,]$deflator_ipa
clean.deflatorIPA <-
  clean.deflatorIPA %>%
  dplyr::mutate(deflator_ipa = deflator_ipa/aux.deflator2017)

rm(aux.deflator2017)


# DEFLATE CATTLE SLAUGHTER 2006 AND CONVERT TO USD
clean.agCensus2006CattleSlaughter <-
  clean.agCensus2006CattleSlaughter %>%
  dplyr::mutate(cattleSlaughter_value_2006 = cattleSlaughter2006_value/clean.deflatorIPA[clean.deflatorIPA$year == 2006,]$deflator_ipa,
                cattleSlaughter_value_2006 = 1000*cattleSlaughter_value_2006/3.192) %>%  # change from thousand BRL to BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
  dplyr::rename(cattleSlaughter_head_2006 = cattleSlaughter2006_head) %>%
  dplyr::select(muni_code, cattleSlaughter_value_2006, cattleSlaughter_head_2006)

# remove object
rm(clean.deflatorIPA)


# CONVERT CATTLE SOLD VALUE TO USD
clean.agCensus2017CattleSold <-
  clean.agCensus2017CattleSold %>%
  dplyr::mutate(cattleSlaughter_value_2017 = 1000*cattleSoldSlaughterLargeProp_value_2017/3.192) %>%  # change from thousand BRL to BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
  dplyr::rename(cattleSlaughter_head_2017 = cattleSoldSlaughterLargeProp_head_2017) %>%
  dplyr::select(muni_code, cattleSlaughter_value_2017, cattleSlaughter_head_2017)


# SELECT VARIABLES
clean.agCensus2017AgUseArea <-
  clean.agCensus2017AgUseArea %>%
  dplyr::select(muni_code, agUse_area_2017, pasture_area_2017, crop_area_2017)

# ADJUST VARIABLE NAMES
clean.agCensus2006AgUseArea <-
  clean.agCensus2006AgUseArea %>%
  dplyr::select(muni_code, agUse_area_2006, pasture_area_2006, crop_area_2006)


# HISTORICAL CLIMATE

# match spatial sample crs with raster
sampleMuniSpatial.prepData <- sf::st_transform(sampleMuniSpatial.prepData, sf::st_crs(raster.precip))

# crop rasters
raster.precip <- terra::crop(raster.precip, sampleMuniSpatial.prepData)
raster.temp <- terra::crop(raster.temp, sampleMuniSpatial.prepData)



# calculate total yearly precipitation
raster.precip <- mean(raster.precip)

# calculate average yearly temperature
raster.temp <- mean(raster.temp)

# extract total precipitation data by muni
sampleMuniSpatial.prepData$historical_precip <- terra::extract(raster.precip, terra::vect(sampleMuniSpatial.prepData), fun = mean, na.rm = T)[,2]


# extract average temperature data by muni
sampleMuniSpatial.prepData$historical_temp <- terra::extract(raster.temp, terra::vect(sampleMuniSpatial.prepData), fun = mean, na.rm = T)[,2]

# reproject spatial sample
sampleMuniSpatial.prepData <- sf::st_transform(sampleMuniSpatial.prepData, sf::st_crs(5880))

# clean environment
rm(raster.precip, raster.temp)



# ADD LON LAT (MUNI CENTROIDS)
aux.centroids <- sampleMuniSpatial.prepData %>% sf::st_centroid() %>% sf::st_coordinates()
sampleMuniSpatial.prepData$lon <- aux.centroids[,"X"]
sampleMuniSpatial.prepData$lat <- aux.centroids[,"Y"]

# clear environment
rm(aux.centroids)





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------


# MERGE ALL DATASETS AND CREATE VARIABLES OF INTEREST
muniTheta.prepData <-
  sampleMuniSpatial.prepData %>%
  dplyr::left_join(clean.agCensus2017CattleSold, by = c("muni_code")) %>%
  dplyr::left_join(clean.agCensus2017AgUseArea, by = c("muni_code")) %>%
  dplyr::left_join(clean.agCensus2006CattleSlaughter, by = c("muni_code")) %>%
  dplyr::left_join(clean.agCensus2006AgUseArea, by = c("muni_code")) %>%
  dplyr::left_join(gamma_merge_df,by=c("muni_code"))%>%
  dplyr::filter(!is.na(biomeAmazon_share)) %>% # remove municipalities outside amazon biome
  # CREATE AGRICULTURAL CENSUS VARIABLES
  dplyr::mutate(cattleSlaughter_valuePerHa_2017 = dplyr::if_else(pasture_area_2017 == 0 & !is.na(cattleSlaughter_value_2017), 0, cattleSlaughter_value_2017/pasture_area_2017)) %>%
  dplyr::mutate(cattleSlaughter_valuePerHa_2006 = dplyr::if_else(pasture_area_2006 == 0 & !is.na(cattleSlaughter_value_2006), 0, cattleSlaughter_value_2006/pasture_area_2006)) %>%
  dplyr::mutate(cattleSlaughter_carcassWeightPerHa_2017 = if_else(pasture_area_2017  == 0, as.numeric(NA), 225*cattleSlaughter_head_2017/pasture_area_2017), # average cattle weight 225 kg
                cattleSlaughter_farmGatePrice_2017 = if_else(cattleSlaughter_head_2017  == 0, as.numeric(NA), cattleSlaughter_value_2017/(cattleSlaughter_head_2017*15))) %>%  # USD/@ (average cattle weight 15@)
  dplyr::mutate(cattleSlaughter_carcassWeightPerHa_2006 = if_else(pasture_area_2006  == 0, as.numeric(NA), 225*cattleSlaughter_head_2006/pasture_area_2006), # average cattle weight 225 kg
                cattleSlaughter_farmGatePrice_2006 = if_else(cattleSlaughter_head_2006  == 0, as.numeric(NA), cattleSlaughter_value_2006/(cattleSlaughter_head_2006*15))) %>%  # USD/@ (average cattle weight 15@)
  # dplyr::mutate(                forestArea_2017_ha_muni = forest_2017*muni_area,
  #                               z_2017_muni = agriculturaluse_2017*muni_area,
  #                               otherArea_2017_ha_muni = other_2017*muni_area,
  #                               zbar_2017_muni = forestArea_2017_ha_muni + z_2017_muni,
  #                               forestArea_1995_ha_muni = forest_1995*muni_area,
  #                               z_1995_muni = agriculturaluse_1995*muni_area,
  #                               otherArea_1995_ha_muni = other_1995*muni_area,
  #                               zbar_1995_muni = forestArea_1995_ha_muni + z_1995_muni,
  #                               forestArea_2008_ha_muni = forest_2008*muni_area,
  #                               z_2008_muni = agriculturaluse_2008*muni_area,
  #                               otherArea_2008_ha_muni = other_2008*muni_area,
  #                               zbar_2008_muni = forestArea_2008_ha_muni + z_2008_muni)%>%
  # SELECT VARIABLES OF INTEREST
  dplyr::select(muni_code, muni_area, biomeAmazon_share,
                agUse_area_2017, agUse_area_2006,
                cattleSlaughter_valuePerHa_2017, cattleSlaughter_carcassWeightPerHa_2017, cattleSlaughter_farmGatePrice_2017, pasture_area_2017, cattleSlaughter_head_2017,
                cattleSlaughter_valuePerHa_2006, cattleSlaughter_carcassWeightPerHa_2006, cattleSlaughter_farmGatePrice_2006, pasture_area_2006, cattleSlaughter_head_2006,
                lon, lat, historical_precip, historical_temp, geometry,agb_2017)
                # forestArea_2017_ha_muni,z_2017_muni,otherArea_2017_ha_muni,zbar_2017_muni,agb_2017,
                # forestArea_2008_ha_muni,z_2008_muni,otherArea_2008_ha_muni,zbar_2008_muni,
                # forestArea_1995_ha_muni,z_1995_muni,otherArea_1995_ha_muni,zbar_1995_muni)


# clear environment
rm(clean.agCensus2017AgUseArea, sampleMuniSpatial.prepData,
   clean.agCensus2017CattleSold,
   clean.agCensus2006AgUseArea, clean.agCensus2006CattleSlaughter)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(muniTheta.prepData$lon) <- "longitude of municipality centroid (calculate under EPSG:5880)"
sjlabelled::set_label(muniTheta.prepData$lat) <- "latitude  of municipality centroid (calculate under EPSG:5880)"
sjlabelled::set_label(muniTheta.prepData$historical_precip) <- "historical (1970-2000) total annual precipitation (mm)"
sjlabelled::set_label(muniTheta.prepData$historical_temp) <- "historical (1970-2000) mean annual average temperature (celsius degrees) "
sjlabelled::set_label(muniTheta.prepData$agUse_area_2017) <- "agricultural use area (ha, 2017 Ag Census)"
sjlabelled::set_label(muniTheta.prepData$cattleSlaughter_carcassWeightPerHa_2017) <- "total carcass weigth of cattle sold for slaughter per pasture area (kg/ha, 2017 Ag Census)"
sjlabelled::set_label(muniTheta.prepData$cattleSlaughter_valuePerHa_2017) <- "value of cattle sold for slaughter per pasture area (2017 constantUSD/ha, 2017 Ag Census)"
sjlabelled::set_label(muniTheta.prepData$cattleSlaughter_farmGatePrice_2017) <- "farm gate price cattle sold for slaughter (2017 constant USD/@, 2017 Ag Census)"
sjlabelled::set_label(muniTheta.prepData$agUse_area_2006) <- "agricultural use area (ha)"
sjlabelled::set_label(muniTheta.prepData$cattleSlaughter_carcassWeightPerHa_2006) <- "total carcass weigth of cattle sold for slaughter per pasture area (kg/ha, 2017 Ag Census)"
sjlabelled::set_label(muniTheta.prepData$cattleSlaughter_valuePerHa_2006) <- "value of cattle sold for slaughter per pasture area (2017 constantUSD/ha, 2017 Ag Census)"
sjlabelled::set_label(muniTheta.prepData$cattleSlaughter_farmGatePrice_2006) <- "farm gate price cattle sold for slaughter (2017 constant USD/@, 2017 Ag Census)"


# POST-TREATMENT OVERVIEW
# summary(muniTheta.prepData)
# View(muniTheta.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(muniTheta.prepData,
     file = "data/calibration/prepData/muniTheta_prepData.Rdata")



# # END TIMER
# tictoc::toc(log = T)

# # export time to csv table
# ExportTimeProcessing("code/calibration")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
