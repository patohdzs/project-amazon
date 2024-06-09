
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




# START TIMER
tictoc::tic(msg = "muniTheta_prepData.R script", log = TRUE)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------


# Gamma MUNI LEVEL
load("data/prepData/muni_Biomass2017_prepData.Rdata")

gamma_merge<-gamma_muni_2017 %>%
  group_by(muni_code)%>%
  summarise(
    agb_2017=mean(agb_2017,na.rm = TRUE)
  )

gamma_merge_df <- st_set_geometry(gamma_merge, NULL)


# SAMPLE SPATIAL MUNI LEVEL
load("data/prepData/sampleMuniSpatial_prepData.Rdata")


# 2017 AG CENSUS CATTLE SOLD
load("data/clean/agcensus2017_cattlesold.Rdata")


# 2017 AG CENSUS AGRICULTURAL USE AREA
load("data/clean/agcensus2017_ag_usearea.Rdata")


# LAND COVER AND USE (MAPBIOMAS - MUNI LEVEL)
load("data/clean/land_use_cover_muni.Rdata")


# 2006 AG CENSUS CATTLE FOR SLAUGHTER
load("data/clean/agcensus2006_cattlesold.Rdata")


# 2006 AG CENSUS AGRICULTURAL USE AREA
load("data/clean/agcensus2006_ag_usearea.Rdata")


# DEFLATOR (IPA-EP-DI)
load("data/clean/deflatorIPA.Rdata")



# HISTORICAL TEMPERATURE
raster_temp <- terra::rast("data/clean/temperature.tif")


# HISTORICAL PRECIPITATION
raster_precip <- terra::rast("data/clean/precipitation.tif")




# DATA PREP ------------------------------------------------------------------------------------------------------------------------------------------

# AGGREGATE DEFLATOR BY YEAR
clean_deflatorIPA <-
  clean_deflatorIPA %>%
  dplyr::mutate(year = lubridate::year(date)) %>% # construct year, month and trimester variables
  dplyr::group_by(year) %>%
  dplyr::summarise(deflator_ipa = mean(deflator_ipa)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(year >= 2000) # select period of interest

# CHANGE DEFLATOR BASE TO 2017
aux_deflator2017 <- clean_deflatorIPA[clean_deflatorIPA$year == 2017,]$deflator_ipa
clean_deflatorIPA <-
  clean_deflatorIPA %>%
  dplyr::mutate(deflator_ipa = deflator_ipa/aux_deflator2017)

rm(aux_deflator2017)


# DEFLATE CATTLE SLAUGHTER 2006 AND CONVERT TO USD
clean_agCensus2006CattleSlaughter <-
  clean_agCensus2006CattleSlaughter %>%
  dplyr::mutate(cattleSlaughter_value_2006 = cattleSlaughter2006_value/clean_deflatorIPA[clean_deflatorIPA$year == 2006,]$deflator_ipa,
                cattleSlaughter_value_2006 = 1000*cattleSlaughter_value_2006/3.192) %>%  # change from thousand BRL to BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
  dplyr::rename(cattleSlaughter_head_2006 = cattleSlaughter2006_head) %>%
  dplyr::select(muni_code, cattleSlaughter_value_2006, cattleSlaughter_head_2006)

# remove object
rm(clean_deflatorIPA)


# CONVERT CATTLE SOLD VALUE TO USD
clean_agCensus2017CattleSold <-
  clean_agCensus2017CattleSold %>%
  dplyr::mutate(cattleSlaughter_value_2017 = 1000*cattleSoldSlaughterLargeProp_value_2017/3.192) %>%  # change from thousand BRL to BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
  dplyr::rename(cattleSlaughter_head_2017 = cattleSoldSlaughterLargeProp_head_2017) %>%
  dplyr::select(muni_code, cattleSlaughter_value_2017, cattleSlaughter_head_2017)


# SELECT VARIABLES
clean_agCensus2017AgUseArea <-
  clean_agCensus2017AgUseArea %>%
  dplyr::select(muni_code, agUse_area_2017, pasture_area_2017, crop_area_2017)

# ADJUST VARIABLE NAMES
clean_agCensus2006AgUseArea <-
  clean_agCensus2006AgUseArea %>%
  dplyr::select(muni_code, agUse_area_2006, pasture_area_2006, crop_area_2006)


# HISTORICAL CLIMATE

# match spatial sample crs with raster
sampleMuniSpatial_prepData <- sf::st_transform(sampleMuniSpatial_prepData, sf::st_crs(raster_precip))

# crop rasters
raster_precip <- terra::crop(raster_precip, sampleMuniSpatial_prepData)
raster_temp <- terra::crop(raster_temp, sampleMuniSpatial_prepData)



# calculate total yearly precipitation
raster_precip <- mean(raster_precip)

# calculate average yearly temperature
raster_temp <- mean(raster_temp)

# extract total precipitation data by muni
sampleMuniSpatial_prepData$historical_precip <- terra::extract(raster_precip, terra::vect(sampleMuniSpatial_prepData), fun = mean, na.rm = T)[,2]


# extract average temperature data by muni
sampleMuniSpatial_prepData$historical_temp <- terra::extract(raster_temp, terra::vect(sampleMuniSpatial_prepData), fun = mean, na.rm = T)[,2]

# reproject spatial sample
sampleMuniSpatial_prepData <- sf::st_transform(sampleMuniSpatial_prepData, sf::st_crs(5880))

# clean environment
rm(raster_precip, raster_temp)



# ADD LON LAT (MUNI CENTROIDS)
aux_centroids <- sampleMuniSpatial_prepData %>% sf::st_centroid() %>% sf::st_coordinates()
sampleMuniSpatial_prepData$lon <- aux_centroids[,"X"]
sampleMuniSpatial_prepData$lat <- aux_centroids[,"Y"]

# clear environment
rm(aux_centroids)





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------


# MERGE ALL DATASETS AND CREATE VARIABLES OF INTEREST
muniTheta_prepData <-
  sampleMuniSpatial_prepData %>%
  dplyr::left_join(clean_agCensus2017CattleSold, by = c("muni_code")) %>%
  dplyr::left_join(clean_agCensus2017AgUseArea, by = c("muni_code")) %>%
  dplyr::left_join(clean_agCensus2006CattleSlaughter, by = c("muni_code")) %>%
  dplyr::left_join(clean_agCensus2006AgUseArea, by = c("muni_code")) %>%
  dplyr::left_join(gamma_merge_df,by=c("muni_code"))%>%
  dplyr::filter(!is.na(biomeAmazon_share)) %>% # remove municipalities outside amazon biome
  # CREATE AGRICULTURAL CENSUS VARIABLES
  dplyr::mutate(cattleSlaughter_valuePerHa_2017 = dplyr::if_else(pasture_area_2017 == 0 & !is.na(cattleSlaughter_value_2017), 0, cattleSlaughter_value_2017/pasture_area_2017)) %>%
  dplyr::mutate(cattleSlaughter_valuePerHa_2006 = dplyr::if_else(pasture_area_2006 == 0 & !is.na(cattleSlaughter_value_2006), 0, cattleSlaughter_value_2006/pasture_area_2006)) %>%
  dplyr::mutate(cattleSlaughter_carcassWeightPerHa_2017 = if_else(pasture_area_2017  == 0, as.numeric(NA), 225*cattleSlaughter_head_2017/pasture_area_2017), # average cattle weight 225 kg
                cattleSlaughter_farmGatePrice_2017 = if_else(cattleSlaughter_head_2017  == 0, as.numeric(NA), cattleSlaughter_value_2017/(cattleSlaughter_head_2017*15))) %>%  # USD/@ (average cattle weight 15@)
  dplyr::mutate(cattleSlaughter_carcassWeightPerHa_2006 = if_else(pasture_area_2006  == 0, as.numeric(NA), 225*cattleSlaughter_head_2006/pasture_area_2006), # average cattle weight 225 kg
                cattleSlaughter_farmGatePrice_2006 = if_else(cattleSlaughter_head_2006  == 0, as.numeric(NA), cattleSlaughter_value_2006/(cattleSlaughter_head_2006*15))) %>%  # USD/@ (average cattle weight 15@)
  dplyr::select(muni_code, muni_area, biomeAmazon_share,
                agUse_area_2017, agUse_area_2006,
                cattleSlaughter_valuePerHa_2017, cattleSlaughter_carcassWeightPerHa_2017, cattleSlaughter_farmGatePrice_2017, pasture_area_2017, cattleSlaughter_head_2017,
                cattleSlaughter_valuePerHa_2006, cattleSlaughter_carcassWeightPerHa_2006, cattleSlaughter_farmGatePrice_2006, pasture_area_2006, cattleSlaughter_head_2006,
                lon, lat, historical_precip, historical_temp, geometry,agb_2017)
                # forestArea_2017_ha_muni,z_2017_muni,otherArea_2017_ha_muni,zbar_2017_muni,agb_2017,
                # forestArea_2008_ha_muni,z_2008_muni,otherArea_2008_ha_muni,zbar_2008_muni,
                # forestArea_1995_ha_muni,z_1995_muni,otherArea_1995_ha_muni,zbar_1995_muni)


# clear environment
rm(clean_agCensus2017AgUseArea, sampleMuniSpatial_prepData,
   clean_agCensus2017CattleSold,
   clean_agCensus2006AgUseArea, clean_agCensus2006CattleSlaughter)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(muniTheta_prepData$lon) <- "longitude of municipality centroid (calculate under EPSG:5880)"
sjlabelled::set_label(muniTheta_prepData$lat) <- "latitude  of municipality centroid (calculate under EPSG:5880)"
sjlabelled::set_label(muniTheta_prepData$historical_precip) <- "historical (1970-2000) total annual precipitation (mm)"
sjlabelled::set_label(muniTheta_prepData$historical_temp) <- "historical (1970-2000) mean annual average temperature (celsius degrees) "
sjlabelled::set_label(muniTheta_prepData$agUse_area_2017) <- "agricultural use area (ha, 2017 Ag Census)"
sjlabelled::set_label(muniTheta_prepData$cattleSlaughter_carcassWeightPerHa_2017) <- "total carcass weigth of cattle sold for slaughter per pasture area (kg/ha, 2017 Ag Census)"
sjlabelled::set_label(muniTheta_prepData$cattleSlaughter_valuePerHa_2017) <- "value of cattle sold for slaughter per pasture area (2017 constantUSD/ha, 2017 Ag Census)"
sjlabelled::set_label(muniTheta_prepData$cattleSlaughter_farmGatePrice_2017) <- "farm gate price cattle sold for slaughter (2017 constant USD/@, 2017 Ag Census)"
sjlabelled::set_label(muniTheta_prepData$agUse_area_2006) <- "agricultural use area (ha)"
sjlabelled::set_label(muniTheta_prepData$cattleSlaughter_carcassWeightPerHa_2006) <- "total carcass weigth of cattle sold for slaughter per pasture area (kg/ha, 2017 Ag Census)"
sjlabelled::set_label(muniTheta_prepData$cattleSlaughter_valuePerHa_2006) <- "value of cattle sold for slaughter per pasture area (2017 constantUSD/ha, 2017 Ag Census)"
sjlabelled::set_label(muniTheta_prepData$cattleSlaughter_farmGatePrice_2006) <- "farm gate price cattle sold for slaughter (2017 constant USD/@, 2017 Ag Census)"




# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(muniTheta_prepData,
     file = "data/prepData/muniTheta_prepData.Rdata")



# END TIMER
tictoc::toc(log = TRUE)





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
