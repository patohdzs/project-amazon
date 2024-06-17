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

library(sf)
library(tictoc)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter())


# Start timer
tic(msg = "productivity_data.R script", log = TRUE)

# Load municipality-level biomass data
load("data/processed/muni_biomass_2017.Rdata")

# Load spatial municipality-level sample
load("data/processed/spatial_muni_sample.Rdata")

# 2006 AG CENSUS CATTLE FOR SLAUGHTER
load("data/clean/cattle_slaughter_2006.Rdata")

# 2017 AG CENSUS CATTLE SOLD
load("data/clean/cattle_sold_2017.Rdata")

# 2006 AG CENSUS AGRICULTURAL USE AREA
load("data/clean/use_area_2006.Rdata")

# 2017 AG CENSUS AGRICULTURAL USE AREA
load("data/clean/use_area_2017.Rdata")

# LAND COVER AND USE (MAPBIOMAS - MUNI LEVEL)
load("data/clean/land_use_cover_muni.Rdata")

# DEFLATOR (IPA-EP-DI)
load("data/clean/deflator.Rdata")

# HISTORICAL TEMPERATURE
raster_temp <- rast("data/clean/temperature.tif")

# HISTORICAL PRECIPITATION
raster_precip <- rast("data/clean/precipitation.tif")

# Get mean biomass per municipality
muni_biomass_2017 <- muni_biomass_2017 %>%
  group_by(muni_code) %>%
  summarise(agb_2017 = mean(agb_2017, na.rm = TRUE)) %>%
  st_set_geometry(NULL)

# AGGREGATE DEFLATOR BY YEAR
deflator <-
  deflator %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarise(deflator_ipa = mean(deflator_ipa)) %>%
  ungroup() %>%
  filter(year >= 2000)

# CHANGE DEFLATOR BASE TO 2017
aux_deflator_2017 <- deflator[deflator$year == 2017, ]$deflator_ipa
deflator <-
  deflator %>%
  mutate(deflator_ipa = deflator_ipa / aux_deflator_2017)

rm(aux_deflator_2017)



# Commercial exchange rate - selling - average - annual - 2017 - ipeadata
brl_to_usd <- 3.192

# DEFLATE CATTLE SLAUGHTER 2006 AND CONVERT TO USD
# Change from thousand BRL to BRL to USD,
cattle_slaughter_2006 <-
  cattle_slaughter_2006 %>%
  mutate(
    cattleSlaughter_value_2006 = cattleSlaughter2006_value / deflator[deflator$year == 2006, ]$deflator_ipa,
    cattleSlaughter_value_2006 = 1000 * cattleSlaughter_value_2006 / brl_to_usd
  ) %>%
  rename(cattleSlaughter_head_2006 = cattleSlaughter2006_head) %>%
  select(muni_code, cattleSlaughter_value_2006, cattleSlaughter_head_2006)

# Remove deflator
rm(deflator)

# CONVERT CATTLE SOLD VALUE TO USD
# Change from thousand BRL to BRL to USD
cattle_sold_2017 <-
  cattle_sold_2017 %>%
  mutate(cattleSlaughter_value_2017 = 1000 * cattleSoldSlaughterLargeProp_value_2017 / brl_to_usd) %>%
  rename(cattleSlaughter_head_2017 = cattleSoldSlaughterLargeProp_head_2017) %>%
  select(muni_code, cattleSlaughter_value_2017, cattleSlaughter_head_2017)


# Select land use area columns
use_area_2017 <-
  use_area_2017 %>%
  select(
    muni_code,
    agUse_area_2017,
    pasture_area_2017,
    crop_area_2017
  )

use_area_2006 <-
  use_area_2006 %>%
  select(
    muni_code,
    agUse_area_2006,
    pasture_area_2006,
    crop_area_2006
  )


# Match spatial sample crs with raster
spatial_muni_sample <- st_transform(spatial_muni_sample, st_crs(raster_precip))

# Crop rasters
raster_precip <- crop(raster_precip, spatial_muni_sample)
raster_temp <- crop(raster_temp, spatial_muni_sample)

# Calculate average yearly precipitation
raster_precip <- mean(raster_precip)

# Calculate average yearly temperature
raster_temp <- mean(raster_temp)

# Extract total precipitation data by muni
spatial_muni_sample$historical_precip <- extract(raster_precip, vect(spatial_muni_sample), fun = mean, na.rm = TRUE)[, 2]

# Extract average temperature data by muni
spatial_muni_sample$historical_temp <- extract(raster_temp, vect(spatial_muni_sample), fun = mean, na.rm = TRUE)[, 2]

# Reproject spatial sample
spatial_muni_sample <- st_transform(spatial_muni_sample, st_crs(5880))

# Clean environment
rm(raster_precip, raster_temp)


# ADD LON LAT (MUNI CENTROIDS)
aux_centroids <-
  spatial_muni_sample %>%
  st_centroid() %>%
  st_coordinates()

spatial_muni_sample$lon <- aux_centroids[, "X"]
spatial_muni_sample$lat <- aux_centroids[, "Y"]

# Clear environment
rm(aux_centroids)


# MERGE ALL DATASETS AND CREATE VARIABLES OF INTEREST
productivity_data <-
  spatial_muni_sample %>%
  left_join(cattle_sold_2017, by = c("muni_code")) %>%
  left_join(use_area_2017, by = c("muni_code")) %>%
  left_join(cattle_slaughter_2006, by = c("muni_code")) %>%
  left_join(use_area_2006, by = c("muni_code")) %>%
  left_join(muni_biomass_2017, by = c("muni_code")) %>%
  filter(!is.na(biomeAmazon_share)) %>% # remove municipalities outside amazon biome
  # CREATE AGRICULTURAL CENSUS VARIABLES
  mutate(cattleSlaughter_valuePerHa_2017 = if_else(pasture_area_2017 == 0 & !is.na(cattleSlaughter_value_2017), 0, cattleSlaughter_value_2017 / pasture_area_2017)) %>%
  mutate(cattleSlaughter_valuePerHa_2006 = if_else(pasture_area_2006 == 0 & !is.na(cattleSlaughter_value_2006), 0, cattleSlaughter_value_2006 / pasture_area_2006)) %>%
  mutate(
    cattleSlaughter_carcassWeightPerHa_2017 = if_else(pasture_area_2017 == 0, as.numeric(NA), 225 * cattleSlaughter_head_2017 / pasture_area_2017), # average cattle weight 225 kg
    cattleSlaughter_farmGatePrice_2017 = if_else(cattleSlaughter_head_2017 == 0, as.numeric(NA), cattleSlaughter_value_2017 / (cattleSlaughter_head_2017 * 15))
  ) %>% # USD/@ (average cattle weight 15@)
  mutate(
    cattleSlaughter_carcassWeightPerHa_2006 = if_else(pasture_area_2006 == 0, as.numeric(NA), 225 * cattleSlaughter_head_2006 / pasture_area_2006), # average cattle weight 225 kg
    cattleSlaughter_farmGatePrice_2006 = if_else(cattleSlaughter_head_2006 == 0, as.numeric(NA), cattleSlaughter_value_2006 / (cattleSlaughter_head_2006 * 15))
  ) %>% # USD/@ (average cattle weight 15@)
  select(
    lon,
    lat,
    muni_code,
    muni_area,
    biomeAmazon_share,
    agUse_area_2017,
    agUse_area_2006,
    cattleSlaughter_valuePerHa_2017,
    cattleSlaughter_carcassWeightPerHa_2017,
    cattleSlaughter_farmGatePrice_2017,
    pasture_area_2017,
    cattleSlaughter_head_2017,
    cattleSlaughter_valuePerHa_2006,
    cattleSlaughter_carcassWeightPerHa_2006,
    cattleSlaughter_farmGatePrice_2006,
    pasture_area_2006,
    cattleSlaughter_head_2006,
    historical_precip,
    historical_temp,
    geometry,
    agb_2017
  )

# Clear environment
rm(
  use_area_2006,
  use_area_2017,
  cattle_sold_2017,
  cattle_slaughter_2006,
  spatial_muni_sample,
)

# Set labels
set_label(productivity_data$lon) <- "longitude of municipality centroid (calculate under EPSG:5880)"
set_label(productivity_data$lat) <- "latitude  of municipality centroid (calculate under EPSG:5880)"
set_label(productivity_data$historical_precip) <- "historical (1970-2000) total annual precipitation (mm)"
set_label(productivity_data$historical_temp) <- "historical (1970-2000) mean annual average temperature (celsius degrees) "
set_label(productivity_data$agUse_area_2017) <- "agricultural use area (ha, 2017 Ag Census)"
set_label(productivity_data$cattleSlaughter_carcassWeightPerHa_2017) <- "total carcass weigth of cattle sold for slaughter per pasture area (kg/ha, 2017 Ag Census)"
set_label(productivity_data$cattleSlaughter_valuePerHa_2017) <- "value of cattle sold for slaughter per pasture area (2017 constantUSD/ha, 2017 Ag Census)"
set_label(productivity_data$cattleSlaughter_farmGatePrice_2017) <- "farm gate price cattle sold for slaughter (2017 constant USD/@, 2017 Ag Census)"
set_label(productivity_data$agUse_area_2006) <- "agricultural use area (ha)"
set_label(productivity_data$cattleSlaughter_carcassWeightPerHa_2006) <- "total carcass weigth of cattle sold for slaughter per pasture area (kg/ha, 2017 Ag Census)"
set_label(productivity_data$cattleSlaughter_valuePerHa_2006) <- "value of cattle sold for slaughter per pasture area (2017 constantUSD/ha, 2017 Ag Census)"
set_label(productivity_data$cattleSlaughter_farmGatePrice_2006) <- "farm gate price cattle sold for slaughter (2017 constant USD/@, 2017 Ag Census)"


# Save data set
save(productivity_data, zfile = "data/prepData/productivity_data.Rdata")


# End timer
toc(log = TRUE)
