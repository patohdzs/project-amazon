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
library(terra)
library(nngeo)
library(sjlabelled)
library(conflicted)

conflicts_prefer(dplyr::filter())
conflicts_prefer(terra::extract)

# Start timer
tic(msg = "prep_muni_data.R script", log = TRUE)

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

# Deflate cattle slaughter value and
# change from thousand-BRL to USD
aux_deflator_2006 <- deflator[deflator$year == 2006, ]$deflator_ipa

cattle_slaughter_2006 <-
  cattle_slaughter_2006 %>%
  mutate(slaughter_value_2006 = cattleSlaughter2006_value / aux_deflator_2006) %>%
  mutate(slaughter_value_2006 = 1000 * slaughter_value_2006 / brl_to_usd) %>%
  rename(slaughter_head_2006 = cattleSlaughter2006_head) %>%
  select(muni_code, slaughter_value_2006, slaughter_head_2006)

# Remove deflator
rm(deflator)

# CONVERT CATTLE SOLD VALUE TO USD
# change from thousand BRL to BRL to USD
cattle_sold_2017 <-
  cattle_sold_2017 %>%
  mutate(slaughter_value_2017 = 1000 * cattleSoldSlaughterLargeProp_value_2017 / brl_to_usd) %>%
  rename(slaughter_head_2017 = cattleSoldSlaughterLargeProp_head_2017) %>%
  select(muni_code, slaughter_value_2017, slaughter_head_2017)


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
spatial_muni_sample <-
  spatial_muni_sample %>%
  st_transform(st_crs(raster_precip))

# Crop rasters
raster_precip <- crop(raster_precip, spatial_muni_sample)
raster_temp <- crop(raster_temp, spatial_muni_sample)

# Calculate average yearly precipitation and temperature
raster_precip <- mean(raster_precip)
raster_temp <- mean(raster_temp)

# Extract total precipitation data by muni
spatial_muni_sample$historical_precip <- extract(raster_precip, vect(spatial_muni_sample), fun = mean, na.rm = TRUE)[, 2]

# Extract average temperature data by muni
spatial_muni_sample$historical_temp <- extract(raster_temp, vect(spatial_muni_sample), fun = mean, na.rm = TRUE)[, 2]

# Reproject spatial sample
spatial_muni_sample <-
  spatial_muni_sample %>%
  st_transform(st_crs(5880))

# Clean environment
rm(raster_precip, raster_temp)


# Add municipality centroid coordinates
aux_centroids <-
  spatial_muni_sample %>%
  st_centroid() %>%
  st_coordinates()

spatial_muni_sample$lon <- aux_centroids[, "X"]
spatial_muni_sample$lat <- aux_centroids[, "Y"]

# Clear environment
rm(aux_centroids)


# Merge all municipal data sets
muni_data <-
  spatial_muni_sample %>%
  left_join(cattle_sold_2017, by = c("muni_code")) %>%
  left_join(use_area_2017, by = c("muni_code")) %>%
  left_join(cattle_slaughter_2006, by = c("muni_code")) %>%
  left_join(use_area_2006, by = c("muni_code")) %>%
  left_join(muni_biomass_2017, by = c("muni_code"))

# Remove municipalities outside Amazon biome
muni_data <- muni_data %>% filter(!is.na(biomeAmazon_share))


# Create agricultural census variables
# NOTE: average cattle weights 225 kg
# NOTE 2: USD/@ (average cattle weight 15@)
muni_data <-
  muni_data %>%
  mutate(
    slaughter_value_per_ha_2017 = if_else(pasture_area_2017 == 0 & !is.na(slaughter_value_2017), 0, slaughter_value_2017 / pasture_area_2017),
    slaughter_value_per_ha_2006 = if_else(pasture_area_2006 == 0 & !is.na(slaughter_value_2006), 0, slaughter_value_2006 / pasture_area_2006),
    carcass_weight_per_ha_2017 = if_else(pasture_area_2017 == 0, as.numeric(NA), 225 * slaughter_head_2017 / pasture_area_2017),
    carcass_weight_per_ha_2006 = if_else(pasture_area_2006 == 0, as.numeric(NA), 225 * slaughter_head_2006 / pasture_area_2006),
    farm_gate_price_2017 = if_else(slaughter_head_2017 == 0, as.numeric(NA), slaughter_value_2017 / (slaughter_head_2017 * 15)),
    farm_gate_price_2006 = if_else(slaughter_head_2006 == 0, as.numeric(NA), slaughter_value_2006 / (slaughter_head_2006 * 15))
  )

# Clear environment
rm(use_area_2006)
rm(use_area_2017)
rm(cattle_sold_2017)
rm(cattle_slaughter_2006)
rm(spatial_muni_sample)

# Set labels
set_label(muni_data$lon) <- "longitude of municipality centroid (calculate under EPSG:5880)"
set_label(muni_data$lat) <- "latitude  of municipality centroid (calculate under EPSG:5880)"
set_label(muni_data$historical_precip) <- "historical (1970-2000) total annual precipitation (mm)"
set_label(muni_data$historical_temp) <- "historical (1970-2000) mean annual average temperature (celsius degrees) "
set_label(muni_data$agUse_area_2017) <- "agricultural use area (ha, 2017 Ag Census)"
set_label(muni_data$carcass_weight_per_ha_2017) <- "total carcass weigth of cattle sold for slaughter per pasture area (kg/ha, 2017 Ag Census)"
set_label(muni_data$slaughter_value_per_ha_2017) <- "value of cattle sold for slaughter per pasture area (2017 constantUSD/ha, 2017 Ag Census)"
set_label(muni_data$farm_gate_price_2017) <- "farm gate price cattle sold for slaughter (2017 constant USD/@, 2017 Ag Census)"
set_label(muni_data$agUse_area_2006) <- "agricultural use area (ha)"
set_label(muni_data$carcass_weight_per_ha_2006) <- "total carcass weigth of cattle sold for slaughter per pasture area (kg/ha, 2017 Ag Census)"
set_label(muni_data$slaughter_value_per_ha_2006) <- "value of cattle sold for slaughter per pasture area (2017 constantUSD/ha, 2017 Ag Census)"
set_label(muni_data$farm_gate_price_2006) <- "farm gate price cattle sold for slaughter (2017 constant USD/@, 2017 Ag Census)"


# Interpolate missing values of farm gate price (alternative approach)

# Select municipality codes for missing observations
missing <- muni_data %>%
  filter(is.na(farm_gate_price_2017)) %>%
  select(muni_code)

# Select prices and weights of complete observations
complete <- muni_data %>%
  filter(!is.na(farm_gate_price_2017)) %>%
  rename(
    price = farm_gate_price_2017,
    weight = slaughter_value_per_ha_2017
  ) %>%
  select(price, weight)

# Compute weighted average price of four nearest neighbours
imputed <- st_join(missing, complete, st_nn, k = 4) %>%
  group_by(muni_code) %>%
  summarize(price = sum(price * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE))

# Assign imputed prices
missing <- is.na(muni_data$farm_gate_price_2017)
muni_data$farm_gate_price_2017[missing] <- imputed$price

# Save data set
save(muni_data, file = "data/processed/muni_data.Rdata")

# End timer
toc(log = TRUE)
