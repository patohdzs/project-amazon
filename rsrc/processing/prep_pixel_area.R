# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CALCULATE AREA OF MAPBIOMAS 30M-PIXELS
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(sf)
library(tictoc)
library(terra)
library(tidyverse)
library(conflicted)
library(sjlabelled)

conflicts_prefer(dplyr::filter())
conflicts_prefer(terra::extract())

# Start timer
tic(msg = "pixel_area.R script", log = TRUE)

# Load raster data
clean_mapbiomas <- rast("data/clean/land_use_cover_2000.tif")

# Load MapBiomas 30m pixel sample
load("data/processed/pixel_sample.Rdata")

# Calculate area of each pixel (sq km)
mapbiomas_area <- cellSize(clean_mapbiomas, unit = "ha")

# Change raster layer name
names(mapbiomas_area) <- "pixel_area"

# Clean environment
rm(clean_mapbiomas)

# Extract one year from sample to extract pixel area
pixel_sample_2000 <-
  pixel_sample %>%
  filter(year == 2000) %>%
  select(-year, -mapbiomas_class) %>%
  st_as_sf(
    coords = c("lon", "lat"),
    remove = FALSE,
    crs = st_crs(4326)
  )

# Extract pixel area raster data for sample points
aux_pixel_areas <- extract(mapbiomas_area, pixel_sample_2000)

# Merge pixel area variable with sample 2000
pixel_sample_2000$pixel_area <- aux_pixel_areas

# Merge sample 2000 with panel sample
pixel_areas <-
  pixel_sample %>%
  left_join(pixel_sample_2000)

# Clear environmnet
rm(pixel_sample, pixel_sample_2000, aux_pixel_areas)

# Set labels
set_label(pixel_areas$pixel_area) <- "pixel area (ha)"

# Save data set
save(pixel_areas, file = "data/processed/pixel_areas.Rdata")

# End timer
toc(log = TRUE)
