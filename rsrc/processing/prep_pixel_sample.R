# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: EXTRACT RANDOM SAMPLE OF MAPBIOMAS 30M-PIXELS AND RECOVER FULL PANEL (1985-2019)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(sf)
library(units)
library(tictoc)
library(terra)
library(tidyverse)
library(sjlabelled)
library(conflicted)

conflicts_prefer(terra::extract)

# Start timer
tic(msg = "prep_pixel_sample.R script", log = TRUE)

# Load mapbiomas data
clean_mapbiomas <- rast("data/clean/land_use_cover_2000.tif")

# Set random seed to guarantee reproducibility
set.seed(123)

# Extract random sample cells and keep only coordinate info
pixel_sample <- data.frame(spatSample(clean_mapbiomas, 1200000, na.rm = T, xy = T)[, 1:2])

# Read all MapBiomas layers (one for each year)
aux_mapbiomas <- rast(
  list.files("data/raw/mapbiomas/land_use_cover/",
    pattern = "COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-\\d{4}",
    full.names = TRUE
  )
)

# Change raster stack layer names
names(aux_mapbiomas) <- c(1985:2019)

# Use points sample to extract complete panel
pixel_sample <- cbind(pixel_sample, extract(aux_mapbiomas, pixel_sample, df = T)[, -1])

# Reshape
pixel_sample <- pivot_longer(
  pixel_sample,
  cols = matches("\\d{4}"),
  names_to = "year",
  values_to = "mapbiomas_class"
)

# Minor fix year column
pixel_sample$year <- pixel_sample$year %>%
  str_extract("\\d{4}") %>%
  as.numeric()

# Rename columns
pixel_sample <-
  pixel_sample %>%
  rename(lon = x, lat = y)

# Set labels
set_label(pixel_sample$lon) <- "longitude of the pixel centroid (degrees)"
set_label(pixel_sample$lat) <- "latitude of the pixel centroid (degrees)"
set_label(pixel_sample$year) <- "year of reference"
set_label(pixel_sample$mapbiomas_class) <- "mapbiomas land use/cover classification (id)"

# Save data set
save(pixel_sample, file = "data/processed/pixel_sample.Rdata")

# End timer
toc(log = TRUE)
