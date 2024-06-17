# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: ADD 2018 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
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
tic(msg = "pixel_biomass_2018.R script", log = TRUE)

# Load pixel data for primary forest
load("data/processed/pixel_primary_forest_2018.Rdata")

# Aboveground biomass rasters
agb_raster <- map(
  list.files(
    "data/raw/esa/above_ground_biomass/",
    pattern = "_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2018-fv3.0.tif",
    full.names = TRUE
  ),
  rast
)

# Transform to sf
pixel_primary_forest <- st_as_sf(
  x = pixel_primary_forest,
  coords = c("lon", "lat"),
  crs = st_crs(4326),
  remove = FALSE
)

# Merge AGB data with primary forest sample
pixel_biomass_2018 <-
  map_df(
    .x = seq_along(agb_raster),
    .f = function(.x) {
      # Transform to spatVector
      aux_polygons <- vect(pixel_primary_forest)

      # Crop spatial points to raster extent
      aux_polygons <- crop(aux_polygons, agb_raster[[.x]])

      # Skip raster with no intersection with the sample
      if (is.null(aux_polygons)) {
        next()
      }

      # Change raster layer name
      names(agb_raster[[.x]]) <- "agb_2018"

      # Extract agb raster data by spatial point
      aux_agb <- extract(agb_raster[[.x]], aux_polygons, xy = TRUE)

      # Adjust and select column names
      aux_agb <- aux_agb %>%
        rename(lon = x, lat = y) %>%
        select(-ID)
    }
  )

# Clear environment
rm(agb_raster)

# Transform to sf
pixel_biomass_2018 <- st_as_sf(
  x = pixel_biomass_2018,
  coords = c("lon", "lat"),
  crs = st_crs(4326),
  remove = FALSE
)

# Check if all points were extracted
if (length(pixel_biomass_2018$lon) != length(pixel_primary_forest$lon)) {
  print("Number of spatial points does not match!")
  print(length(pixel_biomass_2018$lon))
  print(length(pixel_primary_forest$lon))
}

# Set labels
set_label(pixel_biomass_2018$agb_2018) <- "aboveground biomass in 2018 (Mg per ha)"

# Save data set
save(pixel_biomass_2018, file = "data/processed/pixel_biomass_2018.Rdata")

# End timer
toc(log = TRUE)
