# Load necessary libraries
library(sf)
library(tictoc)
library(terra)
library(tidyverse)
library(conflicted)
library(sjlabelled)

# Resolve conflicts
conflicts_prefer(dplyr::filter())
conflicts_prefer(terra::extract())

# Start timer
tic(msg = "merge_agb_sec_veg.R script", log = TRUE)

# Load secondary vegetation age data for 2017
sec_veg_age_rst <- rast(
  list.files(
    "data/raw/mapbiomas/secondary_vegetation_age/",
    pattern = "2017",
    full.names = TRUE
  )
)

# Load aboveground biomass raster data
agb_rasters  <- map(
  list.files(
    "data/raw/esa/above_ground_biomass/",
    pattern = "_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv3.0.tif",
    full.names = TRUE
  ),
  rast
)

# Combine resampled AGB rasters into a single raster stack
combined_agb_raster <- do.call(merge, agb_rasters)

# resample to match sec_veg resolution 
resampled_agb_raster <- resample(combined_agb_raster, sec_veg_age_rst, method = "near")


writeRaster(resampled_agb_raster, "data/calibration/resampled_agb_2017.tif", overwrite = TRUE)


load("data/calibration/calibration_1043_sites.Rdata")
calib_df <- st_transform(calib_df, crs = st_crs(sec_veg_age_rst))






raster_sf <- st_as_sf(
  as.data.frame(cbind(
    terra::xyFromCell(resampled_agb_raster, seq_len(ncell(resampled_agb_raster))), # Extract coordinates
    agb_2017 = terra::values(resampled_agb_raster) # Extract AGB values
  ), na.rm = TRUE), # Remove NA values
  coords = c("x", "y"), # Specify coordinate columns
  crs = st_crs(sec_veg_age_rst), # Set CRS
  remove = FALSE # Retain coordinate columns
)











