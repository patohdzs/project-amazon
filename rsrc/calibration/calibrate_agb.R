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

# Load raster data
clean_mapbiomas <- rast("data/clean/land_use_cover_2000.tif")

# Load secondary vegetation age data for 2017
sec_veg_age_rst <- rast(
  list.files(
    "data/raw/mapbiomas/secondary_vegetation_age/",
    pattern = "2017",
    full.names = TRUE
  )
)

# Load aboveground biomass raster data
agb_rasters <- map(
  list.files(
    "data/raw/esa/above_ground_biomass/",
    pattern = "_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv3.0.tif",
    full.names = TRUE
  ),
  rast
)

# Merge AGB reginonal rasters into a single raster
agb_rst <- do.call(merge, agb_rasters)

# Crop the the Amazonia subset
sec_veg_age_rst <- sec_veg_age_rst |>
  crop(clean_mapbiomas) |>
  mask(clean_mapbiomas)

# Resample to match sec_veg resolution and mask values outside Amazonia
agb_rst <- resample(agb_rst, sec_veg_age_rst, method = "near") |>
  mask(clean_mapbiomas)

rsts <- c(agb_rst, sec_veg_age_rst)
writeRaster(agb_rst, "data/clean/clean_agb_2017.tif", overwrite = TRUE)
writeRaster(sec_veg_age_rst, "data/clean/clean_sec_veg_age_2017.tif", overwrite = TRUE)



# Remove NA values
amazon_mask <- !is.na(agb_rst) & !is.na(sec_veg_age_rst)

# Get values
agb <- values(agb_rst)[amazon_mask]
sec_veg_age <- values(sec_veg_age_rst)[amazon_mask]



# Get pixels into sf
df <- c(agb_rst, sec_veg_age_rst) %>%
  as.polygons() %>%
  st_as_sf()


# agb_rst <- agb_rst |>
#  crop(clean_mapbiomas) |>
#  mask(clean_mapbiomas)
#

# writeRaster(resampled_agb_raster, "data/calibration/resampled_agb_2017.tif", overwrite = TRUE)
#
#
# load("data/calibration/calibration_1043_sites.Rdata")
# calib_df <- st_transform(calib_df, crs = st_crs(sec_veg_age_rst))
#
# raster_sf <- st_as_sf(
#  as.data.frame(cbind(
#    terra::xyFromCell(resampled_agb_raster, seq_len(ncell(resampled_agb_raster))), # Extract coordinates
#    agb_2017 = terra::values(resampled_agb_raster) # Extract AGB values
#  ), na.rm = TRUE), # Remove NA values
#  coords = c("x", "y"), # Specify coordinate columns
#  crs = st_crs(sec_veg_age_rst), # Set CRS
#  remove = FALSE # Retain coordinate columns
# )
#
