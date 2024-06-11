# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: UNIFY MONTHLY TIF FILES - HISTORICAL PRECIPITATION (WORLD CLIM)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(terra)

# START TIMER
tic(msg = "precipitation_raw2clean.R script", log = TRUE)

# Read and merge all layers (12) of the year
raw_raster <- rast(list.files(
  "data/raw/worldclim/precipitation/",
  full.names = TRUE,
  pattern = ".tif"
))

# Change layer names
names(raw_raster) <- paste0("historicalPrecip_", 1:12)

# Save unified raster
out_path <- "data/clean/precipitation.tif"
writeRaster(raw_raster, out_path, overwrite = TRUE)

# CLEAN TEMP DIR
tmpFiles(current = TRUE, remove = TRUE)
gc()

# END TIMER
toc(log = TRUE)
