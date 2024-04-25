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

# SETUP
# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")

# START TIMER
tictoc::tic(msg = "precipitation_raw2clean.R script", log = T)

# RASTER OPTIONS
terra::terraOptions(
  tmpdir = here::here("data/_temp"),
  timer = T
)

# DATA INPUT
# RAW DATA
# read all parts (12) of the year and merge them together
raw.raster <- terra::rast(list.files(here::here("data/raw2clean/precipitation_worldClim/input"),
  full.names = T,
  pattern = ".tif"
))

# DATASET CLEANUP AND PREP
# change layer names
names(raw.raster) <- paste0("historicalPrecip_", 1:12)

# EXPORT
# save unified tif
terra::writeRaster(raw.raster, here::here("data/raw2clean/precipitation_worldClim/output/clean_precipitation.tif"), overwrite = T)

# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()

# END TIMER
tictoc::toc(log = T)

# export time to csv table
# ExportTimeProcessing("code/raw2clean")
