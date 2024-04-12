# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: RECLASSIFY PIXELS OUTSIDE BIOME BOUNDARY - LAND USE AND COVER (MAPBIOMAS - COL5)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

# SETUP

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")

# START TIMER
tictoc::tic(msg = "landUseCover_raw2clean.R script", log = T)

# RASTER OPTIONS
terra::terraOptions(
  tmpdir = here::here("data/_temp"),
  timer = T
)

# DATA INPUT
# RAW DATA
# read only year 2000 - used as the base for random sample extraction
raw.raster <- terra::rast(here::here("data/raw2clean/landUseCover_mapbiomas/input/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-2000.tif"))

# DATASET CLEANUP AND PREP
# RECLASSIFY
# change 0s to NAs so that any value represents only areas inside the Amazon Biome
raw.raster <- terra::subst(raw.raster, from = 0, to = as.numeric(NA))

# EXPORT
# save reclassified tif
terra::writeRaster(raw.raster, here::here("data/raw2clean/landUseCover_mapbiomas/output/clean_landUseCover_2000.tif"), overwrite = T)

# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()

# END TIMER
tictoc::toc(log = T)

# export time to csv table
# ExportTimeProcessing("code/raw2clean")
