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


# START TIMER
tictoc::tic(msg = "landUseCover_raw2clean.R script", log = TRUE)
tictoc::tic(msg = "landUseCover_raw2clean.R script", log = TRUE)

# RASTER OPTIONS
terra::terraOptions(
  tmpdir = "data/_temp",
  timer = T
)

# DATA INPUT
# RAW DATA
# read only year 2000 - used as the base for random sample extraction
raw_raster <- terra::rast("data/raw/mapbiomas/land_use_cover/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-2000.tif")

# DATASET CLEANUP AND PREP
# RECLASSIFY
# change 0s to NAs so that any value represents only areas inside the Amazon Biome
raw_raster <- terra::subst(raw_raster, from = 0, to = as.numeric(NA))

# EXPORT
# save reclassified tif
terra::writeRaster(raw_raster, "data/clean/landusecover_2000.tif", overwrite = T)

# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()

# END TIMER
tictoc::toc(log = TRUE)


