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

# RASTER OPTIONS
terra::terraOptions(
  tmpdir = here::here("data/_temp"),
  timer = TRUE
)

# Read only year 2000 - used as the base for random sample extraction
raw_raster <- terra::rast("data/raw/mapbiomas/land_use_cover/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-2000.tif")

# Change 0s to NAs so that any value represents only areas inside the Amazon Biome
raw_raster <- terra::subst(raw.raster, from = 0, to = as.numeric(NA))

# Save reclassified tif
terra::writeRaster(raw.raster, here::here("data/clean/land_use_cover_2000.tif"), overwrite = TRUE)

# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()

# END TIMER
tictoc::toc(log = T)
