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
  timer = TRUE
)

# Read only year 2000 - used as the base for random sample extraction
in_path <- "data/raw/mapbiomas/land_use_cover/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-2000.tif"
raw_raster <- terra::rast(in_path)

# Change 0s to NAs to represent only areas inside the Amazon Biome
raw_raster <- terra::subst(raw_raster, from = 0, to = as.numeric(NA))

# Save reclassified tif
out_path <- "data/clean/land_use_cover_2000.tif"
terra::writeRaster(raw_raster, out_path, overwrite = TRUE)

# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()

# END TIMER
tictoc::toc(log = TRUE)
