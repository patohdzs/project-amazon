
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CALCULATE SITE AREA - 1055 SITES
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "pixelArea_1043SitesModel.R script", log = T)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA
raster.biome <- terra::rast(here::here("data/calibration/1043SitesModel/aux_tifs/raster_amazonBiome_1043SitesModel.tif"))





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# rasterize amazon biome
raster.pixelArea <- terra::cellSize(raster.biome, unit = "ha")

# clean environment
rm(raster.biome)

# add name
names(raster.pixelArea) <- "pixelArea_ha"



# EXPORT
# save unified tif
terra::writeRaster(raster.pixelArea, here::here("data/calibration/1043SitesModel/aux_tifs/raster_pixelArea_1043SitesModel.tif"), overwrite = T)

# clean environment
rm(raster.pixelArea)



# CLEAN TEMP DIR
terra::tmpFiles(current = T, remove = T)
gc()



# END TIMER
tictoc::toc(log = T)

# # export time to csv table
# ExportTimeProcessing("code/calibration")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------