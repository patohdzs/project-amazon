
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




# START TIMER
tictoc::tic(msg = "pixelArea_1043SitesModel.R script", log = TRUE)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tmpdir = "data/_temp",
                      timer  = T)





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA
raster_biome <- terra::rast("data/calibration/1043SitesModel/aux_tifs/raster_amazonBiome_1043SitesModel.tif")





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# rasterize amazon biome
raster_pixelArea <- terra::cellSize(raster_biome, unit = "ha")

# clean environment
rm(raster_biome)

# add name
names(raster_pixelArea) <- "pixelArea_ha"



# EXPORT
# save unified tif
terra::writeRaster(raster_pixelArea, "data/calibration/1043SitesModel/aux_tifs/raster_pixelArea_1043SitesModel.tif", overwrite = T)

# clean environment
rm(raster_pixelArea)



# CLEAN TEMP DIR
terra::tmpFiles(current = T, remove = T)
gc()



# END TIMER
tictoc::toc(log = TRUE)



# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------