
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: GENERATE AGGREGATED SAMPLE OF INTEREST (DIVIDE AMAZON BIOME INTO 1055 CELLS)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -




# START TIMER
tictoc::tic(msg = "amazonBiome_1043SitesModel.R script", log = TRUE)

# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tmpdir = "data/_temp",
                      timer  = T)

# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA - TO USE AS MASK (30M MINICELLS)
raster_biome <- terra::rast("data/raw/mapbiomas/land_use_cover/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-2017.tif")



# AMAZON BIOME VECTOR DATA
load("data/clean/amazon_biome.Rdata")




# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# change projection to match raster
amazon_biome <- sf::st_transform(amazon_biome, crs(raster_biome))


# rasterize amazon biome into 30m raster resolution to minimize area distortion
raster_biome <- terra::rasterize(terra::vect(amazon_biome), raster_biome, field = 1)

# clean environment
rm(amazon_biome)

# aggregate raster calculating the share of minicells that are in the biome
raster_biome <- terra::aggregate(raster_biome, fact = 2250, fun = sum, na.rm = T)/(2250^2) # (2250^2) is the total number of minicells

# add name
names(raster_biome) <- "share_amazonBiome"


if (!dir.exists("data/calibration/1043SitesModel/aux_tifs")) {
    dir.create("data/calibration/1043SitesModel/aux_tifs", recursive = TRUE)
}
# EXPORT
# save unified tif
terra::writeRaster(raster_biome, "data/calibration/1043SitesModel/aux_tifs/raster_amazonBiome_1043SitesModel.tif", overwrite = T)

# clean environment
rm(raster_biome)



# CLEAN TEMP DIR
terra::tmpFiles(current = T, remove = T)
gc()



# END TIMER
tictoc::toc(log = TRUE)






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
