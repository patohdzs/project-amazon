
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




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("code/setup.R")


# START TIMER
tictoc::tic(msg = "amazonBiome_1055SitesModel.R script", log = T)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA - TO USE AS MASK (30M MINICELLS)
raster.biome <- terra::rast(here::here("data/raw2clean/landUseCover_mapbiomas/input/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-2017.tif"))



# AMAZON BIOME VECTOR DATA
load(here::here("data/raw2clean/amazonBiome_ibge/output/clean_amazonBiome.Rdata"))





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# change projection to match raster
clean.amazonBiome <- sf::st_transform(clean.amazonBiome, crs(raster.biome))


# rasterize amazon biome into 30m raster resolution to minimize area distortion
raster.biome <- terra::rasterize(terra::vect(clean.amazonBiome), raster.biome, field = 1)

# clean environment
rm(clean.amazonBiome)

# aggregate raster calculating the share of minicells that are in the biome
raster.biome <- terra::aggregate(raster.biome, fact = 2250, fun = sum, na.rm = T)/(2250^2) # (2250^2) is the total number of minicells

# add name
names(raster.biome) <- "share_amazonBiome"


# EXPORT
# save unified tif
terra::writeRaster(raster.biome, here::here("data/calibration/1055SitesModel/aux_tifs/raster_amazonBiome_1055SitesModel.tif"), overwrite = T)

# clean environment
rm(raster.biome)



# CLEAN TEMP DIR
terra::tmpFiles(current = T, remove = T)
gc()



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("code/calibration")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------