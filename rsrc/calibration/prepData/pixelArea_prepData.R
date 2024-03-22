
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CALCULATE AREA OF MAPBIOMAS 30M-PIXELS
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "pixelArea_prepData.R script", log = T)


# RASTER OPTIONS
terra::terraOptions(tmpdir = here::here("data/_temp"),
                      timer  = T)






# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA
clean.mapbiomas <- terra::rast(here::here("data/raw2clean/landUseCover_mapbiomas/output/clean_landUseCover_2000.tif"))



# MAPBIOMAS 30M-PIXELS SAMPLE
load(here::here("data/calibration/prepData/samplePixel_prepData.Rdata"))





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# calculate area of each pixel (sq km)
clean.mapbiomasArea <- terra::cellSize(clean.mapbiomas, unit = "ha")

# change raster layer name
names(clean.mapbiomasArea) <- "pixel_area"

# clean environment
rm(clean.mapbiomas)

# extract one year from sample to extract pixel area
samplePixel.prepData.2000 <-
  samplePixel.prepData %>%
  dplyr::filter(year == 2000) %>%
  dplyr::select(-year, -mapbiomas_class) # remove unncessary columns

# transform to spatialPoints
samplePixel.prepData.2000.sf <- sf::st_as_sf(x = samplePixel.prepData.2000,
                                            coords = c("lon", "lat"),
                                            remove = FALSE,
                                            crs = sf::st_crs(4326))

# extract pixel area raster data for sample points
aux.pixelArea <- terra::extract(clean.mapbiomasArea, samplePixel.prepData.2000.sf)

# clean environment
rm(samplePixel.prepData.2000.sf)

# merge pixel area variable with sample 2000
samplePixel.prepData.2000$pixel_area <- aux.pixelArea

# merge sample 2000 with panel sample
pixelArea.prepData <-
  samplePixel.prepData %>%
  dplyr::left_join(samplePixel.prepData.2000)

# clear environmnet
rm(samplePixel.prepData, samplePixel.prepData.2000, aux.pixelArea)



# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(pixelArea.prepData$pixel_area) <- "pixel area (ha)"



# POST-TREATMENT OVERVIEW
# summary(pixelArea.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(pixelArea.prepData,
     file = here::here("data/calibration/prepData",
                      paste0("pixelArea_prepData", ".Rdata")))



# # END TIMER
# tictoc::toc(log = T)

# # export time to csv table
# ExportTimeProcessing("code/calibration")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------