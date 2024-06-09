
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



# START TIMER
tictoc::tic(msg = "pixelArea_prepData.R script", log = TRUE)


# RASTER OPTIONS
terra::terraOptions(tmpdir = "data/_temp",
                      timer  = T)


# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA
clean_mapbiomas <- terra::rast("data/clean/landusecover_2000.tif")


# MAPBIOMAS 30M-PIXELS SAMPLE
load("data/prepData/samplePixel_prepData.Rdata")



# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# calculate area of each pixel (sq km)
clean_mapbiomasArea <- terra::cellSize(clean_mapbiomas, unit = "ha")

# change raster layer name
names(clean_mapbiomasArea) <- "pixel_area"

# clean environment
rm(clean_mapbiomas)

# extract one year from sample to extract pixel area
samplePixel_prepData_2000 <-
  samplePixel_prepData %>%
  dplyr::filter(year == 2000) %>%
  dplyr::select(-year, -mapbiomas_class) # remove unncessary columns

# transform to spatialPoints
samplePixel_prepData_2000.sf <- sf::st_as_sf(x = samplePixel_prepData_2000,
                                            coords = c("lon", "lat"),
                                            remove = FALSE,
                                            crs = sf::st_crs(4326))

# extract pixel area raster data for sample points
aux_pixelArea <- terra::extract(clean_mapbiomasArea, samplePixel_prepData_2000.sf)

# clean environment
rm(samplePixel_prepData_2000.sf)

# merge pixel area variable with sample 2000
samplePixel_prepData_2000$pixel_area <- aux_pixelArea

# merge sample 2000 with panel sample
pixelArea_prepData <-
  samplePixel_prepData %>%
  dplyr::left_join(samplePixel_prepData_2000)

# clear environmnet
rm(samplePixel_prepData, samplePixel_prepData_2000, aux_pixelArea)



# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()



# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(pixelArea_prepData$pixel_area) <- "pixel area (ha)"



# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(pixelArea_prepData,
     file = "data/prepData/pixelArea_prepData.Rdata")

# END TIMER
tictoc::toc(log = TRUE)



# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------