
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: ADD 2018 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("code/setup.R")


# START TIMER
tictoc::tic(msg = "pixelBiomass2018_prepData.R script", log = T)


# RASTER OPTIONS
terra::terraOptions(tmpdir = here::here("data/_temp"),
                     timer  = T)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# MAPBIOMAS SAMPLE
load(here::here("data/calibration/prepData/pixelPrimaryForest2018_prepData.Rdata"))



# ABOVEGROUND BIOMASS RASTERS
agb.raster <- purrr::map(list.files(here::here("data/raw2clean/abovegroundBiomass_esa/input"), pattern = "_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2018-fv3.0.tif", full.names = T),
                         terra::rast)




# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# transform to sf
pixelPrimaryForest.prepData <- sf::st_as_sf(x = pixelPrimaryForest.prepData,
                                            coords = c("lon", "lat"),
                                            remove = FALSE,
                                            crs = sf::st_crs(4326))



# MERGE AGB DATA WITH MAPBIOMAS SAMPLE
pixelBiomass2018.prepData <-
  purrr::map_df(.x = seq_along(agb.raster),
                .f = function(.x) {

                # transform to spatVector
                aux.polygons <- terra::vect(pixelPrimaryForest.prepData)

                # crop spatial points to raster extent
                aux.polygons <- terra::crop(aux.polygons, agb.raster[[.x]])

                # skip raster with no intersection with the sample
                if (is.null(aux.polygons)) {
                  next()
                }

                # change raster layer name
                names(agb.raster[[.x]]) <- "agb_2018"

                # extract agb raster data by spatial point
                aux.agb <- terra::extract(agb.raster[[.x]], aux.polygons, xy = TRUE)

                # adjust and select column names
                aux.agb <- aux.agb %>% dplyr::rename(lon = x, lat = y) %>% dplyr::select(-ID)

              })


# clear environment
rm(agb.raster)

# transform to sf
pixelBiomass2018.prepData <-
  sf::st_as_sf(x = pixelBiomass2018.prepData,
               coords = c("lon", "lat"),
               crs = sf::st_crs(4326))



# check if all points were extracted
if (length(pixelBiomass2018.prepData$lon) != length(pixelPrimaryForest.prepData$lon)) {
  print("Different number of spatial points between original and final sample! Check script!")
}


# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(pixelBiomass2018.prepData$agb_2018) <- "aboveground biomass in 2018 (Mg per ha)"



# POST-TREATMENT OVERVIEW
# summary(pixelBiomass2018.prepData)
# View(pixelBiomass2018.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(pixelBiomass2018.prepData,
     file = here::here("data/calibration/prepData",
                       paste0("pixelBiomass2018_prepData", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("code/calibration")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
