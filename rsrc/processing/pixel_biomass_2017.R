
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: ADD 2017 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -



# START TIMER
tictoc::tic(msg = "pixelBiomass2017_prepData.R script", log = TRUE)


# RASTER OPTIONS
terra::terraOptions(tmpdir = "data/_temp",
                     timer  = T)


# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# MAPBIOMAS SAMPLE
load("data/prepData/pixelPrimaryForest2018_prepData.Rdata")



# ABOVEGROUND BIOMASS RASTERS
agb_raster <- purrr::map(list.files("data/raw/esa/above_ground_biomass/", pattern = "_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv3.0.tif", full.names = T),
                         terra::rast)




# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# transform to sf
pixelPrimaryForest_prepData <- sf::st_as_sf(x = pixelPrimaryForest_prepData,
                                            coords = c("lon", "lat"),
                                            remove = FALSE,
                                            crs = sf::st_crs(4326))



# MERGE AGB DATA WITH MAPBIOMAS SAMPLE
pixelBiomass2017_prepData <-
  purrr::map_df(.x = seq_along(agb_raster),
                .f = function(.x) {

                # transform to spatVector
                aux_polygons <- terra::vect(pixelPrimaryForest_prepData)

                # crop spatial points to raster extent
                aux_polygons <- terra::crop(aux_polygons, agb_raster[[.x]])

                # skip raster with no intersection with the sample
                if (is.null(aux_polygons)) {
                  next()
                }

                # change raster layer name
                names(agb_raster[[.x]]) <- "agb_2017"

                # extract agb raster data by spatial point
                aux_agb <- terra::extract(agb_raster[[.x]], aux_polygons, xy = TRUE)

                # adjust and select column names
                aux_agb <- aux_agb %>% dplyr::rename(lon = x, lat = y) %>% dplyr::select(-ID)

              })


# clear environment
rm(agb_raster)

# transform to sf
pixelBiomass2017_prepData <-
  sf::st_as_sf(x = pixelBiomass2017_prepData,
               coords = c("lon", "lat"),
               crs = sf::st_crs(4326))



# check if all points were extracted
if (length(pixelBiomass2017_prepData$lon) != length(pixelPrimaryForest_prepData$lon)) {
  print("Different number of spatial points between original and final sample! Check script!")
}

# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(pixelBiomass2017_prepData$agb_2017) <- "aboveground biomass in 2017 (Mg per ha)"




# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(pixelBiomass2017_prepData,
     file = "data/prepData/pixelBiomass2017_prepData.Rdata")



# END TIMER
tictoc::toc(log = TRUE)


# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
