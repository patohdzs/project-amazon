
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


setwd("C:/Users/pengyu/Desktop/code_data_20230628")



# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("code/setup.R")


# START TIMER
tictoc::tic(msg = "pixelBiomass2017_prepData.R script", log = T)


# RASTER OPTIONS
terra::terraOptions(tmpdir = here::here("data/_temp"),
                     timer  = T)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------



# MAPBIOMAS SAMPLE
load("data/calibration/prepData/pixelPrimaryForest2018_prepData.Rdata")



# ABOVEGROUND BIOMASS RASTERS
path <- "data/raw2clean/abovegroundBiomass_esa/input"


agb.raster <- purrr::map(list.files(path, pattern = "(_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv3.0.tif)", full.names = T),
                         terra::rast)

agb_sd_raster <- purrr::map(list.files(path, pattern = "(_ESACCI-BIOMASS-L4-AGB_SD-MERGED-100m-2017-fv3.0.tif)", full.names = T),
                            terra::rast)


# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# transform to sf
pixelPrimaryForest.prepData <- sf::st_as_sf(x = pixelPrimaryForest.prepData,
                                            coords = c("lon", "lat"),
                                            remove = FALSE,
                                            crs = sf::st_crs(4326))


extract_and_merge <- function(raster_list, column_name) {
  data <- purrr::map_df(.x = seq_along(raster_list),
                        .f = function(.x) {
                          
                          aux.polygons <- terra::vect(pixelPrimaryForest.prepData)
                          aux.polygons <- terra::crop(aux.polygons, raster_list[[.x]])
                          
                          if (is.null(aux.polygons)) {
                            return(NULL)  # Changed from 'next()' for safe data.frame return
                          }
                          
                          names(raster_list[[.x]]) <- column_name
                          aux.data <- terra::extract(raster_list[[.x]], aux.polygons, xy = TRUE)
                          aux.data <- aux.data %>% dplyr::rename(lon = x, lat = y) %>% dplyr::select(-ID)
                          
                          return(aux.data)
                        })
  
  return(data)
}


pixelBiomass2017.prepData <- extract_and_merge(agb.raster, "agb_2017")
pixelBiomass2017_sd.prepData <- extract_and_merge(agb_sd_raster, "agbSD_2017")
pixelBiomass2017.prepData <- dplyr::left_join(pixelBiomass2017.prepData, pixelBiomass2017_sd.prepData, by = c("lon", "lat"))

pixelBiomass2017.prepData <-
  sf::st_as_sf(x = pixelBiomass2017.prepData,
               coords = c("lon", "lat"),
               crs = sf::st_crs(4326))


# check if all points were extracted
if (length(pixelBiomass2017.prepData$lon) != length(pixelPrimaryForest.prepData$lon)) {
  print("Different number of spatial points between original and final sample! Check script!")
}

# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(pixelBiomass2017.prepData$agb_2017) <- "aboveground biomass in 2017 (Mg per ha)"



# POST-TREATMENT OVERVIEW
# summary(pixelBiomass2017.prepData)
# View(pixelBiomass2017.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

file_path <- paste(getwd(), "data/calibration/prepData", paste0("pixelBiomass2017_prepData", ".Rdata"), sep = "/")

save(pixelBiomass2017.prepData, file = file_path)


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("code/calibration")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
