
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: ADD ABOVEGROUND BIOMASS DATA (KAYRROS) TO MAPBIOMAS 30M-PIXELS
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# GROUNDHOG (REPRODUCIBILITY SOLUTION TO HANDLING DIFFERENT VERSIONS OF R AND ITS PACKAGES)

# check if groundhog is installed and load it
if ("groundhog" %in% installed.packages()) {
  library("groundhog")
} else {
  install.packages("groundhog")
  library("groundhog")
}

# define date of reference to load all packages
groundhog.date <- "2022-04-01"

# guarantee version 1.5 of groundhog is being used
groundhog::meta.groundhog(date = "2022-04-01")


# HERE
groundhog::groundhog.library("here", groundhog.date) # load package here


# TICTOC
groundhog::groundhog.library("tictoc", groundhog.date) # load package tictoc


# DECLARE LOCATION OF CURRENT SCRIPT TO SET UP PROJECT ROOT CORRECTLY
here::i_am("code/projectSpecific/prepData/pixelBiomass2021_projectSpecific_prepData.R", uuid = "19c04f37-96dc-4e8d-a2ae-ade294205749")


# START TIME
tictoc::tic(msg = "pixelBiomass2021_projectSpecific_prepData script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary
groundhog::groundhog.library("terra", groundhog.date)  # manipulate spatial data (raster format)
groundhog::groundhog.library("sf", groundhog.date)  # manipulate spatial data (vector format)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# MAPBIOMAS PRIMARY FOREST SAMPLE
load(here::here("data/projectSpecific/prepData/pixelPrimaryForest2019_prepData.Rdata"))



# ABOVEGROUND BIOMASS RASTERS
agb.raster <- purrr::map(list.files(here::here("data/raw2clean/abovegroundBiomass_kayrros/input/2021"), pattern = ".tif", full.names = T),
                         terra::rast)


# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# transform to sf
pixelPrimaryForest.prepData <- sf::st_as_sf(x = pixelPrimaryForest.prepData,
                                            coords = c("lon", "lat"),
                                            remove = TRUE,
                                            crs = sf::st_crs(4326))

# add id
pixelPrimaryForest.prepData$id <- 1:nrow(pixelPrimaryForest.prepData)

# remove unnecessary columns
pixelPrimaryForest.prepData <- pixelPrimaryForest.prepData %>% dplyr::select(id, geometry)


# MERGE AGB DATA WITH MAPBIOMAS SAMPLE
pixelBiomass2021.prepData <-
  purrr::map_df(.x = seq_along(agb.raster),
                .f = function(.x) {

                  # adjust projection to match with raster
                  aux.polygons <- sf::st_transform(x = pixelPrimaryForest.prepData, crs = sf::st_crs(agb.raster[[.x]]))

                  # transform to spatVector
                  aux.polygons <- terra::vect(aux.polygons)

                  # crop spatial points to raster extent
                  aux.polygons <- terra::crop(aux.polygons, agb.raster[[.x]])

                  # skip raster with no intersection with the sample
                  if (nrow(aux.polygons) == 0) {
                    return(NULL)
                  }

                  # change raster layer name
                  names(agb.raster[[.x]]) <- "agb_2021"

                  # extract agb raster data by spatial point
                  aux.agb <- dplyr::bind_cols(id = aux.polygons$id, terra::extract(agb.raster[[.x]], aux.polygons, xy = TRUE))

                  # adjust and select column names
                  aux.agb <- aux.agb %>% dplyr::rename(lon = x, lat = y) %>% dplyr::select(-ID)


                  # transform to sf
                  aux.agb <- sf::st_as_sf(x = aux.agb,
                                          coords = c("lon", "lat"),
                                          crs = sf::st_crs(agb.raster[[.x]]))

                  # adjust projection
                  aux.agb <- sf::st_transform(aux.agb, sf::st_crs(4326))

                  return(aux.agb)
                })

# clear environment
rm(agb.raster)

# remove duplicates
pixelBiomass2021.prepData <- dplyr::distinct(pixelBiomass2021.prepData, id, .keep_all = TRUE)

# check if all points were extracted
if (length(pixelBiomass2021.prepData$id) != length(pixelPrimaryForest.prepData$id)) {
  print("Different number of spatial points between original and final sample! Check script!")
}





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(pixelBiomass2021.prepData$agb_2021) <- "aboveground biomass in 2021 (Mg per ha)"



# POST-TREATMENT OVERVIEW
# summary(pixelBiomass2021.prepData)
# View(pixelBiomass2021.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(pixelBiomass2021.prepData,
     file = here::here("data/projectSpecific/prepData",
                       paste0("pixelBiomass2021_prepData", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/prepData")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
