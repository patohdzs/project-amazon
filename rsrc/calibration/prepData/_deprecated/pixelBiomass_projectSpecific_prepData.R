
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: ADD ABOVEGROUND BIOMASS DATA (BACCINI) TO MAPBIOMAS 30M-PIXELS
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
here::i_am("code/projectSpecific/prepData/pixelBiomass_projectSpecific_prepData.R", uuid = "19c04f37-96dc-4e8d-a2ae-ade294205749")


# START TIME
tictoc::tic(msg = "pixelBiomass_projectSpecific_prepData script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary
groundhog::groundhog.library("sp", groundhog.date) # to manipulate spatial data (vector format) interacting with raster data
groundhog::groundhog.library("raster", groundhog.date)  # manipulate spatial data (raster format)
groundhog::groundhog.library("rgdal", groundhog.date)  # to read raster
groundhog::groundhog.library("rgeos", groundhog.date)  # to support spatial manipulation
groundhog::groundhog.library("sf", groundhog.date)  # manipulate spatial data (vector format)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# MAPBIOMAS SAMPLE
load(here::here("data/projectSpecific/prepData/pixelCategories_prepData_2000.Rdata"))



# ABOVEGROUND BIOMASS RASTERS
agb.raster <- lapply(list.files(here::here("data/raw2clean/abovegroundBiomass_gfw/input"), pattern = "aboveground_biomass_ha_2000.tif", full.names = T),
                     raster::raster)





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# transform matrix to spatialPoints
pixelCategories.prepData.2000 <- sp::SpatialPointsDataFrame(coords = pixelCategories.prepData.2000[, 1:2],
                                                                         data = data.frame(pixelCategories.prepData.2000),
                                                      proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))

# create empty sf to be filled in the loop
pixelBiomass.prepData <- sf::st_as_sf(SpatialPointsDataFrame(data.frame(x = 0, y = 0), data.frame(x = 0, y = 0, agb = 0),
                              proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))[-1,])



# MERGE AGB DATA WITH MAPBIOMAS SAMPLE
for (i in seq_along(agb.raster)) {

  # crop spatial points to raster extent
  aux.sample <- raster::crop(pixelCategories.prepData.2000, agb.raster[[i]])

  # skip raster with no intersection with the sample
  if (is.null(aux.sample)) {
    next()
  }

  # change raster layer name
  names(agb.raster[[i]]) <- "agb"

  # extract agb raster data by spatial point and transform to sf
  aux.agb <- sf::st_as_sf(raster::extract(agb.raster[[i]], aux.sample, sp = T))

  # merge spatial points (sf)
  pixelBiomass.prepData <- rbind(pixelBiomass.prepData, aux.agb)

}

# clear environment
rm(agb.raster, aux.agb, aux.sample)

# check if all points were extracted
if (length(pixelBiomass.prepData$lon) != length(pixelCategories.prepData.2000$lon)) {
  print("Different number of spatial points between original and final sample! Check script!")
}





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(pixelBiomass.prepData$agb) <- "aboveground biomass in 2000 (Mg per ha)"



# POST-TREATMENT OVERVIEW
# summary(pixelBiomass.prepData)
# View(pixelBiomass.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(pixelBiomass.prepData,
     file = here::here("data/projectSpecific/prepData",
                      paste0("pixelBiomass_prepData", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/prepData")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
