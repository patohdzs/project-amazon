
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: ADD ABOVEGROUND BIOMASS DATA (BACCINI)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "agb_projectSpecific_minicellFloreser script", log = T)



# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(sp)     # manipulate spatial data
library(sf)     # manipulate spatial data
library(raster) # manipulate raster data
library(rgdal)  # manipulate spatial data
library(rgeos)  # manipulate spatial data
library(tidyverse)  # manipulate tables, works with sf





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# FLORESER SAMPLE
load("data/projectSpecific/minicellFloreser/pixelArea_minicellFloreser.Rdata")



# ABOVEGROUND BIOMASS RASTERS
agb.raster <- lapply(list.files("data/raw2clean/abovegroundBiomass_gfw/input", pattern = "aboveground_biomass_ha_2000.tif", full.names = T),
                     raster::raster)





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# transform matrix to spatialPoints
pixelArea.minicellFloreser <- sp::SpatialPointsDataFrame(coords = data.frame(pixelArea.minicellFloreser[, 1:2]), data = data.frame(pixelArea.minicellFloreser),
                                                      proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))

# create empty sf to be filled in the loop
agb <- sf::st_as_sf(SpatialPointsDataFrame(data.frame(x = 0, y = 0), data.frame(x = 0, y = 0, clean_floreser = 0, agb = 0, pixel_area = 0),
                              proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))[-1,])



# MERGE AGB DATA WITH FLORESER SAMPLE
for (i in seq_along(agb.raster)) {

  # crop spatial points to raster extent
  aux.sample <- raster::crop(pixelArea.minicellFloreser, agb.raster[[i]])

  # skip raster with no intersection with the sample
  if (is.null(aux.sample)) {
    next()
  }

  # change raster layer name
  names(agb.raster[[i]]) <- "agb"

  # extract agb raster data by spatial point and transform to sf
  aux.agb <- sf::st_as_sf(raster::extract(agb.raster[[i]], aux.sample, sp = T))

  # merge spatial points (sf)
  agb <- rbind(agb, aux.agb)

}

# clear environment
rm(agb.raster, aux.agb, aux.sample)

# check if all points were extracted
if (length(agb) != length(pixelArea.minicellFloreser)) {
  print("Different number of spatial points between original and final sample! Check script!")
}

# adjust column names
agb.projectSpecific.minicellFloreser <-
  agb %>%
  dplyr::rename(lon = x, lat = y, age_secveg = clean_floreser)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(agb.projectSpecific.minicellFloreser$lon)        <- "longitude of the cell centroid (degrees)"
sjlabelled::set_label(agb.projectSpecific.minicellFloreser$lat)        <- "latitude of the cell centroid (degrees)"
sjlabelled::set_label(agb.projectSpecific.minicellFloreser$age_secveg) <- "secondary vegetation age in 2000 (year)"
sjlabelled::set_label(agb.projectSpecific.minicellFloreser$agb)        <- "aboveground biomass in 2000 (Mg per ha)"



# POST-TREATMENT OVERVIEW
# summary(agb.projectSpecific.minicellFloreser)
# View(agb.projectSpecific.minicellFloreser)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(agb.projectSpecific.minicellFloreser,
     file = file.path("data/projectSpecific/minicellFloreser",
                      paste0("agb_minicellFloreser", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/minicellFloreser")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
