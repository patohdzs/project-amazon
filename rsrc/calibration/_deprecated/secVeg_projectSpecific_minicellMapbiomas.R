
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: ADD SECONDARY VEGETATION DATA (FLORESER - 2000)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "secVeg_projectSpecific_minicellMapbiomas script", log = T)



# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(sp) # manipulate spatial data
library(sf) # manipulate spatial data
library(raster) # manipulate raster data
library(rgdal) # manipulate spatial data
library(rgeos) # manipulate spatial data





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# MAPBIOMAS SAMPLE
load("data/projectSpecific/minicellMapbiomas/mapbiomasCategories_minicellMapbiomas_2000.Rdata")



# SEC VEG RASTER
clean.floreser <- raster::raster("data/raw2clean/floreser_imazon/output/2000/clean_floreser.tif")






# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# transform matrix to spatialPoints
mapbiomasCategories.minicellMapbiomas.2000 <- sp::SpatialPointsDataFrame(coords = mapbiomasCategories.minicellMapbiomas.2000[, 1:2], data = data.frame(mapbiomasCategories.minicellMapbiomas.2000),
                                                      proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))



# MERGE SECVEG DATA WITH MAPBIOMAS SAMPLE

# change raster layer name
names(clean.floreser) <- "secVeg_age"

# extract secVeg raster data by spatial point and transform to sf
secVeg.projectSpecific.minicellMapbiomas.2000 <- sf::st_as_sf(raster::extract(clean.floreser, mapbiomasCategories.minicellMapbiomas.2000, sp = T))

# check if all points were extracted
if (length(secVeg.projectSpecific.minicellMapbiomas.2000$lon) != length(mapbiomasCategories.minicellMapbiomas.2000$lon)) {
  print("Different number of spatial points between original and final sample! Check script!")
}

# clear environment
rm(clean.floreser, mapbiomasCategories.minicellMapbiomas.2000)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(secVeg.projectSpecific.minicellMapbiomas.2000$secVeg_age)         <- "secondary forest age (years)"



# POST-TREATMENT OVERVIEW
# summary(secVeg.projectSpecific.minicellMapbiomas)
# View(secVeg.projectSpecific.minicellMapbiomas)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(secVeg.projectSpecific.minicellMapbiomas.2000,
     file = file.path("data/projectSpecific/minicellMapbiomas",
                      paste0("secVeg_minicellMapbiomas_2000", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/minicellMapbiomas")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
