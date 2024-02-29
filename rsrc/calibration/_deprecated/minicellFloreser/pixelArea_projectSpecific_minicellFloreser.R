
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CALCULATE AREA OF MINICELLS
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "pixelArea_projectSpecific_minicellFloreser script", log = T)



# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(raster) # for raster manipulation
library(tidyverse) # for data frame manipulation



# RASTER OPTIONS
raster::rasterOptions(tmpdir = file.path("data/_temp"),
                      timer  = T)





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA
clean.floreser <- raster::raster("data/raw2clean/floreser_imazon/output/2000/clean_floreser.tif")



# FLORESER SAMPLE
load("data/projectSpecific/minicellFloreser/sample_minicellFloreser.Rdata")





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# calculate area of each pixel (sq km)
clean.floreserArea <- raster::area(clean.floreser)

# convert pixel area from sq km to ha
clean.floreserArea <- clean.floreserArea*100

# change raster layer name
names(clean.floreserArea) <- "pixel_area"

# clean environment
rm(clean.floreser)

# transform to spatialPoints
sample.minicellFloreser.sp <- sp::SpatialPointsDataFrame(coords = data.frame(sample.minicellFloreser[, 1:2]),
                                                               data = data.frame(sample.minicellFloreser[, 1:2]),
                                                               proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))


# extract pixel area raster data for sample points
aux.pixelArea <- raster::extract(clean.floreserArea, sample.minicellFloreser.sp)

# clean environment
rm(sample.minicellFloreser.sp)

# merge pixel area variable with sample 2000
sample.minicellFloreser$pixel_area <- aux.pixelArea

# change final object name
pixelArea.minicellFloreser <- sample.minicellFloreser

# clear environmnet
rm(sample.minicellFloreser, aux.pixelArea)



# CLEAN TEMP DIR
# showTmpFiles()
raster::removeTmpFiles(h = 2)
gc()





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(pixelArea.minicellFloreser$pixel_area) <- "pixel area (ha)"



# POST-TREATMENT OVERVIEW
# summary(pixelArea.minicellFloreser)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(pixelArea.minicellFloreser,
     file = file.path("data/projectSpecific/minicellFloreser",
                      paste0("pixelArea_minicellFloreser", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/minicellFloreser")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
