
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: ADD NUMBER OF PIXELS (FROM THE PIXEL SAMPLE) IN EACH MUNI
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
here::i_am("code/projectSpecific/prepData/muniPixelCount_projectSpecific_prepData.R", uuid = "feeeefe6-67dd-4e8a-8b8e-73d8f33e4a4e")


# START TIME
tictoc::tic(msg = "muniPixelCount_projectSpecific_prepData script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sf", groundhog.date)  # manipulate spatial data (vector format)
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# SPATIAL MUNI SAMPLE
load(here::here("data/projectSpecific/prepData/sampleMuniSpatial_prepData.Rdata"))



# SPATIAL MINICELL SAMPLE
load(here::here("data/projectSpecific/prepData/samplePixel_prepData.Rdata"))





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# DATA PREP

# select cross-section of points
samplePixel.prepData <- samplePixel.prepData %>% dplyr::group_by(lon, lat) %>% dplyr::slice(1) %>% dplyr::select(lon, lat)

# transform to spatial object
samplePixel.prepData <- sf::st_as_sf(coords = c("lon", "lat"),
                                         x = samplePixel.prepData,
                                         crs = sf::st_crs("+proj=longlat +datum=WGS84 +no_defs"))

# adjust projection
samplePixel.prepData <- sf::st_transform(samplePixel.prepData, sf::st_crs(sampleMuniSpatial.prepData))

# add number of pixels in each muni
sampleMuniSpatial.prepData$count_pixels <- lengths(sf::st_intersects(sampleMuniSpatial.prepData, samplePixel.prepData))

# change final object name
muniPixelCount.prepData <- sampleMuniSpatial.prepData

# clear environment
rm(sampleMuniSpatial.prepData, samplePixel.prepData)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(muniPixelCount.prepData$count_pixels) <- "number of sampled pixels inside the muni"



# POST-TREATMENT OVERVIEW
# summary(muniPixelCount.prepData)
# View(muniPixelCount.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(muniPixelCount.prepData,
     file = file.path("data/projectSpecific/prepData",
                      paste0("muniPixelCount_prepData", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/prepData")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
