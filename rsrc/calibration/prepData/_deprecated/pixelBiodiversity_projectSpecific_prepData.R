
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: ADD DUMMY IF THE PIXEL IS INSIDE A PRIORITY AREA FOR BIODIVERSITY
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
here::i_am("code/projectSpecific/prepData/pixelBiodiversity_projectSpecific_prepData.R", uuid = "feeeefe6-67dd-4e8a-8b8e-73d8f33e4a4e")


# START TIME
tictoc::tic(msg = "pixelBiodiversity_projectSpecific_prepData script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sf", groundhog.date)  # manipulate spatial data (vector format)
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# SPATIAL BIOVERSITY PRIORITY AREAS
load(here::here("data/raw2clean/biodiversity_imazon/output",
                "clean_biodiversity_2007.Rdata"))
load(here::here("data/raw2clean/biodiversity_imazon/output",
                "clean_biodiversity_2018.Rdata"))



# SPATIAL MINICELL SAMPLE
load(here::here("data/projectSpecific/prepData/samplePixel_prepData.Rdata"))





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# DATA PREP

# select cross-section of points
samplePixel.prepData <- samplePixel.prepData %>% dplyr::group_by(lon, lat) %>% dplyr::slice(1) %>% dplyr::select(lon, lat)

# transform to spatial object
samplePixel.prepData <- sf::st_as_sf(coords = c("lon", "lat"),
                                     x = samplePixel.prepData,
                                     crs = sf::st_crs("+proj=longlat +datum=WGS84 +no_defs"),
                                     remove = FALSE)

# adjust projection
samplePixel.prepData <- sf::st_transform(samplePixel.prepData, sf::st_crs(clean.biodiversity.2007))

# add dummy =1 if pixel is inside a priority area for biodiversitt
samplePixel.prepData <-
  samplePixel.prepData %>%
  dplyr::mutate(d_biodiversity2007 = if_else(lengths(sf::st_intersects(samplePixel.prepData, clean.biodiversity.2007)) > 0, 1 , 0),
                d_biodiversity2018 = if_else(lengths(sf::st_intersects(samplePixel.prepData, clean.biodiversity.2018)) > 0, 1 , 0))

# change final object name
pixelBiodiversity.prepData <- samplePixel.prepData %>% sf::st_drop_geometry()

# clear environment
rm(clean.biodiversity.2007, clean.biodiversity.2018, samplePixel.prepData)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(pixelBiodiversity.prepData$d_biodiversity2007) <- "=1 if pixel is inside a 2007 priority area for biodiversity"
sjlabelled::set_label(pixelBiodiversity.prepData$d_biodiversity2018) <- "=1 if pixel is inside a 2018 priority area for biodiversity"



# POST-TREATMENT OVERVIEW
# summary(pixelBiodiversity.prepData)
# View(pixelBiodiversity.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(pixelBiodiversity.prepData,
     file = file.path("data/projectSpecific/prepData",
                      paste0("pixelBiodiversity_prepData", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/prepData")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
