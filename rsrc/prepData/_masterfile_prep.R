
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: MASTERFILE SCRIPT TO SOURCE ALL CALIBRATION SCRIPTS
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "_masterfile_calibration.R script", log = T)





# PREP DATA ------------------------------------------------------------------------------------------------------------------------------------------

# SERIES

# CONSTRUCT MONTHLY COMMODITY REAL PRICES INDICES
source(here::here("rsrc/calibration/prepData/seriesPriceCattle_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())



# PIXEL LEVEL

# EXTRACT RANDOM SAMPLE OF MAPBIOMAS 30M-PIXELS AND RECOVER FULL PANEL (1985-2019)
source(here::here("rsrc/calibration/prepData/sampleConstructionPixel_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# CALCULATE AREA OF MAPBIOMAS 30M-PIXELS
source(here::here("rsrc/calibration/prepData/pixelArea_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# CREATE AGGREGATED CATEGORIES OF INTEREST BASED ON MAPBIOMAS 30M-PIXELS VALUES
source(here::here("rsrc/calibration/prepData/pixelCategories_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# ADD 2010 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
source(here::here("rsrc/calibration/prepData/pixelBiomass2010_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# ADD 2017 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
source(here::here("rsrc/calibration/prepData/pixelBiomass2017_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# ADD 2018 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
source(here::here("rsrc/calibration/prepData/pixelBiomass2018_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())



# MUNI LEVEL

# DEFINE MUNI-LEVEL SAMPLE AND SAVE IT IN SPATIAL, CROSS-SECTION AND PANEL FORMATS
source(here::here("rsrc/calibration/prepData/sampleConstructionMuni_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# Prepare gamma muni
source(here::here("rsrc/calibration/prepData/merge_muni_gamma.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# COMBINE VARIABLES RELEVANT FOR THETA CALIBRATION AT THE MUNI LEVEL
source(here::here("rsrc/calibration/prepData/muniTheta_prepData_gamma.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# STATE LEVEL

# PREPATE DATA TO ESTIMATE PARAMETER K (EMISSION FACTOR OF AGRICULTURAL SECTOR)
source("rsrc/calibration/prepData/stateEmission_prepData.R", encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())





# EXPORT TIME PROCESSING -----------------------------------------------------------------------------------------------------------------------------

# END TIMER
tictoc::toc(log = T)


# # SOURCE FUNCTIONS - after scripts to avoid rm(list = ls())
# source(here::here("rsrc/_functions/ExportTimeProcessing.R"))

# # export time to csv table
# ExportTimeProcessing("rsrc/calibration")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------