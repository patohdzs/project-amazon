
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


# START TIMER
tictoc::tic(msg = "_masterfile_prep.R script", log = TRUE)


if (!dir.exists("data/prepData")) {
    dir.create("data/prepData", recursive = TRUE)
}



# PREP DATA ------------------------------------------------------------------------------------------------------------------------------------------

# SERIES

# CONSTRUCT MONTHLY COMMODITY REAL PRICES INDICES
source(here::here("rsrc/prepData/seriesPriceCattle_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())



# PIXEL LEVEL

# EXTRACT RANDOM SAMPLE OF MAPBIOMAS 30M-PIXELS AND RECOVER FULL PANEL (1985-2019)
source(here::here("rsrc/prepData/sampleConstructionPixel_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# CALCULATE AREA OF MAPBIOMAS 30M-PIXELS
source(here::here("rsrc/prepData/pixelArea_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# CREATE AGGREGATED CATEGORIES OF INTEREST BASED ON MAPBIOMAS 30M-PIXELS VALUES
source(here::here("rsrc/prepData/pixelCategories_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# ADD 2010 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
source(here::here("rsrc/prepData/pixelBiomass2010_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# ADD 2017 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
source(here::here("rsrc/prepData/pixelBiomass2017_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# ADD 2018 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
source(here::here("rsrc/prepData/pixelBiomass2018_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())



# MUNI LEVEL

# DEFINE MUNI-LEVEL SAMPLE AND SAVE IT IN SPATIAL, CROSS-SECTION AND PANEL FORMATS
source(here::here("rsrc/prepData/sampleConstructionMuni_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# Prepare gamma muni
source(here::here("rsrc/prepData/merge_muni_gamma.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# COMBINE VARIABLES RELEVANT FOR THETA CALIBRATION AT THE MUNI LEVEL
source(here::here("rsrc/prepData/muniTheta_prepData_gamma.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# STATE LEVEL

# PREPATE DATA TO ESTIMATE PARAMETER K (EMISSION FACTOR OF AGRICULTURAL SECTOR)
source("rsrc/prepData/stateEmission_prepData.R", encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())





# EXPORT TIME PROCESSING -----------------------------------------------------------------------------------------------------------------------------

# END TIMER
tictoc::toc(log = TRUE)




# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------