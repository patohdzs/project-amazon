
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
source("code/setup.R")


# START TIMER
tictoc::tic(msg = "_masterfile_calibration.R script", log = T)





# PREP DATA ------------------------------------------------------------------------------------------------------------------------------------------

# SERIES

# CONSTRUCT MONTHLY COMMODITY REAL PRICES INDICES
source(here::here("code/calibration/prepData/seriesPriceCattle_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())



# PIXEL LEVEL

# EXTRACT RANDOM SAMPLE OF MAPBIOMAS 30M-PIXELS AND RECOVER FULL PANEL (1985-2019)
source(here::here("code/calibration/prepData/sampleConstructionPixel_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# CALCULATE AREA OF MAPBIOMAS 30M-PIXELS
source(here::here("code/calibration/prepData/pixelArea_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# CREATE AGGREGATED CATEGORIES OF INTEREST BASED ON MAPBIOMAS 30M-PIXELS VALUES
source(here::here("code/calibration/prepData/pixelCategories_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# ADD 2010 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
source(here::here("code/calibration/prepData/pixelBiomass2010_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# ADD 2017 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
source(here::here("code/calibration/prepData/pixelBiomass2017_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# ADD 2018 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
source(here::here("code/calibration/prepData/pixelBiomass2018_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())



# MUNI LEVEL

# DEFINE MUNI-LEVEL SAMPLE AND SAVE IT IN SPATIAL, CROSS-SECTION AND PANEL FORMATS
source(here::here("code/calibration/prepData/sampleConstructionMuni_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# COMBINE VARIABLES RELEVANT FOR THETA CALIBRATION AT THE MUNI LEVEL
source(here::here("code/calibration/prepData/muniTheta_prepData.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())



# STATE LEVEL

# PREPATE DATA TO ESTIMATE PARAMETER K (EMISSION FACTOR OF AGRICULTURAL SECTOR)
source("code/calibration/prepData/stateEmission_prepData.R", encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())





# 1055 SITES MODEL -------------------------------------------------------------------------------------------------------------------------------

# GENERATE AGGREGATED SAMPLE OF INTEREST (DIVIDE AMAZON BIOME INTO 1055 CELLS)
source(here::here("code/calibration/amazonBiome_1055SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# CALCULATE SITE AREA - 1055 SITES
source(here::here("code/calibration/pixelArea_1055SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# GENERATE AGGREGATED MAPBIOMAS VARIABLES (FOREST, AGRICULTURAL USE, OTHER) - 1055 SITES
source(here::here("code/calibration/mapbiomas_1055SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# PARAMETERS CALIBRATION (1055 SITES MODEL)
source(here::here("code/calibration/calibration_1055SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())





# 81 SITES MODEL -------------------------------------------------------------------------------------------------------------------------------

# PARAMETERS CALIBRATION (81 SITES MODEL)
source(here::here("code/calibration/calibration_81SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())





# 41 SITES MODEL -------------------------------------------------------------------------------------------------------------------------------

# PARAMETERS CALIBRATION (41 SITES MODEL)
source(here::here("code/calibration/calibration_41SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())






# 25 SITES MODEL -------------------------------------------------------------------------------------------------------------------------------

# PARAMETERS CALIBRATION (25 SITES MODEL)
source(here::here("code/calibration/calibration_25SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())





# GLOBAL MODEL ---------------------------------------------------------------------------------------------------------------------------------------

# PARAMETERS CALIBRATION AND INITIAL CONDITIONS
source(here::here("code/calibration/calibration_globalModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())





# EXPORT TIME PROCESSING -----------------------------------------------------------------------------------------------------------------------------

# END TIMER
tictoc::toc(log = T)


# SOURCE FUNCTIONS - after scripts to avoid rm(list = ls())
source(here::here("code/_functions/ExportTimeProcessing.R"))

# export time to csv table
ExportTimeProcessing("code/calibration")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------