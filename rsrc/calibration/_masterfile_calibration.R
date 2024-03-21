
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




# 1055 SITES MODEL -------------------------------------------------------------------------------------------------------------------------------

# GENERATE AGGREGATED SAMPLE OF INTEREST (DIVIDE AMAZON BIOME INTO 1055 CELLS)
source(here::here("rsrc/calibration/amazonBiome_1043SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# CALCULATE SITE AREA - 1055 SITES
source(here::here("rsrc/calibration/pixelArea_1043SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# GENERATE AGGREGATED MAPBIOMAS VARIABLES (FOREST, AGRICULTURAL USE, OTHER) - 1055 SITES
source(here::here("rsrc/calibration/mapbiomas_1043SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# PARAMETERS CALIBRATION (1055 SITES MODEL)
source(here::here("rsrc/hmc_calibration/hmc_1043SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())





# 81 SITES MODEL -------------------------------------------------------------------------------------------------------------------------------

# PARAMETERS CALIBRATION (81 SITES MODEL)
source(here::here("rsrc/hmc_calibration/hmc_78SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())





# 41 SITES MODEL -------------------------------------------------------------------------------------------------------------------------------

# PARAMETERS CALIBRATION (41 SITES MODEL)
source(here::here("rsrc/hmc_calibration/hmc_40SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())






# 25 SITES MODEL -------------------------------------------------------------------------------------------------------------------------------

# PARAMETERS CALIBRATION (25 SITES MODEL)
source(here::here("rsrc/hmc_calibration/hmc_24SitesModel.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())





# GLOBAL MODEL ---------------------------------------------------------------------------------------------------------------------------------------

# PARAMETERS CALIBRATION AND INITIAL CONDITIONS
source(here::here("rsrc/hmc_calibration/hmc_GlobalSitesModel.R"), encoding = "UTF-8", echo = T)

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