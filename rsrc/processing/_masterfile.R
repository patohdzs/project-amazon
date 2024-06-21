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


if (!dir.exists("data/processed")) {
  dir.create("data/processed", recursive = TRUE)
}


# CONSTRUCT MONTHLY COMMODITY REAL PRICES INDICES
source(here::here("rsrc/processing/prep_cattle_price_index.R"), encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())



# PIXEL LEVEL

# EXTRACT RANDOM SAMPLE OF MAPBIOMAS 30M-PIXELS AND RECOVER FULL PANEL (1985-2019)
source(here::here("rsrc/processing/prep_pixel_sample.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# CALCULATE AREA OF MAPBIOMAS 30M-PIXELS
source(here::here("rsrc/processing/prep_pixel_area.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# CREATE AGGREGATED CATEGORIES OF INTEREST BASED ON MAPBIOMAS 30M-PIXELS VALUES
source(here::here("rsrc/processing/prep_pixel_categories.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())



# ADD 2017 ABOVEGROUND BIOMASS DATA (ESA) TO MAPBIOMAS 30M-PIXELS
source(here::here("rsrc/processing/prep_pixel_biomass.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())




# MUNI LEVEL

# DEFINE MUNI-LEVEL SAMPLE AND SAVE IT IN SPATIAL, CROSS-SECTION AND PANEL FORMATS
source(here::here("rsrc/processing/prep_muni_sample.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# Prepare gamma muni
source(here::here("rsrc/processing/prep_muni_biomass.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())


# COMBINE VARIABLES RELEVANT FOR THETA CALIBRATION AT THE MUNI LEVEL
source(here::here("rsrc/processing/prep_muni_data.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# STATE LEVEL

# PREPATE DATA TO ESTIMATE PARAMETER K (EMISSION FACTOR OF AGRICULTURAL SECTOR)
source("rsrc/processing/prep_state_emissions.R", encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# PREPATE GENERATE AGGREGATED MAPBIOMAS VARIABLES (FOREST, AGRICULTURAL USE, OTHER) - 1055 SITES
source("rsrc/processing/prep_land_use_rasters.R", encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# PREPATE GENERATE AGGREGATED SAMPLE OF INTEREST (DIVIDE AMAZON BIOME INTO 1055 CELLS)
source("rsrc/processing/prep_biome_rasters.R", encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())



# EXPORT TIME PROCESSING -----------------------------------------------------------------------------------------------------------------------------

# END TIMER
tictoc::toc(log = TRUE)




# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
