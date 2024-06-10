# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: MASTERFILE SCRIPT TO SOURCE ALL RAW2CLEAN SCRIPTS
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -



# START TIMER
tictoc::tic(msg = "_masterfile_raw2clean.R script", log = TRUE)

# Create folder
if (!dir.exists("data/clean")) {
  dir.create("data/clean", recursive = TRUE)
}
if (!dir.exists("data/calibration")) {
  dir.create("data/calibration", recursive = TRUE)
}
if (!dir.exists("data/calibration/hmc")) {
  dir.create("data/calibration/hmc", recursive = TRUE)
}

# SOURCE
# TREAT RAW DATA - AMAZON BIOME BOUNDARY (IBGE - 2019)
source(here::here("rsrc/raw2clean/amazonBiome_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# CLEAN RAW DATA SHAPEFILE OF MUNI DIVISION (2015)
source(here::here("rsrc/raw2clean/muniDivision2015_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# CLEAN RAW AGRICULTURAL COMMODITY PRICES FROM THE PARANA SECRETARIAT OF SUPPLY AND AGRICULTURE (SEAB-PR)
source(here::here("rsrc/raw2clean/commodityPrices_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# CLEAN RAW DEFLATOR IPA-DI - FGV
source(here::here("rsrc/raw2clean/deflatorIPA_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# CLEAN RAW CATTLE SOLD - AGRICULTURAL CENSUS 2017 (IBGE)
source(here::here("rsrc/raw2clean/agCensus2017CattleSold_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# CLEAN RAW AGRICULTURAL USE AREA - AGRICULTURAL CENSUS 2017 (IBGE)
source(here::here("rsrc/raw2clean/agCensus2017AgUseArea_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# CLEAN RAW CATTLE SOLD FOR SLAUGHTER - AGRICULTURAL CENSUS 2006 (IBGE)
source(here::here("rsrc/raw2clean/agCensus2006CattleSlaughter_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# CLEAN RAW AGRICULTURAL USE AREA - AGRICULTURAL CENSUS 2006 (IBGE)
source(here::here("rsrc/raw2clean/agCensus2006AgUseArea_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# CLEAN RAW SEEG AGRICULTURAL EMISSION DATA - LEGAL AMAZON STATES 1990-2019
source(here::here("rsrc/raw2clean/emission_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# TREAT RAW DATA EMISSIONS AND GDP BY COUNTRY (WORLD BANK)
source(here::here("rsrc/raw2clean/emissionKuznets_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# TREAT RAW DATA CARBON PRICE (WORLD BANK)
source(here::here("rsrc/raw2clean/carbonPrice_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# CLEAN RAW BRAZILIAN AMAZON LAND USE DATA TO PROVIDE MUNICIPALITY-BY-YEAR PANEL (MAPBIOMAS COLLECTION 5)
source(here::here("rsrc/raw2clean/landUseCoverMuni_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# UNIFY MONTHLY TIF FILES - HISTORICAL TEMPERATURE (WORLD CLIM)
source(here::here("rsrc/raw2clean/temperature_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# UNIFY MONTHLY TIF FILES - HISTORICAL PRECIPITATION (WORLD CLIM)
source(here::here("rsrc/raw2clean/precipitation_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# RECLASSIFY PIXELS OUTSIDE BIOME BOUNDARY - LAND USE AND COVER (MAPBIOMAS - COL5)
source(here::here("rsrc/raw2clean/landUseCover_raw2clean.R"), encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

# DOWNLOAD ABOVERGROUND BIOMASS/CARBON DATA (ESA BIOMASS - 2010, 2017, 2018)
# DOWNLOAD PROCESS IS OPTIONAL GIVEN THAT THE DATA IS PROVIDED
# source(here::here("rsrc/raw2clean/abovegroundBiomassESA_download.R"), encoding = "UTF-8", echo = T)


# # EXPORT TIME PROCESSING

# END TIMER
tictoc::toc(log = TRUE)


# END OF SCRIPT
