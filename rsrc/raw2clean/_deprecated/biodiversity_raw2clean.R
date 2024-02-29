
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW DATA SHAPEFILE OF BIODIVERSITY (2007 AND 2018)
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
here::i_am("code/raw2clean/biodiversity_raw2clean.R", uuid = "9e480471-3c96-4ac9-81e8-690223817d27")


# START TIME
tictoc::tic(msg = "biodiversity_raw2clean script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sf", groundhog.date)  # manipulate vector data
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read shapefile
raw.biodiversity.2007 <- sf::st_read(dsn   = here::here("data/raw2clean/biodiversity_imazon/input/2007"),
                           layer = "AP2007_Areas_Prioritarias")

raw.biodiversity.2018 <- sf::st_read(dsn   = here::here("data/raw2clean/biodiversity_imazon/input/2018"),
                                     layer = "AP2018_Areas_Prioritarias_Amazonia")



# DATA EXPLORATION [disabled for speed]
# summary(raw.biodiversity.2007)
# View(raw.biodiversity.2007)
# plot(raw.biodiversity.2007$geometry)
# summary(raw.biodiversity.2018)
# View(raw.biodiversity.2018)
# plot(raw.biodiversity.2018$geometry)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# COLUMN CLEANUP
# names
colnames(raw.biodiversity.2007)
colnames(raw.biodiversity.2018)



# SELECT AMAZON BIOME
raw.biodiversity.2007 <- raw.biodiversity.2007 %>% dplyr::filter(BIOMA == "Am")


# DROP Z RANGE
raw.biodiversity.2018 <- sf::st_zm(raw.biodiversity.2018)



# PROJECTION
raw.biodiversity.2007 <- sf::st_transform(x = raw.biodiversity.2007, crs = 5880) # SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)
raw.biodiversity.2018 <- sf::st_transform(x = raw.biodiversity.2018, crs = 5880) # SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)



# GEOMETRY CLEANUP
raw.biodiversity.2007 <- sf::st_make_valid(raw.biodiversity.2007)
raw.biodiversity.2018 <- sf::st_make_valid(raw.biodiversity.2018)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS



# change object name for exportation
clean.biodiversity.2007 <- raw.biodiversity.2007
clean.biodiversity.2018 <- raw.biodiversity.2018



# POST-TREATMENT OVERVIEW
# summary(clean.biodiversity.2007)
# View(clean.biodiversity.2007)
# plot(clean.biodiversity.2007$geometry)
# summary(clean.biodiversity.2018)
# View(clean.biodiversity.2018)
# plot(clean.biodiversity.2018$geometry)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.biodiversity.2007,
     file = here::here("data/raw2clean/biodiversity_imazon/output",
                      "clean_biodiversity_2007.Rdata"))

save(clean.biodiversity.2018,
     file = here::here("data/raw2clean/biodiversity_imazon/output",
                       "clean_biodiversity_2018.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------