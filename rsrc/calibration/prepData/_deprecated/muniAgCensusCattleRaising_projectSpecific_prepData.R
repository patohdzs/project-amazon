
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: COMBINE VARIABLES RELEVANT FOR CATTLE RAISING WORKERS PER HECTARE CALIBRATION AT THE MUNI LEVEL
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
here::i_am("code/projectSpecific/prepData/muniAgCensusCattleRaising_projectSpecific_prepData.R", uuid = "a0a69fcb-bf58-4f1f-a0d5-41283814f248")


# START TIME
tictoc::tic(msg = "muniAgCensusCattleRaising_projectSpecific_prepData script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sf", groundhog.date)  # manipulate spatial data (vector format)
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# INPUT VARIABLES TO CALCULATE THETA
load(here::here("data/projectSpecific/prepData/muniTheta_prepData.Rdata"))


# AG CENSUS CATTLE RAISING FARMS VARIABLES
load(here::here("data/raw2clean/agCensusCattleRaising_ibge/output", "clean_agCensusCattleRaising.Rdata"))




# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# MERGE ALL DATASETS
muniAgCensusCattleRaising.prepData <-
  muniTheta.prepData %>%
  dplyr::left_join(clean.agCensusCattleRaising %>% dplyr::rename_with(.fn = ~paste0(., "_cattleRaising"), .cols = -muni_code), by = c("muni_code")) %>%
  dplyr::mutate(wage_expenditure_cattleRaising  = 1000*wage_expenditure_cattleRaising /3.192) %>%  # change from thousand BRL to BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
  dplyr::mutate(workers_ha = dplyr::if_else(pastureArea_value_cattleRaising == 0 & !is.na(workers_number_cattleRaising), 0, workers_number_cattleRaising/pastureArea_value_cattleRaising),
                wage_ha = dplyr::if_else(pastureArea_value_cattleRaising == 0 & !is.na(wage_expenditure_cattleRaising), 0, wage_expenditure_cattleRaising/pastureArea_value_cattleRaising),
                wage_worker = dplyr::if_else(workers_number_cattleRaising == 0 & !is.na(wage_expenditure_cattleRaising), 0, wage_expenditure_cattleRaising/workers_number_cattleRaising),
                addedValue = if_else(cattleSlaughter_value_ha == 0 | is.na(cattleSlaughter_value_ha),
                                     as.numeric(NA),
                                     (cattleSlaughter_value_ha - wage_ha)/cattleSlaughter_value_ha),
                addedValue = if_else(addedValue < 0, 0, addedValue)) %>%
  dplyr::select(muni_code, muni_area, biomeAmazon_share,
                wage_expenditure_cattleRaising, workers_ha,  wage_worker, workers_number_cattleRaising,
                wage_ha, cattleSlaughter_value_ha, pastureArea_value_cattleRaising, pastureArea_value, addedValue,
                lon, lat, historical_precip, historical_temp, geometry)

summary(muniAgCensusCattleRaising.prepData)
quantile(muniAgCensusCattleRaising.prepData$addedValue, na.rm = T, seq(0,1,0.01))
weighted.mean(muniAgCensusCattleRaising.prepData$addedValue, w = muniAgCensusCattleRaising.prepData$pastureArea_value, na.rm = TRUE)
weighted.mean(muniAgCensusCattleRaising.prepData$addedValue, w = muniAgCensusCattleRaising.prepData$pastureArea_value_cattleRaising, na.rm = TRUE)

# clear environment
rm(clean.agCensusCattleRaising, muniTheta.prepData)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(muniAgCensusCattleRaising.prepData$lon) <- "longitude of municipality centroid (calculate under EPSG:5880)"
sjlabelled::set_label(muniAgCensusCattleRaising.prepData$lat) <- "latitude  of municipality centroid (calculate under EPSG:5880)"
sjlabelled::set_label(muniAgCensusCattleRaising.prepData$historical_precip) <- "historical (1970-2000) total annual precipitation (mm)"
sjlabelled::set_label(muniAgCensusCattleRaising.prepData$historical_temp) <- "historical (1970-2000) mean annual average temperature (celsius degrees) "
sjlabelled::set_label(muniAgCensusCattleRaising.prepData$workers_ha) <- "number of workers per pasture area (people/ha)"
sjlabelled::set_label(muniAgCensusCattleRaising.prepData$wage_worker) <- "wage expenditure per worker (USD/people)"



# POST-TREATMENT OVERVIEW
# summary(muniAgCensusCattleRaising.prepData)
# View(muniAgCensusCattleRaising.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(muniAgCensusCattleRaising.prepData,
     file = here::here("data/projectSpecific/prepData",
                      paste0("muniAgCensusCattleRaising_prepData", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/prepData")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
