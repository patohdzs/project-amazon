
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: COMBINE VARIABLES RELEVANT FOR CATTLE SLAUGHTER BY HECTARE ANALYSIS
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "cattleSlaughter_projectSpecific_muniLevel script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns




# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# SAMPLE CROSS SECTION MUNI LEVEL
load(file.path("data/projectSpecific/muniLevel/sampleCrossSection_muniLevel.Rdata"))


# AG CENSUS CATTLE FOR SLAUGHTER
load(file.path("data/raw2clean/agCensusCattleSlaughter_ibge/output", "clean_agCensusCattleSlaughter.Rdata"))


# AG CENSUS PASTURE AREA
load(file.path("data/raw2clean/agCensusPastureArea_ibge/output", "clean_agCensusPastureArea.Rdata"))





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# MERGE ALL DATASETS
crossSection.cattleSlaughter <-
  clean.agCensusPastureArea %>%
  dplyr::left_join(clean.agCensusCattleSlaughter, by = c("muni_code")) %>%
  dplyr::left_join(sampleCrossSection.muniLevel, by = c("muni_code")) %>%
  dplyr::select(muni_code, pastureArea_value, cattleSlaughter_head, muni_area) %>%
  dplyr::mutate(d_sample = if_else(is.na(muni_area), 0, 1)) %>% # create dummy =1 if muni is in the sample and 0 otherwise
  dplyr::mutate(kg_carcass_per_ha = if_else(pastureArea_value  == 0, as.numeric(NA), 225*cattleSlaughter_head/pastureArea_value))



# clear environment
rm(clean.agCensusCattleSlaughter, clean.agCensusPastureArea, sampleCrossSection.muniLevel)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(crossSection.cattleSlaughter$d_sample) <- "=1 if muni is in the sample and 0 otherwise"
sjlabelled::set_label(crossSection.cattleSlaughter$kg_carcass_per_ha) <- "total carcass weigth of cattle sold for slaughter"



# POST-TREATMENT OVERVIEW
# summary(crossSection.cattleSlaughter)
# View(crossSection.cattleSlaughter)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(crossSection.cattleSlaughter,
     file = file.path("data/projectSpecific/muniLevel",
                      paste0("crossSection_cattleSlaughter_muniLevel", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/muniLevel")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
