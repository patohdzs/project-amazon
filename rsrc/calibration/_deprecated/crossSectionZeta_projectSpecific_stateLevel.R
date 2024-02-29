
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: COMBINE VARIABLES RELEVANT FOR ZETA CALIBRATION
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "crossSectionZeta_projectSpecific_stateLevel script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns
library(sf) # manipulate spatial data





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# SAMPLE SPATIAL MUNI LEVEL
load(file.path("data/projectSpecific/muniLevel/sampleSpatial_muniLevel.Rdata"))



# CONSUMPTION EXPENDITURE
load(file.path("data/raw2clean/pof_ibge/output",
               "clean_pof.Rdata"))


# DEFLATOR (IPA-EP-DI)
load("data/raw2clean/deflatorIPA_fgv/output/clean_deflatorIPA.Rdata")





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

aux.01.2017 <- clean.deflatorIPA[clean.deflatorIPA$date == "2017-01-01",]$deflator_ipa

aux.deflator <-
  clean.deflatorIPA %>%
  dplyr::mutate(deflator_ipa = deflator_ipa/aux.01.2017) %>% # change base year to january 2017
  dplyr::filter(date >= '2017-07-01' & date <= '2018-07-01') %>%
  dplyr::pull(deflator_ipa) %>%
  mean()

# MERGE ALL DATASETS
crossSection.zeta <-
  sampleSpatial.muniLevel %>%
  dplyr::mutate(UF = as.numeric(str_sub(muni_code, 1, 2))) %>%
  dplyr::filter(!is.na(biomeAmazon_share)) %>% # remove municipalities outside amazon biome
  dplyr::group_by(UF) %>%
  dplyr::summarise(biomeAmazon_share = weighted.mean(biomeAmazon_share, muni_area)) %>%
  dplyr::left_join(clean.pof, by = c("UF")) %>%
  dplyr::mutate(pof_consumptionExpenditure_livestockRelated = pof_consumptionExpenditure_livestockRelated/aux.deflator) %>%  # inflation adjustment jan/2017 BRL
  dplyr::mutate(pof_consumptionExpenditure_livestockRelated = pof_consumptionExpenditure_livestockRelated/3.634267) # change from BRL to USD (commercial exchange rate - selling - average - monthly - 07/2017-07/2018 - ipeadata))

# clear environment
rm(aux.01.2017, aux.deflator, sampleSpatial.muniLevel, clean.pof)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(crossSection.zeta$UF) <- "state code identifier"
sjlabelled::set_label(crossSection.zeta$biomeAmazon_share) <- "fraction of the state inside the Amazon Biome"
sjlabelled::set_label(crossSection.zeta$pof_consumptionExpenditure_livestockRelated) <- "total annual consumption expediture (USD 2017) of families with at least one member working in a livestock related job per state (UF)"




# POST-TREATMENT OVERVIEW
# summary(crossSection.zeta)
# View(crossSection.zeta)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(crossSection.zeta,
     file = file.path("data/projectSpecific/stateLevel",
                      paste0("crossSection_zeta_stateLevel", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/stateLevel")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
