
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CONSTRUCT MONTHLY COMMODITY REAL PRICES INDICES
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "pricesMonthly_projectSpecific_series script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns, prefer than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# AGRICUTLURAL COMMODITY PRICES
load("data/raw2clean/commodityPrices_seabpr/output/clean_commodityPrices.Rdata")



# DEFLATOR (IPA-EP-DI)
load("data/raw2clean/deflatorIPA_fgv/output/clean_deflatorIPA.Rdata")





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# TRANSFORM NOMINAL MONTHLY COMMODITY PRICES TO REAL STANDARDIZED PRICES
aux.realPrices <-
  clean.commodityPrices %>%
  dplyr::left_join(clean.deflatorIPA) %>% # merge commodity nominal prices with deflator index
  dplyr::mutate(price_r_rice      = (price_rice/deflator_ipa)    * (1000/50)) %>% # transform to real price, standardize unit to 1t
  dplyr::mutate(price_r_cassava   = price_cassava/deflator_ipa              ) %>% # transform to real price
  dplyr::mutate(price_r_corn      = (price_corn/deflator_ipa)    * (1000/60)) %>% # transform to real price, standardize unit  to 1t
  dplyr::mutate(price_r_soybean   = (price_soybean/deflator_ipa) * (1000/60)) %>% # transform to real price, standardize unit to 1t
  dplyr::mutate(price_r_sugarcane = price_sugarcane/deflator_ipa            ) %>% # transform to real price
  dplyr::mutate(price_r_cattle    = price_cattle/deflator_ipa               ) # transform to real price





# TRANSFORM REAL PRICES TO MONTHLY INDICES (BASELINE 2000)

# extract baseline for 1990-04
aux.baseline.1 <- aux.realPrices[4, ]


# construct annual indices
pricesMonthly.series <-
  aux.realPrices %>%
  dplyr::mutate(price_index_cattle      = (price_r_cattle/aux.baseline.1$price_r_cattle)) %>%
  dplyr::mutate(price_index_sugarcane   = (price_r_sugarcane/aux.baseline.1$price_r_sugarcane)) %>%
  dplyr::mutate(price_index_corn        = (price_r_corn/aux.baseline.1$price_r_corn)) %>%
  dplyr::mutate(price_index_soybean     = (price_r_soybean/aux.baseline.1$price_r_soybean)) %>%
  dplyr::mutate(price_index_cassava     = (price_r_cassava/aux.baseline.1$price_r_cassava)) %>%
  dplyr::mutate(price_index_rice        = (price_r_rice/aux.baseline.1$price_r_rice)) %>%
  dplyr::select(-starts_with("price_r"))





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(pricesMonthly.series$date)                     <- "date of reference"
sjlabelled::set_label(pricesMonthly.series$price_index_cassava)      <- "real price index, cassava (monthly, 2000-01-01 = 1)"
sjlabelled::set_label(pricesMonthly.series$price_index_sugarcane)    <- "real price index, sugarcane (monthly, 2000-01-01 = 1)"
sjlabelled::set_label(pricesMonthly.series$price_index_cattle)       <- "real price index, cattle (monthly, 2000-01-01 = 1)"
sjlabelled::set_label(pricesMonthly.series$price_index_corn)         <- "real price index, corn (monthly, 2000-01-01 = 1)"
sjlabelled::set_label(pricesMonthly.series$price_index_soybean)      <- "real price index, soybean (monthly, 2000-01-01 = 1)"
sjlabelled::set_label(pricesMonthly.series$price_index_rice)         <- "real price index, rice (monthly, 2000-01-01 = 1)"



# POST-TREATMENT OVERVIEW
# summary(pricesMonthly.series)
# View(pricesMonthly.series)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(pricesMonthly.series,
     file = file.path("data/projectSpecific/series",
                      paste0("pricesMonthly_series", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/series")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
