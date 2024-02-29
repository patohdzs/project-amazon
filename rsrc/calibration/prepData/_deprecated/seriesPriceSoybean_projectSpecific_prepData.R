
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
here::i_am("code/projectSpecific/prepData/seriesPriceSoybean_projectSpecific_prepData.R", uuid = "25d56594-ebdd-4f66-a2c5-94645d681c84")


# START TIME
tictoc::tic(msg = "seriesPriceSoybean_projectSpecific_prepData script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# AGRICUTLURAL COMMODITY PRICES
load(here::here("data/raw2clean/commodityPrices_seabpr/output/clean_commodityPrices.Rdata"))



# DEFLATOR (IPA-EP-DI)
load(here::here("data/raw2clean/deflatorIPA_fgv/output/clean_deflatorIPA.Rdata"))





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

aux.01.2017 <- clean.deflatorIPA[clean.deflatorIPA$date == "2017-01-01",]$deflator_ipa

# TRANSFORM NOMINAL MONTHLY COMMODITY PRICES TO REAL STANDARDIZED PRICES
seriesPriceSoybean.prepData <-
  clean.commodityPrices %>%
  dplyr::left_join(clean.deflatorIPA) %>% # merge commodity nominal prices with deflator index
  dplyr::mutate(deflator_ipa = deflator_ipa/aux.01.2017) %>% # change base year to january 2017
  dplyr::mutate(price_soybean = price_soybean/60) %>% # trasform unit to price per kg
  dplyr::mutate(price_real_mon_soybean = price_soybean/deflator_ipa) %>% # transform to real price (constant january 2017)
  dplyr::mutate(month = lubridate::month(date),
                year = lubridate::year(date)) %>% # construct year, month and trimester variables
  dplyr::select(date, year, month, starts_with("price_real_"))





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(seriesPriceSoybean.prepData$date)                   <- "monthly date"
sjlabelled::set_label(seriesPriceSoybean.prepData$year)                   <- "calendar year"
sjlabelled::set_label(seriesPriceSoybean.prepData$month)                  <- "month indicator"
sjlabelled::set_label(seriesPriceSoybean.prepData$price_real_mon_soybean)  <- "real average monthly price, soybean (R$ per kg, constant 01/2017)"



# POST-TREATMENT OVERVIEW
# summary(seriesPriceSoybean.prepData)
# View(seriesPriceSoybean.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(seriesPriceSoybean.prepData,
     file = here::here("data/projectSpecific/prepData",
                      paste0("seriesPriceSoybean_prepData", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/prepData")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
