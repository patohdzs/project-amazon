
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CONSTRUCT TRIMESTRAL COMMODITY REAL PRICES INDICES
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("code/setup.R")


# START TIMER
tictoc::tic(msg = "seriesPriceCattle_prepData.R script", log = T)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# AGRICUTLURAL COMMODITY PRICES
load(here::here("data/raw2clean/commodityPrices_seabpr/output/clean_commodityPrices.Rdata"))



# DEFLATOR (IPA-EP-DI)
load(here::here("data/raw2clean/deflatorIPA_fgv/output/clean_deflatorIPA.Rdata"))





# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

aux.01.2017 <- clean.deflatorIPA[clean.deflatorIPA$date == "2017-01-01",]$deflator_ipa

# TRANSFORM NOMINAL MONTHLY COMMODITY PRICES TO REAL STANDARDIZED PRICES
seriesPriceCattle.prepData <-
  clean.commodityPrices %>%
  dplyr::left_join(clean.deflatorIPA) %>% # merge commodity nominal prices with deflator index
  dplyr::mutate(deflator_ipa = deflator_ipa/aux.01.2017) %>% # change base year to january 2017
  dplyr::mutate(price_real_mon_cattle = price_cattle/deflator_ipa) %>%  # transform to real price (constant january 2017)
  dplyr::mutate(month = lubridate::month(date),
                year = lubridate::year(date)) %>% # construct year, month and trimester variables
  dplyr::select(date, year, month, starts_with("price_real"), price_cattle)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(seriesPriceCattle.prepData$date)                   <- "monthly date"
sjlabelled::set_label(seriesPriceCattle.prepData$year)                   <- "calendar year"
sjlabelled::set_label(seriesPriceCattle.prepData$month)                  <- "month indicator"
sjlabelled::set_label(seriesPriceCattle.prepData$price_real_mon_cattle)  <- "real average monthly price, cattle (R$ per @, constant 01/2017)"



# POST-TREATMENT OVERVIEW
# summary(seriesPriceCattle.prepData)
# View(seriesPriceCattle.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(seriesPriceCattle.prepData,
     file = here::here("data/calibration/prepData",
                      paste0("seriesPriceCattle_prepData", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("code/calibration")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
