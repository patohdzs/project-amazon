
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



# START TIMER
tictoc::tic(msg = "seriesPriceCattle_prepData.R script", log = TRUE)



# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# AGRICUTLURAL COMMODITY PRICES
load("data/clean/commodity_prices.Rdata")

# DEFLATOR (IPA-EP-DI)
load("data/clean/deflatorIPA.Rdata")


# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

aux_01_2017 <- clean_deflatorIPA[clean_deflatorIPA$date == "2017-01-01",]$deflator_ipa

# TRANSFORM NOMINAL MONTHLY COMMODITY PRICES TO REAL STANDARDIZED PRICES
seriesPriceCattle_prepData <-
  clean_commodityPrices %>%
  dplyr::left_join(clean_deflatorIPA) %>% # merge commodity nominal prices with deflator index
  dplyr::mutate(deflator_ipa = deflator_ipa/aux_01_2017) %>% # change base year to january 2017
  dplyr::mutate(price_real_mon_cattle = price_cattle/deflator_ipa) %>%  # transform to real price (constant january 2017)
  dplyr::mutate(month = lubridate::month(date),
                year = lubridate::year(date)) %>% # construct year, month and trimester variables
  dplyr::select(date, year, month, starts_with("price_real"), price_cattle)


# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(seriesPriceCattle_prepData$date)                   <- "monthly date"
sjlabelled::set_label(seriesPriceCattle_prepData$year)                   <- "calendar year"
sjlabelled::set_label(seriesPriceCattle_prepData$month)                  <- "month indicator"
sjlabelled::set_label(seriesPriceCattle_prepData$price_real_mon_cattle)  <- "real average monthly price, cattle (R$ per @, constant 01/2017)"



# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(seriesPriceCattle_prepData,
     file = "data/prepData/seriesPriceCattle_prepData.Rdata")

readr::write_csv(seriesPriceCattle_prepData,
  file = "data/calibration/hmc/seriesPriceCattle_prepared.csv"
)

# END TIMER
tictoc::toc(log = TRUE)


# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
