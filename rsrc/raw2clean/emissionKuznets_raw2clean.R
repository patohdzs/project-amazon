
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: TREAT RAW DATA EMISSIONS AND GDP BY COUNTRY (WORLD BANK)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("code/setup.R")


# START TIMER
tictoc::tic(msg = "emissionKuznets_raw2clean.R script", log = T)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read input file
raw.emissionKuznets <- readr::read_csv(file = here::here("data/raw2clean/emissionKuznets_worldbank/input/API_EN.ATM.CO2E.PC_DS2_en_csv_v2_3731558.csv"), skip = 4)
raw.gdpKuznets <- readr::read_csv(file = here::here("data/raw2clean/emissionKuznets_worldbank/input/API_NY.GDP.PCAP.PP.CD_DS2_en_csv_v2_3731320.csv"), skip = 4)


# DATA EXPLORATION [disabled for speed]
# summary(raw.emissionKuznets)
# View(raw.emissionKuznets)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

raw.emissionKuznets <-
  raw.emissionKuznets %>%
  dplyr::select(`Country Name`, emissionPerCapita_2018 = `2018`) %>%
  dplyr::left_join(raw.gdpKuznets) %>%
  dplyr::select(country_name = `Country Name`, gdpPerCapita_2018 = `2018`, emissionPerCapita_2018)




# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(raw.emissionKuznets$country_name)  <- "name of the country"
sjlabelled::set_label(raw.emissionKuznets$gdpPerCapita_2018)  <- "GDP per capita PPP in 2018 (current international $)"
sjlabelled::set_label(raw.emissionKuznets$emissionPerCapita_2018)  <- "Emission per capita in 2018 (metric tons)"




# change object name for exportation
clean.emissionKuznets <- raw.emissionKuznets



# POST-TREATMENT OVERVIEW
# summary(clean.emissionKuznets)
# View(clean.emissionKuznets)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.emissionKuznets,
     file = here::here("data/raw2clean/emissionKuznets_worldbank/output",
                      "clean_emissionKuznets.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("code/raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------