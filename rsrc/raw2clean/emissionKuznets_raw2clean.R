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


# START TIMER
tictoc::tic(msg = "emissionKuznets_raw2clean.R script", log = TRUE)

# DATA INPUT

# read input file
raw_emissionKuznets <- readr::read_csv(file = "data/raw/worldbank/emission_kuznets/API_EN.ATM.CO2E.PC_DS2_en_csv_v2_3731558.csv", skip = 4)
raw_gdpKuznets <- readr::read_csv(file = "data/raw/worldbank/emission_kuznets/API_NY.GDP.PCAP.PP.CD_DS2_en_csv_v2_3731320.csv", skip = 4)


# DATASET CLEANUP AND PREP

raw_emissionKuznets <-
  raw_emissionKuznets %>%
  dplyr::select(`Country Name`, emissionPerCapita_2018 = `2018`) %>%
  dplyr::left_join(raw_gdpKuznets) %>%
  dplyr::select(country_name = `Country Name`, gdpPerCapita_2018 = `2018`, emissionPerCapita_2018)

# EXPORT PREP

# LABELS
sjlabelled::set_label(raw_emissionKuznets$country_name) <- "name of the country"
sjlabelled::set_label(raw_emissionKuznets$gdpPerCapita_2018) <- "GDP per capita PPP in 2018 (current international $)"
sjlabelled::set_label(raw_emissionKuznets$emissionPerCapita_2018) <- "Emission per capita in 2018 (metric tons)"

# change object name for exportation
clean_emissionKuznets <- raw_emissionKuznets

# EXPORT

save(clean_emissionKuznets,
  file = 
    "data/clean/emission_kuznets.Rdata"
)

# END TIMER
tictoc::toc(log = TRUE)


