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
tictoc::tic(msg = "emissionKuznets_raw2clean.R script", log = TRUE)

# Read input file
emission_in_path <- "data/raw/worldbank/emission_kuznets/API_EN.ATM.CO2E.PC_DS2_en_csv_v2_3731558.csv"
gdp_in_path <- "data/raw/worldbank/emission_kuznets/API_NY.GDP.PCAP.PP.CD_DS2_en_csv_v2_3731320.csv"

raw_emission_kuznets <- readr::read_csv(file = emission_in_path, skip = 4)
raw_gdp_kuznets <- readr::read_csv(file = gdp_in_path, skip = 4)

# DATASET CLEANUP AND PREP
raw_emission_kuznets <-
  raw_emission_kuznets %>%
  dplyr::select(`Country Name`, emissionPerCapita_2018 = `2018`) %>%
  dplyr::left_join(raw_gdp_kuznets) %>%
  dplyr::select(
    country_name = `Country Name`,
    gdpPerCapita_2018 = `2018`,
    emissionPerCapita_2018
  )


# LABELS
sjlabelled::set_label(raw_emission_kuznets$country_name) <- "name of the country"
sjlabelled::set_label(raw_emission_kuznets$gdpPerCapita_2018) <- "GDP per capita PPP in 2018 (current international $)"
sjlabelled::set_label(raw_emission_kuznets$emissionPerCapita_2018) <- "Emission per capita in 2018 (metric tons)"

# Change object name before saving
clean_emission_kuznets <- raw_emission_kuznets

# Save object
save(clean_emission_kuznets,
  file = "data/clean/emission_kuznets.Rdata"
)

# END TIMER
tictoc::toc(log = TRUE)
