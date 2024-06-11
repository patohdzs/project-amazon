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
library(tidyverse)
library(tictoc)
library(sjlabelled)
library(conflicted)

# Resolve conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

# START TIMER
tic(msg = "emissionKuznets_raw2clean.R script", log = TRUE)

# Read input file
emission_in_path <- "data/raw/worldbank/emission_kuznets/API_EN.ATM.CO2E.PC_DS2_en_csv_v2_3731558.csv"
gdp_in_path <- "data/raw/worldbank/emission_kuznets/API_NY.GDP.PCAP.PP.CD_DS2_en_csv_v2_3731320.csv"

emission_kuznets <- read_csv(file = emission_in_path, skip = 4)
raw_gdp_kuznets <- read_csv(file = gdp_in_path, skip = 4)

# DATASET CLEANUP AND PREP
emission_kuznets <-
  emission_kuznets %>%
  select(`Country Name`, emissionPerCapita_2018 = `2018`) %>%
  left_join(raw_gdp_kuznets) %>%
  select(
    country_name = `Country Name`,
    gdpPerCapita_2018 = `2018`,
    emissionPerCapita_2018
  )


# LABELS
set_label(emission_kuznets$country_name) <- "name of the country"
set_label(emission_kuznets$gdpPerCapita_2018) <- "GDP per capita PPP in 2018 (current international $)"
set_label(emission_kuznets$emissionPerCapita_2018) <- "Emission per capita in 2018 (metric tons)"

# Save data set
save(emission_kuznets, file = "data/clean/emission_kuznets.Rdata")

# END TIMER
toc(log = TRUE)
