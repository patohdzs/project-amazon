# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CLEAN RAW SEEG AGRICULTURAL EMISSION DATA - LEGAL AMAZON STATES 1990-2019
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# -
library(tidyverse)
library(tictoc)
library(sjlabelled)
library(conflicted)
library(readxl)

# Resolve conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

# START TIMER
tic(msg = "emissions_raw2clean.R script", log = TRUE)

# Read Excel file
emissions <- read_xlsx(
  path = "data/raw/seeg/emission/emission_agriculture_states.xlsx",
  sheet = 1, col_names = c("state_uf", 1990:2019), skip = 1
)

raw_removal <- read_xlsx(
  path = "data/raw/seeg/emission/removalNCI_agriculture_states.xlsx",
  sheet = 1, col_names = c("state_uf", 1990:2019), skip = 1
)

# RESHAPE
emissions <-
  emissions %>%
  pivot_longer(
    -state_uf,
    names_to = "year",
    values_to = "emissions_co2e"
  )

raw_removal <-
  raw_removal %>%
  pivot_longer(
    -state_uf,
    names_to = "year",
    values_to = "removal_co2e"
  )

# MERGE
emissions <-
  emissions %>%
  left_join(raw_removal)

# Add net emissions variable
emissions <-
  emissions %>%
  mutate(net_emissions_co2e = emissions_co2e + removal_co2e) %>%
  mutate(year = as.numeric(year))

# Clean environmnet
rm(raw_removal)

# set_labelS
set_label(emissions$state_uf) <- "state name abbreviation"
set_label(emissions$year) <- "year of reference (calendar or PRODES year)"
set_label(emissions$emissions_co2e) <- "total emissions from agricultural land (CO2e-GWP-AR5)"
set_label(emissions$removal_co2e) <- "total removals from agricultural land (CO2e-GWP-AR5)"
set_label(emissions$net_emissions_co2e) <- "total net emissions from agricultural land (CO2e-GWP-AR5)"

# Save data set
save(emissions, file = "data/clean/emissions.Rdata")

# END TIMER
toc(log = TRUE)
