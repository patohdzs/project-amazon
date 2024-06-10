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
tic(msg = "emission_raw2clean.R script", log = TRUE)

# Read Excel file
emission <- read_xlsx(
  path = "data/raw/seeg/emission/emission_agriculture_states.xlsx",
  sheet = 1, col_names = c("state_uf", 1990:2019), skip = 1
)

raw_removal <- read_xlsx(
  path = "data/raw/seeg/emission/removalNCI_agriculture_states.xlsx",
  sheet = 1, col_names = c("state_uf", 1990:2019), skip = 1
)

# RESHAPE
emission <-
  emission %>%
  pivot_longer(
    -state_uf,
    names_to = "year",
    values_to = "emission_co2e"
  )

raw_removal <-
  raw_removal %>%
  pivot_longer(
    -state_uf,
    names_to = "year",
    values_to = "removal_co2e"
  )

# MERGE
emission <-
  emission %>%
  left_join(raw_removal)

# ADD NET EMISSION VARIABLE
emission <-
  emission %>%
  mutate(netEmission_co2e = emission_co2e + removal_co2e) %>%
  mutate(year = as.numeric(year))

# Clean environmnet
rm(raw_removal)

# set_labelS
set_label(emission$state_uf) <- "state name abbreviation"
set_label(emission$year) <- "year of reference (calendar or PRODES year)"
set_label(emission$emission_co2e) <- "total emissions from agricultural land (CO2e-GWP-AR5)"
set_label(emission$removal_co2e) <- "total removals from agricultural land (CO2e-GWP-AR5)"
set_label(emission$netEmission_co2e) <- "total net emissions from agricultural land (CO2e-GWP-AR5)"

# Save data set
save(emission, file = "data/clean/emission.Rdata")

# END TIMER
toc(log = TRUE)
