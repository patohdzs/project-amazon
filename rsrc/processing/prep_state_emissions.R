# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: PREPATE DATA TO ESTIMATE PARAMETER K (EMISSION FACTOR OF AGRICULTURAL SECTOR)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(tictoc)
library(tidyverse)
library(conflicted)
library(sjlabelled)

conflicts_prefer(dplyr::filter)

# Start timer
tic(msg = "state_emissions.R script", log = TRUE)

# Emissions data
load("data/clean/emission.Rdata")

# Land cover and use
load("data/clean/land_use_cover_muni.Rdata")

# Cross-sectional municipal-level sample
load("data/processed/cross_section_muni_sample.Rdata")

# Column and year selection + columns aggregation
land_use_cover_muni <-
  land_use_cover_muni %>%
  mutate(
    agriculturalUse_area = mapbiomasLandCoverId_15 +
      mapbiomasLandCoverId_39 +
      mapbiomasLandCoverId_41 +
      mapbiomasLandCoverId_20 +
      mapbiomasLandCoverId_21
  ) %>%
  select(muni_code, year, agriculturalUse_area)

# Merge data sets:
# - select period of interest
# - remove municipalities outside amazon biome
# - group and aggregate
state_emissions <-
  land_use_cover_muni %>%
  left_join(cross_section_muni_sample, by = c("muni_code")) %>%
  left_join(emission, by = c("state_uf", "year")) %>%
  filter(year >= 1990) %>%
  filter(!is.na(biomeAmazon_share)) %>%
  group_by(state_uf, year, emission_co2e, netEmission_co2e) %>%
  summarise(
    agriculturalUse_area = sum(agriculturalUse_area),
    biomeAmazon_share = sum(biomeAmazon_share * muni_area) / sum(muni_area)
  ) %>%
  mutate(
    emission_co2e = emission_co2e * biomeAmazon_share,
    netEmission_co2e = netEmission_co2e * biomeAmazon_share,
    agriculturalUse_area = agriculturalUse_area * biomeAmazon_share
  ) %>%
  select(state_uf, year, emission_co2e, netEmission_co2e, agriculturalUse_area)


# Clear environment
rm(emission, land_use_cover_muni, cross_section_muni_sample)

# Set labels
set_label(state_emissions$emission_co2e) <- "agricultural use emission factor adjusted by fraction of state area inside amazon biome (CO2e)"
set_label(state_emissions$netEmission_co2e) <- "agricultural use net emission factor adjusted by fraction of state area inside amazon biome (CO2e"
set_label(state_emissions$agriculturalUse_area) <- "agricultural area adjusted by fraction of state area inside amazon biome (ha)"

# Save data set
save(state_emissions, file = "data/processed/state_emissions.Rdata")

# End timer
toc(log = TRUE)
