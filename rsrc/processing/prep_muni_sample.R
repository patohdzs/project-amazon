# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: DEFINE MUNI-LEVEL SAMPLE AND SAVE IT IN SPATIAL, CROSS-SECTION AND PANEL FORMATS
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -
library(sf)
library(units)
library(tictoc)
library(tidyverse)
library(sjlabelled)

# Start timer
tic(msg = "prep_muni_sample.R script", log = TRUE)

# BRAZILIAN MUNICIPALITIES DIVISION 2015 SHAPEFILE
load("data/clean/muni_division_2015.Rdata")

# AMAZON BIOME BOUNDARY
load("data/clean/amazon_biome.Rdata")

# Calculate municipality area in square kilometers
muni_division_2015$muni_area <-
  st_area(muni_division_2015) %>%
  set_units(ha) %>%
  unclass()

# Select columns of interest
muni_division_2015 <-
  muni_division_2015 %>%
  select(
    muni_code,
    muni_name,
    state_uf,
    muni_area,
    geometry
  )

amazon_biome <-
  amazon_biome %>%
  select(biome_name, geometry)

# Combine municipalities with biomes
aux_biome_muni <- st_intersection(muni_division_2015, amazon_biome)

# Clear environment
rm(amazon_biome)

# Calculate biome areas inside each municipality
aux_biome_muni$biome_area <-
  st_area(aux_biome_muni) %>%
  set_units(ha) %>%
  unclass()

# Select municipalities of interest - in Amazon or Cerrado
# Convert from long to wide - biome areas in columns
# Select municipalities of interest
aux_biome_muni <-
  aux_biome_muni %>%
  st_drop_geometry() %>%
  pivot_wider(
    names_from = biome_name,
    values_from = biome_area,
    values_fill = 0
  ) %>%
  filter(AMAZON > 0) %>%
  mutate(biomeAmazon_share = AMAZON / muni_area) %>%
  select(-AMAZON)

# Merge original municipalities shapefile with biomes information
spatial_muni_sample <-
  muni_division_2015 %>%
  left_join(aux_biome_muni) %>%
  filter(biomeAmazon_share > 0)

# Clear environment
rm(aux_biome_muni)

# Set labels
set_label(spatial_muni_sample$muni_code) <- "municipality code (7-digit, IBGE - 2015)"
set_label(spatial_muni_sample$muni_name) <- "municipality name"
set_label(spatial_muni_sample$state_uf) <- "state name (abbreviation)"
set_label(spatial_muni_sample$muni_area) <- "municipality area (ha, calculated from shapefile under SIRGAS2000 Polyconic projection)"
set_label(spatial_muni_sample$biomeAmazon_share) <- "share of the municipality area in the Amazon biome"


# Extract data frame from spatial data
cross_section_muni_sample <-
  spatial_muni_sample %>%
  st_drop_geometry()

# Create panel data from cross section and add missing label
panel_muni_sample <- expand_grid(cross_section_muni_sample, year = 2000:2019)
set_label(panel_muni_sample$year) <- "year of reference (calendar or PRODES year)"

# Save data sets
save(spatial_muni_sample,
  file = "data/processed/spatial_muni_sample.Rdata"
)

save(cross_section_muni_sample,
  file = "data/processed/cross_section_muni_sample.Rdata"
)

save(panel_muni_sample,
  file = "data/processed/panel_muni_sample.Rdata"
)

# End timer
toc(log = TRUE)
