# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: TREAT RAW DATA - AMAZON BIOME BONUDARY (IBGE - 2019)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(sf)
library(tidyverse)
library(tictoc)
library(sjlabelled)
library(conflicted)

# Resolve conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

# START TIMER
tic(msg = "amazonBiome_raw2clean.R script", log = TRUE)

# Read shapefile
raw_biome <- st_read(
  dsn = "data/raw/ibge/amazon_biome/",
  layer = "lm_bioma_250"
)


# Column names
colnames(raw_biome)

# Translate column names
raw_biome <-
  raw_biome %>%
  rename(
    biome_code = CD_Bioma,
    biome_name = Bioma
  )

# Class - no change needed
lapply(raw_biome, class)

# TRANSLATION
# 'grepl' used to avoid encoding trouble with latin characters
raw_biome$biome_name[which(grepl(pattern = "Amazônia", x = raw_biome$biome_name))] <- "Amazon"
raw_biome$biome_name[which(grepl(pattern = "Mata Atlântica", x = raw_biome$biome_name))] <- "Atlantic Forest"

# LETTERS CAPITALIZATION
raw_biome <-
  raw_biome %>%
  mutate(biome_name = toupper(biome_name))

# FILTER BIOME OF INTEREST (AMAZON)
raw_biome <-
  raw_biome %>%
  filter(biome_name == "AMAZON")


# Project to CRS 4326 and save
output_path <- "data/calibration/hmc/map.geojson"
raw_biome %>%
  st_transform(crs = 4326) %>%
  st_write(output_path, driver = "GeoJSON")

# PROJECTION
# SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)
raw_biome <- st_transform(x = raw_biome, crs = 5880)

# GEOMETRY CLEANUP
raw_biome <- st_make_valid(raw_biome)

# LABELS
set_label(raw_biome$biome_code) <- "biome code"
set_label(raw_biome$biome_name) <- "biome name"

# Change object name before saving
amazon_biome <- raw_biome

# Save data set
save(amazon_biome, file = "data/clean/amazon_biome.Rdata")

# END TIMER
toc(log = TRUE)
