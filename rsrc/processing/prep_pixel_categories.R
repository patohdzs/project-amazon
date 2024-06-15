# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CREATE AGGREGATED CATEGORIES OF INTEREST BASED ON MAPBIOMAS 30M-PIXELS VALUES
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(sf)
library(tictoc)
library(terra)
library(glue)
library(tidyverse)
library(conflicted)
library(sjlabelled)

conflicts_prefer(dplyr::filter())
conflicts_prefer(terra::extract())


# START TIMER
tic(msg = "pixel_categories.R script", log = TRUE)

# Mapbiomas sample
load("data/processed/pixel_areas.Rdata")

# Add aggregated land use/cover categories
# Note: perennial crops not present in amazon, mosaic only small area in 1985
pixel_categories <-
  pixel_areas %>%
  mutate(mapbiomas_classAgg = case_when(
    mapbiomas_class == 3 ~ "forest",
    mapbiomas_class == 15 ~ "pasture",
    mapbiomas_class == 39 ~ "soybean",
    mapbiomas_class %in% c(20, 41, 21) ~ "otherCrops",
    mapbiomas_class %in% c(4, 5, 9, 11:13, 22:27, 29:33) ~ "otherCategories"
  ))

# Clear environment
rm(pixel_areas)

# Save pixels with primary forest by year
map(
  .x = c(2018),
  .f = function(.x) {
    # Filter primary forest pixels
    pixel_primary_forest <-
      pixel_categories %>%
      filter(year <= .x) %>%
      mutate(d_forest = if_else(mapbiomas_classAgg == "forest", 1, 0)) %>%
      group_by(lon, lat) %>%
      mutate(d_primaryForest = if_else(sum(d_forest) == (.x - 1985 + 1), 1, 0)) %>%
      ungroup() %>%
      filter(year == .x, d_primaryForest == 1)

    # Save data set
    out_file <- glue("data/processed/pixel_primary_forest_{.x}.Rdata")
    save(pixel_primary_forest, file = out_file)
  }
)

# Set labels
set_label(pixel_categories$mapbiomas_classAgg) <- "mapbiomas land use/cover aggregated classification"

# Save pixel categories data set
save(pixel_categories, file = "data/processed/pixel_categories.Rdata")

# END TIMER
toc(log = TRUE)
