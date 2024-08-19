library(sf)
library(glue)
library(tictoc)
library(terra)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter)

# Load raster data (amazon biome share, pixel areas, and land uses)
amazon_rsts <- rast(
  list.files(
    "data/processed/",
    pattern = ".tif",
    full.names = TRUE
  )
)

# Extract variables as polygons and project data
calib_df <- as.polygons(amazon_rsts, dissolve = FALSE) %>%
  st_as_sf() %>%
  st_transform(5880)

# Remove sites with less than 3% overlap with the amazon biome
calib_df <- calib_df %>%
  filter(share_amazon_biome >= 0.03)

# Add id variable
calib_df$id <- seq_len(nrow(calib_df))

# Transform share variables into area (ha)
calib_df <- calib_df %>%
  mutate(across(
    starts_with("share_"),
    .names = "area_{.col}",
    ~ . * pixel_area_ha
  )) %>%
  rename_with(
    ~ str_replace(., "share_", ""),
    starts_with("area_")
  )

# REMOVE
calib_df <- calib_df %>%
  rename_with(~ gsub("2008", "2017", .x), starts_with("area_pasture_quality_"))



# Compute total pasture area
calib_df <- calib_df %>%
  mutate(
    area_pasture_2017 = area_pasture_quality_1_2017 +
      area_pasture_quality_2_2017 +
      area_pasture_quality_3_2017
  )

# Compute pasture area as share of agricultural area (Z)
calib_df <- calib_df %>%
  mutate(
    share_pasture_z = area_pasture_2017 / area_agricultural_use_2017
  )

# Compute as share of total pasture area
calib_df <- calib_df %>%
  mutate(
    share_pq_1_pasture = area_pasture_quality_1_2017 / area_pasture_2017,
    share_pq_2_pasture = area_pasture_quality_2_2017 / area_pasture_2017,
    share_pq_3_pasture = area_pasture_quality_3_2017 / area_pasture_2017,
  )

# Compute as share of total agricultural area (Z)
calib_df <- calib_df %>%
  mutate(
    share_pq_1_z = area_pasture_quality_1_2017 / area_agricultural_use_2017,
    share_pq_2_z = area_pasture_quality_2_2017 / area_agricultural_use_2017,
    share_pq_3_z = area_pasture_quality_3_2017 / area_agricultural_use_2017,
  )

calib_df %>%
  pull(share_pasture_z) %>%
  summary()

calib_df %>% filter(share_pasture_z > 1) %>% nrow()

calib_df %>%
  select(
    share_pq_1_pasture,
    share_pq_2_pasture,
    share_pq_3_pasture
  ) %>%
  summary()

calib_df %>%
  select(
    share_pq_1_z,
    share_pq_2_z,
    share_pq_3_z
  ) %>%
  summary()
