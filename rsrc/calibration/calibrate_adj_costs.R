library(sf)
library(terra)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

# Load raster data (amazon biome share, pixel areas, and land uses)
raster_variables <- rast(
  list.files(
    "data/processed/",
    pattern = "amazon_",
    full.names = TRUE
  )
)

# Extract variables as polygons and project data
calib_df <- as.polygons(raster_variables, dissolve = FALSE) %>%
  st_as_sf() %>%
  st_transform(5880)

# Remove sites with less than 3% overlap with the amazon biome
calib_df <- calib_df %>%
  filter(share_amazon_biome >= 0.03)

# Add id variable
calib_df$id <- seq_len(nrow(calib_df))

# Transform share variables into area (ha)
calib_df <- calib_df %>%
  rename(site_area_ha = pixel_area_ha) %>%
  mutate(
    amazon_biome_area_ha = share_amazon_biome * site_area_ha,
    forest_area_2017_ha = share_forest_2017 * site_area_ha,
    other_area_2017_ha = share_other_2017 * site_area_ha,
    sec_veg_area_2017_ha = share_secondary_vegetation_2017 * site_area_ha,
    pasture_qual_1_area_2017_ha = share_pasture_quality_1_2017 * site_area_ha,
    pasture_qual_2_area_2017_ha = share_pasture_quality_2_2017 * site_area_ha,
    pasture_qual_3_area_2017_ha = share_pasture_quality_3_2017 * site_area_ha
  ) %>%
  mutate(
    z_2017 = share_agricultural_use_2017 * site_area_ha,
    zbar_2017 = forest_area_2017_ha + z_2017
  ) %>%
  rowwise() %>%
  mutate(sum_past = sum(c_across(starts_with("pasture_qual_"))))



# Check if areas with poor pasture have less secondary vegetation
model_1 <- lm(
  log(1 + sec_veg_area_2017_ha) ~ log(1 + pasture_qual_1_area_2017_ha) +
    log(1 + z_2017),
  data = calib_df
)

model_2 <- lm(
  sec_veg_area_2017_ha ~ pasture_qual_1_area_2017_ha + z_2017,
  data = calib_df
)

model_2 <- lm(
  sec_veg_area_2017_ha ~ pasture_qual_1_area_2017_ha + z_2017,
  data = calib_df
)
