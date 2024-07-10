library(sf)
library(terra)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

in_file <- "data/raw/mapbiomas/land_use_cover/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-2017.tif"
land_use_rst <- rast(in_file)

# Read in secondary vegetation raster
sec_veg_rst <- rast("data/raw/mapbiomas/secondary_vegetation/brasil_desmat_vsec_anual_2017.tif")

# Crop the the Amazonia subset
sec_veg_rst <- sec_veg_rst |>
  crop(land_use_rst) |>
  mask(land_use_rst)

# Create dummy for secondary vegetation class
sec_veg_rst[sec_veg_rst != 303] <- 0
sec_veg_rst[sec_veg_rst == 303] <- 1

# Aggregate into pixel shares
sec_veg_rst <- aggregate(
  sec_veg_rst,
  fact = 2250,
  fun = sum,
  na.rm = TRUE
) / (2250^2)

# Rename raster layer
names(sec_veg_rst) <- "share_secondary_vegetation_2017"

# Write raster
out_file <- "data/processed/amazon_sec_veg_2017_shares_1043_sites.tif"
writeRaster(sec_veg_rst, out_file, overwrite = TRUE)

# Clean environment
rm(sec_veg_rst)

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
    forest_area_1995_ha = share_forest_1995 * site_area_ha,
    forest_area_2008_ha = share_forest_2008 * site_area_ha,
    forest_area_2017_ha = share_forest_2017 * site_area_ha,
    other_area_1995_ha = share_other_1995 * site_area_ha,
    other_area_2008_ha = share_other_2008 * site_area_ha,
    other_area_2017_ha = share_other_2017 * site_area_ha,
    sec_veg_area_2017_ha = share_secondary_vegetation_2017 * site_area_ha
  ) %>%
  mutate(
    z_1995 = share_agricultural_use_1995 * site_area_ha,
    z_2008 = share_agricultural_use_2008 * site_area_ha,
    z_2017 = share_agricultural_use_2017 * site_area_ha,
    zbar_1995 = forest_area_1995_ha + z_1995,
    zbar_2008 = forest_area_2008_ha + z_2008,
    zbar_2017 = forest_area_2017_ha + z_2017,
  ) %>%
  select(id, contains("area"), contains("z"))
