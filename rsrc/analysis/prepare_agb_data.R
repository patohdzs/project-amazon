library(sf)
library(tictoc)
library(terra)
library(tidyverse)
library(conflicted)
library(sjlabelled)

# Resolve conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

# Start timer
tic(msg = "merge_agb_sec_veg.R script", log = TRUE)

# Load calibrated gamma
load("data/calibration/gamma_calibration_1043_sites.Rdata")

# Load clean mapbiomas raster
clean_mapbiomas <- rast("data/clean/land_use_cover_2000.tif")

# Load secondary vegetation age data for 2017
sec_veg_age_rst <- rast(
  list.files(
    "data/raw/mapbiomas/secondary_vegetation_age/",
    pattern = "2017",
    full.names = TRUE
  )
)

# Load aboveground biomass raster data for 2017
agb_rst <- map(
  list.files(
    "data/raw/esa/above_ground_biomass/",
    pattern = "_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv3.0.tif",
    full.names = TRUE
  ),
  rast
)

pq_rst <- rast(
  list.files(
    "data/raw/mapbiomas/pasture_quality/",
    pattern = "2014",
    full.names = TRUE
  )
)

# Combine AGB zonal rasters into a single raster
agb_rst <- do.call(merge, agb_rst)

# Crop the the Amazonia subset
sec_veg_age_rst <- sec_veg_age_rst |>
  crop(clean_mapbiomas) |>
  mask(clean_mapbiomas)

# Aggregate secondary vegetation raster by 4
sec_veg_age_rst <- aggregate(
  sec_veg_age_rst,
  fact = 8,
  fun = median,
  na.rm = TRUE
)

# Resample to match sec_veg resolution and mask values outside Amazonia
agb_rst <- resample(agb_rst, sec_veg_age_rst, method = "near") |>
  mask(sec_veg_age_rst)

# Resample pasture quality to match sec_veg resolution
pq_rst <- resample(pq_rst, sec_veg_age_rst, method = "near") |>
  mask(sec_veg_age_rst)

gamma_rst <- calib_df |>
  rasterize(sec_veg_age_rst, field = "site_reg_gamma")

# Stack
carbon_accumulation_rst <- c(sec_veg_age_rst, agb_rst, pq_rst, gamma_rst)

# Convert to data frame
df <- carbon_accumulation_rst |>
  as.data.frame(xy = TRUE) |>
  as_tibble()

# Rename
names(df) <- c("x", "y", "age", "agb", "pq", "gamma")

# Convert AGB -> CO2e
df <- df |>
  mutate(co2e = (agb / 2) * (44 / 12))

# Save
save(df, file = "data/calibration/carbon_accumulation.Rdata")
