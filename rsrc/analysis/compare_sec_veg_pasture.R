library(sf)
library(terra)
library(stargazer)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

# Load raster data
clean_mapbiomas <- rast("data/clean/land_use_cover_2000.tif")

pq_rst <- rast(
  list.files(
    "data/raw/mapbiomas/pasture_quality/",
    pattern = "2014",
    full.names = TRUE
  )
)

sec_veg_age_rst <- rast(
  list.files(
    "data/raw/mapbiomas/secondary_vegetation_age/",
    pattern = "2016",
    full.names = TRUE
  )
)

# Crop the the Amazonia subset
sec_veg_age_rst <- sec_veg_age_rst |>
  crop(clean_mapbiomas) |>
  mask(clean_mapbiomas)

# Resample pasture quality to match sec_veg resolution
pq_rst <- resample(pq_rst, sec_veg_age_rst, method = "near")

# How much two-year-old sec veg has no pasture?
sec_veg_no_pasture <- (sec_veg_age_rst == 2) & (pq_rst == 0)
names(sec_veg_no_pasture) <- "sec_veg_no_pasture"

sec_veg_two_yo <- (sec_veg_age_rst == 2)
names(sec_veg_two_yo) <- "sec_veg_two_yo"

pixel_areas <- cellSize(sec_veg_age_rst, unit = "ha")

areas_1 <- zonal(pixel_areas, sec_veg_no_pasture, sum, na.rm = TRUE)
print(areas_1)

areas_2 <- zonal(pixel_areas, sec_veg_two_yo, sum, na.rm = TRUE)
print(areas_2)

print(areas_1[2, 2] / areas_2[2, 2])
