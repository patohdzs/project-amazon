library(sf)
library(terra)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter())
conflicts_prefer(terra::extract())

# Load clean mapbiomas raster
clean_mapbiomas <- rast("data/clean/land_use_cover_2000.tif")
n_species_rst <- rast("data/raw/TreeRichness_ha.asc")

# Summary stats
summary(n_species_rst)

# Load pixel data for primary forest
load("data/processed/pixel_primary_forest_2017.Rdata")

# Transform to sf
pixel_primary_forest_2017 <- st_as_sf(
  x = pixel_primary_forest_2017,
  coords = c("lon", "lat"),
  crs = st_crs(4326),
  remove = FALSE
)

# Transform to spatVector
aux_polygons <- vect(pixel_primary_forest_2017)

# Crop spatial points to raster extent
aux_polygons <- crop(aux_polygons, n_species_rst)

# Extract n_species raster data by spatial point
aux_n_species <- extract(n_species_rst, aux_polygons, xy = TRUE)

# Adjust and select column names
tree_richness_2017 <- aux_n_species %>%
  rename(lon = x, lat = y) %>%
  select(-ID) %>%
  st_as_sf(
    coords = c("lon", "lat"),
    crs = st_crs(4326),
    remove = FALSE
  ) %>%
  rename(tree_richness_ha = TreeRichness_ha)



# Save data set
save(tree_richness_2017, file = "data/processed/tree_richness_2017.Rdata")
