# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: GENERATE AGGREGATED SAMPLE OF INTEREST (DIVIDE AMAZON BIOME INTO 1055 CELLS)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(sf)
library(tictoc)
library(terra)
library(conflicted)

conflicts_prefer(terra::extract())

# Start timer
tic(msg = "prep_biome_grid.R script", log = TRUE)

# Load MapBiomas land use data (30m minicells)
raster_biome <- rast("data/raw/mapbiomas/land_use_cover/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-2017.tif")

# Load amazon biome data
load("data/clean/amazon_biome.Rdata")

# Change biome crs to match raster
amazon_biome <- st_transform(amazon_biome, crs(raster_biome))

# Rasterize amazon biome into 30m raster resolution to minimize area distortion
raster_biome <- rasterize(vect(amazon_biome), raster_biome, field = 1)

# Clean environment
rm(amazon_biome)

# Aggregate raster calculating the share of minicells that are in the biome
# 2250^2 is the total number of minicells
raster_biome <- aggregate(
  raster_biome,
  fact = 2250,
  fun = sum,
  na.rm = TRUE
) / (2250^2)

# Add name
names(raster_biome) <- "share_amazon_biome"

# Save raster
out_file <- "data/processed/amazon_biome_1043_sites.tif"
writeRaster(raster_biome, out_file, overwrite = TRUE)

# Convert cell size to hectares
raster_biome_ha <- terra::cellSize(raster_biome, unit = "ha")

# Clean environment
rm(raster_biome)

# Add layer name
names(raster_biome_ha) <- "pixel_area_ha"

# Save raster
out_file <- "data/processed/amazon_biome_areas_1043_sites.tif"
terra::writeRaster(raster_biome_ha, , overwrite = TRUE)

# Clean environment
rm(raster_biome_ha)

# End timer
toc(log = TRUE)
