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
biome_raster <- rast("data/raw/mapbiomas/land_use_cover/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-2017.tif")

# Load Amazon biome data
load("data/clean/amazon_biome.Rdata")

# Change biome crs to match raster
amazon_biome <- st_transform(amazon_biome, crs(biome_raster))

# Rasterize amazon biome into 30m raster resolution to minimize area distortion
biome_raster <- rasterize(vect(amazon_biome), biome_raster, field = 1)

# Clean environment
rm(amazon_biome)

# Aggregate raster calculating the share of minicells that are in the biome
# 2250^2 is the total number of minicells
biome_raster <- aggregate(
  biome_raster,
  fact = 2250,
  fun = sum,
  na.rm = TRUE
) / (2250^2)

# Add layer name
names(biome_raster) <- "share_amazon_biome"

# Save biome shares raster
out_file <- "data/processed/amazon_biome_shares_1043_sites.tif"
writeRaster(biome_raster, out_file, overwrite = TRUE)

# Convert cell size to hectares
biome_raster_ha <- cellSize(biome_raster, unit = "ha")

# Clean environment
rm(biome_raster)

# Add layer name
names(biome_raster_ha) <- "pixel_area_ha"

# Save areas raster
out_file <- "data/processed/amazon_biome_areas_1043_sites.tif"
writeRaster(biome_raster_ha, out_file, overwrite = TRUE)

# Clean environment
rm(biome_raster_ha)

# End timer
toc(log = TRUE)
