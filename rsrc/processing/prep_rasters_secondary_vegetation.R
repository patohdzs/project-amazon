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
