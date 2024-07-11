library(sf)
library(glue)
library(tictoc)
library(terra)
library(conflicted)

conflicts_prefer(terra::extract())

# START TIMER
tic(msg = "prep_biome_class_rasters.R script", log = TRUE)

in_file <- "data/raw/mapbiomas/pasture_quality/pasture-quality-amazonia-2017.tif"

for (level in c(1, 2, 3)) {
  # Read in pasture quality raster
  pasture_quality_rst <- rast(in_file)

  # Create dummy for pasture quality class
  pasture_quality_rst[pasture_quality_rst != level] <- 0
  pasture_quality_rst[pasture_quality_rst == level] <- 1

  # Aggregate into pixel shares
  pasture_quality_rst <- aggregate(
    pasture_quality_rst,
    fact = 675,
    fun = sum,
    na.rm = TRUE
  ) / (675^2)

  # Rename raster layer
  names(pasture_quality_rst) <- glue("share_pasture_quality_{level}_2017")

  # Write raster
  out_file <- glue("data/processed/amazon_pq_{level}_2017_shares_1043_sites.tif")
  writeRaster(pasture_quality_rst, out_file, overwrite = TRUE)
}
