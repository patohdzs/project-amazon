# Load necessary libraries
library(sf)
library(stars)
library(tictoc)
library(terra)
library(tidyverse)
library(conflicted)
library(sjlabelled)

# Resolve conflicts
conflicts_prefer(dplyr::filter())
conflicts_prefer(terra::extract())

# Start timer
tic(msg = "merge_agb_sec_veg.R script", log = TRUE)

# Load secondary vegetation age data for 2017
sec_veg_age_rst <- rast(
  list.files(
    "data/raw/mapbiomas/secondary_vegetation_age/",
    pattern = "2017",
    full.names = TRUE
  )
)

# Load aboveground biomass raster data
agb_rasters  <- map(
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


# Aggregate secondary vegetation raster by 4

sec_veg_age_aggregated <- aggregate(sec_veg_age_rst, fact = 4, fun = mean, na.rm = TRUE)


writeRaster(sec_veg_age_aggregated, "data/calibration/sec_veg_age_aggregated.tif", overwrite = TRUE)


# Combine resampled AGB rasters into a single raster stack
combined_agb_raster <- do.call(merge, agb_rasters)

# resample to match sec_veg resolution 
resampled_agb_raster <- resample(combined_agb_raster, sec_veg_age_aggregated, method = "near")


writeRaster(resampled_agb_raster, "data/calibration/resampled_agb_aggregated.tif", overwrite = TRUE)


# Resample pasture quality to match sec_veg resolution
resampled_pq_rst <- resample(pq_rst, sec_veg_age_aggregated, method = "near")

writeRaster(resampled_pq_rst, "data/calibration/resampled_pq_rst.tif", overwrite = TRUE)

# Resample gamma with secondary vege 
load("data/calibration/gamma_calibration_1043_sites.Rdata")

raster_gamma<-as(st_rasterize(calib_df %>% dplyr::select(site_reg_gamma, geometry)),"SpatRaster")

resampled_gamma <- resample(raster_gamma, sec_veg_age_aggregated, method = "near")

writeRaster(resampled_gamma, "data/calibration/resampled_gamma.tif", overwrite = TRUE)




#### Transform into dataframe


resampled_agb_raster <- rast("data/calibration/resampled_agb_aggregated.tif")
sec_veg_age_aggregated <- rast("data/calibration/sec_veg_age_aggregated.tif")
resampled_gamma <- rast("data/calibration/resampled_gamma.tif")
resampled_pq_rst <- rast("data/calibration/resampled_pq_rst.tif")

aggregate_agb_df <- data.frame(terra::values(resampled_agb_raster))

aggregate_sec_df <- data.frame(terra::values(sec_veg_age_aggregated))

aggregate_gamma_df <- data.frame(terra::values(resampled_gamma))

aggregate_resampled_pq_rst <- data.frame(terra::values(resampled_pq_rst))



combined_df <- data.frame(
  agb = ((aggregate_agb_df$N00W050_ESACCI.BIOMASS.L4.AGB.MERGED.100m.2017.fv3.0 / 2) * (44 / 12)),
  sec = aggregate_sec_df$sec_veg_age_2017,
  gamma = aggregate_gamma_df$site_reg_gamma,
  pq = aggregate_resampled_pq_rst$pasture_quality_2014
) %>%
  filter(agb != 0, sec != 0)

save(combined_df, file = "data/calibration/combined_df.Rdata")






