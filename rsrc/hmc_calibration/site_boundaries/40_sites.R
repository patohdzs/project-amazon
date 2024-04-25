# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (25 Sites MODEL)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(MASS)
library(dplyr)
library(boot)
library(conflicted)
library(tictoc)
library(sf)
library(terra)

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::summarize)

# START TIMER
tictoc::tic(msg = "calibration_40SitesModel.R script", log = TRUE)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))

# DATA INPUT
# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster_40_sites <-
  terra::rast(list.files(
    here::here("data/calibration/1043SitesModel/aux_tifs"),
    pattern = "raster_",
    full.names = TRUE
  ))


# MUNI LEVEL SPATIAL SAMPLE
load(here::here(
  "data/calibration/prepData/sampleMuniSpatial_prepData.Rdata"
))

# INITIAL CONDITIONS Z
# AGGREGATE FROM 1000 SITES TO 43 SITES
# transform shares to areas
raster_40_sites$amazonBiomeArea_ha_40Sites <-
  raster_40_sites$share_amazonBiome * raster_40_sites$pixelArea_ha
raster_40_sites$forestArea_1995_ha_40Sites <-
  raster_40_sites$share_forest_1995 * raster_40_sites$pixelArea_ha
raster_40_sites$agriculturalUseArea_1995_ha_40Sites <-
  raster_40_sites$share_agriculturalUse_1995 * raster_40_sites$pixelArea_ha
raster_40_sites$otherArea_1995_ha_40Sites <-
  raster_40_sites$share_other_1995 * raster_40_sites$pixelArea_ha
raster_40_sites$forestArea_2017_ha_40Sites <-
  raster_40_sites$share_forest_2017 * raster_40_sites$pixelArea_ha
raster_40_sites$agriculturalUseArea_2017_ha_40Sites <-
  raster_40_sites$share_agriculturalUse_2017 * raster_40_sites$pixelArea_ha
raster_40_sites$otherArea_2017_ha_40Sites <-
  raster_40_sites$share_other_2017 * raster_40_sites$pixelArea_ha

# select area variables
raster_40_sites <-
  terra::subset(
    raster_40_sites,
    c(
      "amazonBiomeArea_ha_40Sites",
      "pixelArea_ha",
      "forestArea_1995_ha_40Sites",
      "agriculturalUseArea_1995_ha_40Sites",
      "otherArea_1995_ha_40Sites",
      "forestArea_2017_ha_40Sites",
      "agriculturalUseArea_2017_ha_40Sites",
      "otherArea_2017_ha_40Sites"
    )
  )

# aggregate from 1000 sites to 43
raster_40_sites <-
  terra::aggregate(raster_40_sites,
    fact = 6,
    fun = sum,
    na.rm = TRUE
  )



# extract variables as polygons, transform to sf,
# and project data for faster spatial manipulation
calibration_40_sites_model <-
  terra::as.polygons(raster_40_sites, dissolve = FALSE) %>%
  sf::st_as_sf() %>%
  sf::st_transform(5880)

# add id variable
calibration_40_sites_model$id <- seq_len(nrow(calibration_40_sites_model))

# transform share variables in area (ha)
calibration_40_sites_model <-
  calibration_40_sites_model %>%
  dplyr::mutate(
    zbar_1995_40Sites = agriculturalUseArea_1995_ha_40Sites +
      forestArea_1995_ha_40Sites,
    zbar_2017_40Sites = agriculturalUseArea_2017_ha_40Sites +
      forestArea_2017_ha_40Sites
  ) %>%
  dplyr::select(
    amazonBiomeArea_ha_40Sites,
    siteArea_ha_40Sites = pixelArea_ha,
    forestArea_1995_ha_40Sites,
    z_1995_40Sites = agriculturalUseArea_1995_ha_40Sites,
    zbar_1995_40Sites,
    forestArea_2017_ha_40Sites,
    z_2017_40Sites = agriculturalUseArea_2017_ha_40Sites,
    zbar_2017_40Sites
  )


# remove sites with less than 2% of its are intersecting with the amazon biome
calibration_40_sites_model <-
  calibration_40_sites_model %>%
  dplyr::filter(amazonBiomeArea_ha_40Sites / siteArea_ha_40Sites >= 0.03)

# add id variable
calibration_40_sites_model$id <- seq_len(nrow(calibration_40_sites_model))


id <- calibration_40_sites_model %>%
  select(id)

st_write(id,
  "data/hmc/id_40.geojson",
  driver = "GeoJSON",
  delete_dsn = TRUE
)
