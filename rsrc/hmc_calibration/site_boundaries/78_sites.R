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
tictoc::tic(msg = "calibration_78SitesModel.R script", log = TRUE)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))


# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster_78_sites <-
  terra::rast(list.files(
    here::here("data/calibration/1055SitesModel/aux_tifs"),
    pattern = "raster_",
    full.names = TRUE
  ))

# MUNI LEVEL SPATIAL SAMPLE
load(here::here(
  "data/calibration/prepData/sampleMuniSpatial_prepData.Rdata"
))


# AGGREGATE FROM 1000 SITES TO 43 SITES
# transform shares to areas
raster_78_sites$amazonBiomeArea_ha_78Sites <-
  raster_78_sites$share_amazonBiome * raster_78_sites$pixelArea_ha
raster_78_sites$forestArea_1995_ha_78Sites <-
  raster_78_sites$share_forest_1995 * raster_78_sites$pixelArea_ha
raster_78_sites$agriculturalUseArea_1995_ha_78Sites <-
  raster_78_sites$share_agriculturalUse_1995 * raster_78_sites$pixelArea_ha
raster_78_sites$otherArea_1995_ha_78Sites <-
  raster_78_sites$share_other_1995 * raster_78_sites$pixelArea_ha
raster_78_sites$forestArea_2017_ha_78Sites <-
  raster_78_sites$share_forest_2017 * raster_78_sites$pixelArea_ha
raster_78_sites$agriculturalUseArea_2017_ha_78Sites <-
  raster_78_sites$share_agriculturalUse_2017 * raster_78_sites$pixelArea_ha
raster_78_sites$otherArea_2017_ha_78Sites <-
  raster_78_sites$share_other_2017 * raster_78_sites$pixelArea_ha

# select area variables
raster_78_sites <-
  terra::subset(
    raster_78_sites,
    c(
      "amazonBiomeArea_ha_78Sites",
      "pixelArea_ha",
      "forestArea_1995_ha_78Sites",
      "agriculturalUseArea_1995_ha_78Sites",
      "otherArea_1995_ha_78Sites",
      "forestArea_2017_ha_78Sites",
      "agriculturalUseArea_2017_ha_78Sites",
      "otherArea_2017_ha_78Sites"
    )
  )

# aggregate from 1000 sites to 78
raster_78_sites <-
  terra::aggregate(raster_78_sites,
    fact = 4,
    fun = sum,
    na.rm = TRUE
  )



# extract variables as polygons, transform to sf,
# and project data for faster spatial manipulation
calibration_78_sites_model <-
  terra::as.polygons(raster_78_sites, dissolve = FALSE) %>%
  sf::st_as_sf() %>%
  sf::st_transform(5880)




# add id variable
calibration_78_sites_model$id <- seq_len(nrow(calibration_78_sites_model))

# transform share variables in area (ha)
calibration_78_sites_model <-
  calibration_78_sites_model %>%
  dplyr::mutate(
    zbar_1995_78Sites = agriculturalUseArea_1995_ha_78Sites +
      forestArea_1995_ha_78Sites,
    zbar_2017_78Sites = agriculturalUseArea_2017_ha_78Sites +
      forestArea_2017_ha_78Sites
  ) %>%
  dplyr::select(
    amazonBiomeArea_ha_78Sites,
    siteArea_ha_78Sites = pixelArea_ha,
    forestArea_1995_ha_78Sites,
    z_1995_78Sites = agriculturalUseArea_1995_ha_78Sites,
    zbar_1995_78Sites,
    forestArea_2017_ha_78Sites,
    z_2017_78Sites = agriculturalUseArea_2017_ha_78Sites,
    zbar_2017_78Sites
  )


# remove sites with less than 2% of its are intersecting with the amazon biome
calibration_78_sites_model <-
  calibration_78_sites_model %>%
  dplyr::filter(amazonBiomeArea_ha_78Sites / siteArea_ha_78Sites >= 0.03)

# add id variable
calibration_78_sites_model$id <- seq_len(nrow(calibration_78_sites_model))



id <- calibration_78_sites_model %>%
  select(id)

st_write(id,
  "data/hmc/id_78.geojson",
  driver = "GeoJSON",
  delete_dsn = TRUE
)
