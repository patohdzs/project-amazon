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



# START TIMER
tictoc::tic(msg = "10_sites.R script", log = TRUE)

# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tmpdir = "data/_temp",
                      timer  = T)

# DATA INPUT
# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster_10_sites <- terra::rast(list.files("data/calibration/1043SitesModel/aux_tifs",
                                           pattern = "raster_",
                                           full.names = T))


# MUNI LEVEL SPATIAL SAMPLE
load(
  "data/prepData/sampleMuniSpatial_prepData.Rdata"
)

# INITIAL CONDITIONS Z
# AGGREGATE FROM 1000 Sites TO 10 Sites
# transform shares to areas
raster_10_sites$amazonBiomeArea_ha_10Sites <-
  raster_10_sites$share_amazonBiome * raster_10_sites$pixelArea_ha

raster_10_sites$forestArea_1995_ha_10Sites <-
  raster_10_sites$share_forest_1995 * raster_10_sites$pixelArea_ha

raster_10_sites$agriculturalUseArea_1995_ha_10Sites <-
  raster_10_sites$share_agriculturalUse_1995 * raster_10_sites$pixelArea_ha

raster_10_sites$otherArea_1995_ha_10Sites <-
  raster_10_sites$share_other_1995 * raster_10_sites$pixelArea_ha

raster_10_sites$forestArea_2017_ha_10Sites <-
  raster_10_sites$share_forest_2017 * raster_10_sites$pixelArea_ha

raster_10_sites$agriculturalUseArea_2017_ha_10Sites <-
  raster_10_sites$share_agriculturalUse_2017 * raster_10_sites$pixelArea_ha

raster_10_sites$otherArea_2017_ha_10Sites <-
  raster_10_sites$share_other_2017 * raster_10_sites$pixelArea_ha

raster_10_sites$forestArea_2008_ha_10Sites <-
  raster_10_sites$share_forest_2008 * raster_10_sites$pixelArea_ha

raster_10_sites$agriculturalUseArea_2008_ha_10Sites <-
  raster_10_sites$share_agriculturalUse_2008 * raster_10_sites$pixelArea_ha

raster_10_sites$otherArea_2008_ha_10Sites <-
  raster_10_sites$share_other_2008 * raster_10_sites$pixelArea_ha

# select area variables
raster_10_sites <- terra::subset(
  raster_10_sites,
  c(
    "amazonBiomeArea_ha_10Sites",
    "pixelArea_ha",
    "forestArea_1995_ha_10Sites",
    "agriculturalUseArea_1995_ha_10Sites",
    "otherArea_1995_ha_10Sites",
    "forestArea_2017_ha_10Sites",
    "agriculturalUseArea_2017_ha_10Sites",
    "otherArea_2017_ha_10Sites",
    "forestArea_2008_ha_10Sites",
    "agriculturalUseArea_2008_ha_10Sites",
    "otherArea_2008_ha_10Sites"
  )
)

# aggregate from 1000 Sites to 10
raster_10_sites <-
  terra::aggregate(raster_10_sites,
    fact = 14,
    fun = sum,
    na.rm = TRUE
  )


# MAPBIOMAS VARIABLES + AMAZON BIOME + PIXEL AREA (Z_10Sites CONSTRUCTION)
# extract variables as polygons, transform to sf,
# and project data for faster spatial manipulation
calibration_10_sites_model <-
  terra::as.polygons(raster_10_sites, dissolve = FALSE) %>%
  sf::st_as_sf() %>%
  sf::st_transform(5880)




# add id variable
calibration_10_sites_model$id <- seq_len(nrow(calibration_10_sites_model))

# transform share variables in area (ha)
calibration_10_sites_model <-
  calibration_10_sites_model %>%
  dplyr::mutate(
    zbar_1995_10Sites = agriculturalUseArea_1995_ha_10Sites +
      forestArea_1995_ha_10Sites,
    zbar_2017_10Sites = agriculturalUseArea_2017_ha_10Sites +
      forestArea_2017_ha_10Sites,
    zbar_2008_10Sites = agriculturalUseArea_2008_ha_10Sites +
      forestArea_2008_ha_10Sites
  ) %>%
  dplyr::select(
    amazonBiomeArea_ha_10Sites,
    siteArea_ha_10Sites = pixelArea_ha,
    forestArea_1995_ha_10Sites,
    z_1995_10Sites = agriculturalUseArea_1995_ha_10Sites,
    zbar_1995_10Sites,
    forestArea_2017_ha_10Sites,
    z_2017_10Sites = agriculturalUseArea_2017_ha_10Sites,
    zbar_2017_10Sites,
    forestArea_2008_ha_10Sites,
    z_2008_10Sites = agriculturalUseArea_2008_ha_10Sites,
    zbar_2008_10Sites
  )

calibration_10_sites_model <-
  calibration_10_sites_model %>%
  dplyr::filter(amazonBiomeArea_ha_10Sites / siteArea_ha_10Sites >= 0.03)

# add id variable
calibration_10_sites_model$id <- seq_len(nrow(calibration_10_sites_model))



id <- calibration_10_sites_model %>%
  select(id)

st_write(id,
  "data/calibration/hmc/id_10.geojson",
  driver = "GeoJSON",
  delete_dsn = TRUE
)
