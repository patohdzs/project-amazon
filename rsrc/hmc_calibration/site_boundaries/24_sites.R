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
tictoc::tic(msg = "24_sites.R script", log = TRUE)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tmpdir = "data/_temp",
                      timer  = T)





# DATA INPUT
# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster_24_sites <- terra::rast(list.files("data/calibration/1043SitesModel/aux_tifs",
                                           pattern = "raster_",
                                           full.names = T))



# MUNI LEVEL SPATIAL SAMPLE
load(
  "data/prepData/sampleMuniSpatial_prepData.Rdata"
)




# INITIAL CONDITIONS Z
# AGGREGATE FROM 1000 SITES TO 43 SITES
# transform shares to areas
raster_24_sites$amazonBiomeArea_ha_24Sites <-
  raster_24_sites$share_amazonBiome * raster_24_sites$pixelArea_ha
raster_24_sites$forestArea_1995_ha_24Sites <-
  raster_24_sites$share_forest_1995 * raster_24_sites$pixelArea_ha
raster_24_sites$agriculturalUseArea_1995_ha_24Sites <-
  raster_24_sites$share_agriculturalUse_1995 * raster_24_sites$pixelArea_ha
raster_24_sites$otherArea_1995_ha_24Sites <-
  raster_24_sites$share_other_1995 * raster_24_sites$pixelArea_ha
raster_24_sites$forestArea_2017_ha_24Sites <-
  raster_24_sites$share_forest_2017 * raster_24_sites$pixelArea_ha
raster_24_sites$agriculturalUseArea_2017_ha_24Sites <-
  raster_24_sites$share_agriculturalUse_2017 * raster_24_sites$pixelArea_ha
raster_24_sites$otherArea_2017_ha_24Sites <-
  raster_24_sites$share_other_2017 * raster_24_sites$pixelArea_ha

# select area variables
raster_24_sites <-
  terra::subset(
    raster_24_sites,
    c(
      "amazonBiomeArea_ha_24Sites",
      "pixelArea_ha",
      "forestArea_1995_ha_24Sites",
      "agriculturalUseArea_1995_ha_24Sites",
      "otherArea_1995_ha_24Sites",
      "forestArea_2017_ha_24Sites",
      "agriculturalUseArea_2017_ha_24Sites",
      "otherArea_2017_ha_24Sites"
    )
  )

# aggregate from 1000 sites to 24
raster_24_sites <-
  terra::aggregate(raster_24_sites,
    fact = 8,
    fun = sum,
    na.rm = TRUE
  )



# extract variables as polygons, transform to sf,
# and project data for faster spatial manipulation
calibration_24_sites_model <-
  terra::as.polygons(raster_24_sites, dissolve = FALSE) %>%
  sf::st_as_sf() %>%
  sf::st_transform(5880)




# add id variable
calibration_24_sites_model$id <- seq_len(nrow(calibration_24_sites_model))

# transform share variables in area (ha)
calibration_24_sites_model <-
  calibration_24_sites_model %>%
  dplyr::mutate(
    zbar_1995_24Sites = agriculturalUseArea_1995_ha_24Sites +
      forestArea_1995_ha_24Sites,
    zbar_2017_24Sites = agriculturalUseArea_2017_ha_24Sites +
      forestArea_2017_ha_24Sites
  ) %>%
  dplyr::select(
    amazonBiomeArea_ha_24Sites,
    siteArea_ha_24Sites = pixelArea_ha,
    forestArea_1995_ha_24Sites,
    z_1995_24Sites = agriculturalUseArea_1995_ha_24Sites,
    zbar_1995_24Sites,
    forestArea_2017_ha_24Sites,
    z_2017_24Sites = agriculturalUseArea_2017_ha_24Sites,
    zbar_2017_24Sites
  )


# remove sites with less than 2% of its are intersecting with the amazon biome
calibration_24_sites_model <-
  calibration_24_sites_model %>%
  dplyr::filter(amazonBiomeArea_ha_24Sites / siteArea_ha_24Sites >= 0.03)

# add id variable
calibration_24_sites_model$id <- seq_len(nrow(calibration_24_sites_model))

id <- calibration_24_sites_model %>%
  select(id)

st_write(id,
  "data/calibration/hmc/id_24.geojson",
  driver = "GeoJSON",
  delete_dsn = TRUE
)
