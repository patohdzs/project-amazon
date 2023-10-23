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
tictoc::tic(msg = "calibration_78SitesModel.R script", log = T)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))


# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster.78Sites <-
  terra::rast(list.files(
    here::here("data/calibration/1055SitesModel/aux_tifs"),
    pattern = "raster_",
    full.names = T
  ))

# MUNI LEVEL SPATIAL SAMPLE
load(here::here(
  "data/calibration/prepData/sampleMuniSpatial_prepData.Rdata"
))


# AGGREGATE FROM 1000 SITES TO 43 SITES
# transform shares to areas
raster.78Sites$amazonBiomeArea_ha_78Sites <-
  raster.78Sites$share_amazonBiome * raster.78Sites$pixelArea_ha
raster.78Sites$forestArea_1995_ha_78Sites <-
  raster.78Sites$share_forest_1995 * raster.78Sites$pixelArea_ha
raster.78Sites$agriculturalUseArea_1995_ha_78Sites <-
  raster.78Sites$share_agriculturalUse_1995 * raster.78Sites$pixelArea_ha
raster.78Sites$otherArea_1995_ha_78Sites <-
  raster.78Sites$share_other_1995 * raster.78Sites$pixelArea_ha
raster.78Sites$forestArea_2017_ha_78Sites <-
  raster.78Sites$share_forest_2017 * raster.78Sites$pixelArea_ha
raster.78Sites$agriculturalUseArea_2017_ha_78Sites <-
  raster.78Sites$share_agriculturalUse_2017 * raster.78Sites$pixelArea_ha
raster.78Sites$otherArea_2017_ha_78Sites <-
  raster.78Sites$share_other_2017 * raster.78Sites$pixelArea_ha

# select area variables
raster.78Sites <-
  terra::subset(
    raster.78Sites,
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
raster.78Sites <-
  terra::aggregate(raster.78Sites,
                   fact = 4,
                   fun = sum,
                   na.rm = T)



# extract variables as polygons, transform to sf, and project data for faster spatial manipulation
calibration.78SitesModel <-
  terra::as.polygons(raster.78Sites, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)




# add id variable
calibration.78SitesModel$id <- 1:nrow(calibration.78SitesModel)

# transform share variables in area (ha)
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::mutate(
    zbar_1995_78Sites = agriculturalUseArea_1995_ha_78Sites + forestArea_1995_ha_78Sites,
    zbar_2017_78Sites = agriculturalUseArea_2017_ha_78Sites + forestArea_2017_ha_78Sites
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
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::filter(amazonBiomeArea_ha_78Sites / siteArea_ha_78Sites >= 0.03)

# add id variable
calibration.78SitesModel$id <- 1:nrow(calibration.78SitesModel)



id <- calibration.78SitesModel %>%
  select(id)

st_write(id,
         "data/hmc/id_78.geojson",
         driver = "GeoJSON",
         delete_dsn = TRUE)
