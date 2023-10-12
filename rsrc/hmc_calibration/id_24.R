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
# Install and load dplyr package
if (!"dplyr" %in% installed.packages()) {
  install.packages("dplyr")
}
library(dplyr)

# Install and load boot package
if (!"boot" %in% installed.packages()) {
  install.packages("boot")
}
library(boot)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::summarize)

# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("code/setup.R")


# START TIMER
tictoc::tic(msg = "calibration_24SitesModel.R script", log = T)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster.24Sites <- terra::rast(list.files(here::here("data/calibration/1055SitesModel/aux_tifs"),
                                         pattern = "raster_",
                                         full.names = T))


# MUNI LEVEL SPATIAL SAMPLE
load(here::here("data/calibration/prepData/sampleMuniSpatial_prepData.Rdata"))




# INITIAL CONDITIONS Z -------------------------------------------------------------------------------------------------------------------------------

# AGGREGATE FROM 1000 SITES TO 43 SITES
# transform shares to areas
raster.24Sites$amazonBiomeArea_ha_24Sites <- raster.24Sites$share_amazonBiome*raster.24Sites$pixelArea_ha
raster.24Sites$forestArea_1995_ha_24Sites <- raster.24Sites$share_forest_1995*raster.24Sites$pixelArea_ha
raster.24Sites$agriculturalUseArea_1995_ha_24Sites <- raster.24Sites$share_agriculturalUse_1995*raster.24Sites$pixelArea_ha
raster.24Sites$otherArea_1995_ha_24Sites <- raster.24Sites$share_other_1995*raster.24Sites$pixelArea_ha
raster.24Sites$forestArea_2017_ha_24Sites <- raster.24Sites$share_forest_2017*raster.24Sites$pixelArea_ha
raster.24Sites$agriculturalUseArea_2017_ha_24Sites <- raster.24Sites$share_agriculturalUse_2017*raster.24Sites$pixelArea_ha
raster.24Sites$otherArea_2017_ha_24Sites <- raster.24Sites$share_other_2017*raster.24Sites$pixelArea_ha

# select area variables
raster.24Sites <- terra::subset(raster.24Sites, c("amazonBiomeArea_ha_24Sites", "pixelArea_ha",
                                                  "forestArea_1995_ha_24Sites", "agriculturalUseArea_1995_ha_24Sites", "otherArea_1995_ha_24Sites",
                                                  "forestArea_2017_ha_24Sites", "agriculturalUseArea_2017_ha_24Sites", "otherArea_2017_ha_24Sites"))

# aggregate from 1000 sites to 24
raster.24Sites <- terra::aggregate(raster.24Sites, fact = 8, fun = sum, na.rm = T)



# extract variables as polygons, transform to sf, and project data for faster spatial manipulation
calibration.24SitesModel <- terra::as.polygons(raster.24Sites, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)




# add id variable
calibration.24SitesModel$id <- 1:nrow(calibration.24SitesModel)

# transform share variables in area (ha)
calibration.24SitesModel <-
  calibration.24SitesModel %>%
  dplyr::mutate(zbar_1995_24Sites = agriculturalUseArea_1995_ha_24Sites + forestArea_1995_ha_24Sites,
                zbar_2017_24Sites = agriculturalUseArea_2017_ha_24Sites + forestArea_2017_ha_24Sites) %>%
  dplyr::select(amazonBiomeArea_ha_24Sites, siteArea_ha_24Sites = pixelArea_ha,
                forestArea_1995_ha_24Sites,
                z_1995_24Sites = agriculturalUseArea_1995_ha_24Sites, zbar_1995_24Sites,
                forestArea_2017_ha_24Sites,
                z_2017_24Sites = agriculturalUseArea_2017_ha_24Sites, zbar_2017_24Sites)


# remove sites with less than 2% of its are intersecting with the amazon biome
calibration.24SitesModel <-
  calibration.24SitesModel %>%
  dplyr::filter(amazonBiomeArea_ha_24Sites/siteArea_ha_24Sites >= 0.03)

# add id variable
calibration.24SitesModel$id <- 1:nrow(calibration.24SitesModel)



id <-calibration.24SitesModel %>%
  select(id)

st_write(id, "data/hmc/id_24.geojson", driver = "GeoJSON",delete_dsn = TRUE)
