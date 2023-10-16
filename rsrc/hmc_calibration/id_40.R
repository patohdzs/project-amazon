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
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "calibration_40SitesModel.R script", log = T)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster.40Sites <- terra::rast(list.files(here::here("data/calibration/1055SitesModel/aux_tifs"),
                                         pattern = "raster_",
                                         full.names = T))


# MUNI LEVEL SPATIAL SAMPLE
load(here::here("data/calibration/prepData/sampleMuniSpatial_prepData.Rdata"))




# INITIAL CONDITIONS Z -------------------------------------------------------------------------------------------------------------------------------

# AGGREGATE FROM 1000 SITES TO 43 SITES
# transform shares to areas
raster.40Sites$amazonBiomeArea_ha_40Sites <- raster.40Sites$share_amazonBiome*raster.40Sites$pixelArea_ha
raster.40Sites$forestArea_1995_ha_40Sites <- raster.40Sites$share_forest_1995*raster.40Sites$pixelArea_ha
raster.40Sites$agriculturalUseArea_1995_ha_40Sites <- raster.40Sites$share_agriculturalUse_1995*raster.40Sites$pixelArea_ha
raster.40Sites$otherArea_1995_ha_40Sites <- raster.40Sites$share_other_1995*raster.40Sites$pixelArea_ha
raster.40Sites$forestArea_2017_ha_40Sites <- raster.40Sites$share_forest_2017*raster.40Sites$pixelArea_ha
raster.40Sites$agriculturalUseArea_2017_ha_40Sites <- raster.40Sites$share_agriculturalUse_2017*raster.40Sites$pixelArea_ha
raster.40Sites$otherArea_2017_ha_40Sites <- raster.40Sites$share_other_2017*raster.40Sites$pixelArea_ha

# select area variables
raster.40Sites <- terra::subset(raster.40Sites, c("amazonBiomeArea_ha_40Sites", "pixelArea_ha",
                                                  "forestArea_1995_ha_40Sites", "agriculturalUseArea_1995_ha_40Sites", "otherArea_1995_ha_40Sites",
                                                  "forestArea_2017_ha_40Sites", "agriculturalUseArea_2017_ha_40Sites", "otherArea_2017_ha_40Sites"))

# aggregate from 1000 sites to 43
raster.40Sites <- terra::aggregate(raster.40Sites, fact = 6, fun = sum, na.rm = T)



# extract variables as polygons, transform to sf, and project data for faster spatial manipulation
calibration.40SitesModel <- terra::as.polygons(raster.40Sites, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)




# add id variable
calibration.40SitesModel$id <- 1:nrow(calibration.40SitesModel)

# transform share variables in area (ha)
calibration.40SitesModel <-
  calibration.40SitesModel %>%
  dplyr::mutate(zbar_1995_40Sites = agriculturalUseArea_1995_ha_40Sites + forestArea_1995_ha_40Sites,
                zbar_2017_40Sites = agriculturalUseArea_2017_ha_40Sites + forestArea_2017_ha_40Sites) %>%
  dplyr::select(amazonBiomeArea_ha_40Sites, siteArea_ha_40Sites = pixelArea_ha,
                forestArea_1995_ha_40Sites,
                z_1995_40Sites = agriculturalUseArea_1995_ha_40Sites, zbar_1995_40Sites,
                forestArea_2017_ha_40Sites,
                z_2017_40Sites = agriculturalUseArea_2017_ha_40Sites, zbar_2017_40Sites)


# remove sites with less than 2% of its are intersecting with the amazon biome
calibration.40SitesModel <-
  calibration.40SitesModel %>%
  dplyr::filter(amazonBiomeArea_ha_40Sites/siteArea_ha_40Sites >= 0.03)

# add id variable
calibration.40SitesModel$id <- 1:nrow(calibration.40SitesModel)



id <-calibration.40SitesModel %>%
  select(id)

st_write(id, "data/hmc/id_40.geojson", driver = "GeoJSON",delete_dsn = TRUE)
