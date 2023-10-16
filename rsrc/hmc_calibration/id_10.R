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
tictoc::tic(msg = "calibration_10SitesModel.R script", log = T)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster.10Sites <- terra::rast(list.files(paste(getwd(), "data/calibration/1055SitesModel/aux_tifs", sep = "/"),
                                         pattern = "raster_",
                                         full.names = T))


# MUNI LEVEL SPATIAL SAMPLE
load(here::here("data/calibration/prepData/sampleMuniSpatial_prepData.Rdata"))



# INITIAL CONDITIONS Z -------------------------------------------------------------------------------------------------------------------------------

# AGGREGATE FROM 1000 Sites TO 10 Sites
# transform shares to areas
raster.10Sites$amazonBiomeArea_ha_10Sites <- raster.10Sites$share_amazonBiome*raster.10Sites$pixelArea_ha
raster.10Sites$forestArea_1995_ha_10Sites <- raster.10Sites$share_forest_1995*raster.10Sites$pixelArea_ha
raster.10Sites$agriculturalUseArea_1995_ha_10Sites <- raster.10Sites$share_agriculturalUse_1995*raster.10Sites$pixelArea_ha
raster.10Sites$otherArea_1995_ha_10Sites <- raster.10Sites$share_other_1995*raster.10Sites$pixelArea_ha
raster.10Sites$forestArea_2017_ha_10Sites <- raster.10Sites$share_forest_2017*raster.10Sites$pixelArea_ha
raster.10Sites$agriculturalUseArea_2017_ha_10Sites <- raster.10Sites$share_agriculturalUse_2017*raster.10Sites$pixelArea_ha
raster.10Sites$otherArea_2017_ha_10Sites <- raster.10Sites$share_other_2017*raster.10Sites$pixelArea_ha
raster.10Sites$forestArea_2008_ha_10Sites <- raster.10Sites$share_forest_2008*raster.10Sites$pixelArea_ha
raster.10Sites$agriculturalUseArea_2008_ha_10Sites <- raster.10Sites$share_agriculturalUse_2008*raster.10Sites$pixelArea_ha
raster.10Sites$otherArea_2008_ha_10Sites <- raster.10Sites$share_other_2008*raster.10Sites$pixelArea_ha

# select area variables
raster.10Sites <- terra::subset(raster.10Sites,
                                c("amazonBiomeArea_ha_10Sites", "pixelArea_ha",
                                  "forestArea_1995_ha_10Sites", "agriculturalUseArea_1995_ha_10Sites", "otherArea_1995_ha_10Sites",
                                  "forestArea_2017_ha_10Sites", "agriculturalUseArea_2017_ha_10Sites", "otherArea_2017_ha_10Sites",
                                  "forestArea_2008_ha_10Sites", "agriculturalUseArea_2008_ha_10Sites", "otherArea_2008_ha_10Sites"))

# aggregate from 1000 Sites to 10
raster.10Sites <- terra::aggregate(raster.10Sites, fact = 14, fun = sum, na.rm = T)


# MAPBIOMAS VARIABLES + AMAZON BIOME + PIXEL AREA (Z_10Sites CONSTRUCTION)
# extract variables as polygons, transform to sf, and project data for faster spatial manipulation
calibration.10SitesModel <- terra::as.polygons(raster.10Sites, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)




# add id variable
calibration.10SitesModel$id <- 1:nrow(calibration.10SitesModel)

# transform share variables in area (ha)
calibration.10SitesModel <-
  calibration.10SitesModel %>%
  dplyr::mutate(zbar_1995_10Sites = agriculturalUseArea_1995_ha_10Sites + forestArea_1995_ha_10Sites,
                zbar_2017_10Sites = agriculturalUseArea_2017_ha_10Sites + forestArea_2017_ha_10Sites,
                zbar_2008_10Sites = agriculturalUseArea_2008_ha_10Sites + forestArea_2008_ha_10Sites) %>%
  dplyr::select(amazonBiomeArea_ha_10Sites, siteArea_ha_10Sites = pixelArea_ha,
                forestArea_1995_ha_10Sites,
                z_1995_10Sites = agriculturalUseArea_1995_ha_10Sites, zbar_1995_10Sites,
                forestArea_2017_ha_10Sites,
                z_2017_10Sites = agriculturalUseArea_2017_ha_10Sites, zbar_2017_10Sites,
                forestArea_2008_ha_10Sites,
                z_2008_10Sites = agriculturalUseArea_2008_ha_10Sites, zbar_2008_10Sites)

calibration.10SitesModel <-
  calibration.10SitesModel %>%
  dplyr::filter(amazonBiomeArea_ha_10Sites/siteArea_ha_10Sites >= 0.03)

# add id variable
calibration.10SitesModel$id <- 1:nrow(calibration.10SitesModel)



id <-calibration.10SitesModel %>%
  select(id)

st_write(id, "data/hmc/id_10.geojson", driver = "GeoJSON",delete_dsn = TRUE)
