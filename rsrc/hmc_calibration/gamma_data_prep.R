
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
conflicts_prefer(dplyr::filter)

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

library(MASS)
conflicts_prefer(dplyr::select)
# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "calibration_25SitesModel.R script", log = T)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = paste(getwd(), "data", "_temp", sep = "/"))




# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster.25Sites <- terra::rast(list.files(paste(getwd(), "data/calibration/1055SitesModel/aux_tifs", sep = "/"),
                                         pattern = "raster_",
                                         full.names = T))


# MUNI LEVEL SPATIAL SAMPLE
load(paste(getwd(), "data/calibration/prepData/sampleMuniSpatial_prepData.Rdata", sep = "/"))


my_data <- read_excel("C:/Users/pengyu/Desktop/code_data_20230628/data/ipeadata[21-08-2023-01-28].xls")
my_data$muni_code <- as.numeric(my_data$muni_code)


# INITIAL CONDITIONS Z -------------------------------------------------------------------------------------------------------------------------------

# AGGREGATE FROM 1000 Sites TO 25 Sites
# transform shares to areas
raster.25Sites$amazonBiomeArea_ha_25Sites <- raster.25Sites$share_amazonBiome*raster.25Sites$pixelArea_ha
raster.25Sites$forestArea_1995_ha_25Sites <- raster.25Sites$share_forest_1995*raster.25Sites$pixelArea_ha
raster.25Sites$agriculturalUseArea_1995_ha_25Sites <- raster.25Sites$share_agriculturalUse_1995*raster.25Sites$pixelArea_ha
raster.25Sites$otherArea_1995_ha_25Sites <- raster.25Sites$share_other_1995*raster.25Sites$pixelArea_ha
raster.25Sites$forestArea_2017_ha_25Sites <- raster.25Sites$share_forest_2017*raster.25Sites$pixelArea_ha
raster.25Sites$agriculturalUseArea_2017_ha_25Sites <- raster.25Sites$share_agriculturalUse_2017*raster.25Sites$pixelArea_ha
raster.25Sites$otherArea_2017_ha_25Sites <- raster.25Sites$share_other_2017*raster.25Sites$pixelArea_ha

# select area variables
raster.25Sites <- terra::subset(raster.25Sites,
                                c("amazonBiomeArea_ha_25Sites", "pixelArea_ha",
                                  "forestArea_1995_ha_25Sites", "agriculturalUseArea_1995_ha_25Sites", "otherArea_1995_ha_25Sites",
                                  "forestArea_2017_ha_25Sites", "agriculturalUseArea_2017_ha_25Sites", "otherArea_2017_ha_25Sites"))

# aggregate from 1000 Sites to 25
raster.25Sites <- terra::aggregate(raster.25Sites, fact = 8, fun = sum, na.rm = T)

# extract variables as polygons, transform to sf, and project data for faster spatial manipulation
calibration.25SitesModel <- terra::as.polygons(raster.25Sites, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)

# transform share aggregate in area (ha)
calibration.25SitesModel <-
  calibration.25SitesModel %>%
  dplyr::mutate(zbar_1995_25Sites = agriculturalUseArea_1995_ha_25Sites + forestArea_1995_ha_25Sites,
                zbar_2017_25Sites = agriculturalUseArea_2017_ha_25Sites + forestArea_2017_ha_25Sites) %>%
  dplyr::select(amazonBiomeArea_ha_25Sites, siteArea_ha_25Sites = pixelArea_ha,
                forestArea_1995_ha_25Sites,
                z_1995_25Sites = agriculturalUseArea_1995_ha_25Sites, zbar_1995_25Sites,
                forestArea_2017_ha_25Sites,
                z_2017_25Sites = agriculturalUseArea_2017_ha_25Sites, zbar_2017_25Sites)

# remove Sites with less than 1% of its are intersecting with the amazon biome
calibration.25SitesModel <-
  calibration.25SitesModel %>%
  dplyr::filter(amazonBiomeArea_ha_25Sites/siteArea_ha_25Sites >= 0.01)

# add id variable
calibration.25SitesModel$id <- 1:nrow(calibration.25SitesModel)



# PARAMETER GAMMA ------------------------------------------------------------------------------------------------------------------------------------


# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/calibration/prepData/muniTheta_prepData_gamma.Rdata")




muniTheta.prepData<-muniTheta.prepData %>%
  dplyr::mutate(co2e_ha_2017 = (agb_2017/2)*(44/12))



muniTheta.prepData.filter<- muniTheta.prepData %>%
  filter(!is.na(co2e_ha_2017))

#reg.share.2017 <-
#  muniTheta.prepData  %>%
#  lm(formula = log(share)  ~ log(lat)+log(lon), na.action = na.exclude)

#summary(reg.share.2017)


#residuals_values <- residuals(reg.share.2017)
#residuals_values_df<-data.frame(residuals_values)


#muniTheta.prepData<-muniTheta.prepData%>%
#  dplyr::mutate(residuals=residuals_values_df$residuals_values)


# DATA INPUT (2017)
# load pixel sample with biomass data
# Load pixelBiomass2017_prepData.Rdata




reg.gamma.2017 <-
  muniTheta.prepData.filter  %>%
  lm(formula = log(co2e_ha_2017)  ~ log(historical_precip) + log(historical_temp) +log(lat)+log(lon), na.action = na.exclude)

reg_summary <- summary(reg.gamma.2017)
regressor_df <- as.data.frame(reg.gamma.2017$model[-1])




new_df <- muniTheta.prepData %>%
  select(historical_precip,
         historical_temp, lat,lon)

new_df <- new_df %>%
  mutate(log_historical_precip = log(historical_precip))
new_df <- new_df %>%
  mutate(log_historical_temp = log(historical_temp))
new_df <- new_df %>%
  mutate(log_lat = log(lat))
new_df <- new_df %>%
  mutate(log_lon = log(lon))
new_df<-new_df %>%
  mutate(log_historical_precip = (log_historical_precip-mean(regressor_df$`log(historical_precip)`))/ sd(regressor_df$`log(historical_precip)`)) %>%
  mutate(log_historical_temp = (log_historical_temp-mean(regressor_df$`log(historical_temp)`))/ sd(regressor_df$`log(historical_temp)`)) %>%
  mutate(log_lat = (log_lat-mean(regressor_df$`log(lat)`))/ sd(regressor_df$`log(lat)`)) %>%
  mutate(log_lon = (log_lon-mean(regressor_df$`log(lon)`))/ sd(regressor_df$`log(lon)`))



new_df <- cbind(1, new_df)
#new_df <- as.data.frame(new_df)
new_df <- new_df %>%
  select(1, log_historical_precip,log_historical_temp,log_lat,log_lon)


st_write(new_df, "data/hmc/data_gamma.geojson", driver = "GeoJSON",delete_dsn = TRUE)
