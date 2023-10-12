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
  dplyr::filter(amazonBiomeArea_ha_25Sites/siteArea_ha_25Sites >= 0.03)

# add id variable
calibration.25SitesModel$id <- 1:nrow(calibration.25SitesModel)





# PARAMETER THETA ------------------------------------------------------------------------------------------------------------------------------------

# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/calibration/prepData/muniTheta_prepData.Rdata")

# load cattle price series
load("data/calibration/prepData/seriesPriceCattle_prepData.Rdata")


# DATA MANIPULATION

# EXTRACT AVERAGE 2017 PRICE (use real prices because it is normalized to 2017 )
aux.price.2017 <-
  seriesPriceCattle.prepData %>%
  dplyr::filter(year == 2017) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(mean_price_2017 = mean(price_real_mon_cattle)/3.192) %>% # BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
  dplyr::pull(mean_price_2017)


# REGRESSION - CATTLE VALUE (2017)

# cattle value per ha
#reg.cattleValueperHa.2017 <-
# muniTheta.prepData %>%
#lm(formula = cattleSlaughter_valuePerHa_2017 ~ pasture_area_2017 + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
#    lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = pasture_area_2017)

#muniTheta.prepData <- muniTheta.prepData[-c(142, 106, 112), ]


my_data <- read_excel("C:/Users/pengyu/Desktop/code_data_20230628/data/ipeadata[21-08-2023-01-28].xls")
my_data$muni_code <- as.numeric(my_data$muni_code)

# Remove rows from attribute data
a<-muniTheta.prepData
muniTheta.prepData_data <- as.data.frame(muniTheta.prepData)  # Convert to regular dataframe
muniTheta.prepData_data <- muniTheta.prepData_data[-c(142, 106, 112), ]

# Remove geometries
geo_backup <- st_geometry(muniTheta.prepData)
geo_backup <- geo_backup[-c(142, 106, 112)]


predicted_values <- read_excel("C:/Users/pengyu/Desktop/code_data_20230628/data/farm_gate_price.xlsx")

# Combine back into an sf object
muniTheta.prepData <- st_sf(muniTheta.prepData_data, geometry = geo_backup)

# 2. Merging the cleaned muniTheta.prepData with my_data

# Convert to non-spatial dataframe for the merge
muniTheta_no_geo <- as.data.frame(muniTheta.prepData)

# Perform the merge
merged_data <- left_join(muniTheta_no_geo, my_data, by = "muni_code")

# Reattach the geometry
merged_data_sf <- st_sf(merged_data, geometry = geo_backup)

muniTheta.prepData<-merged_data_sf



merged_data <- muniTheta.prepData %>%
  left_join(predicted_values, by = "muni_code") %>%
  mutate(cattleSlaughter_farmGatePrice_2017 = ifelse(is.na(cattleSlaughter_farmGatePrice_2017),
                                                     average_weighted_price,
                                                     cattleSlaughter_farmGatePrice_2017))
# select(-predicted_value_column_name)  # Remove the additional column from the result


muniTheta.prepData<-merged_data


muniTheta.prepData<- muniTheta.prepData %>%
  filter(!is.na(distance))







muniTheta.prepData_filtered <- muniTheta.prepData %>%
  filter(cattleSlaughter_valuePerHa_2017 > 0) # Exclude zeros
#muniTheta.prepData_filtered <- na.omit(muniTheta.prepData_filtered)


b<-muniTheta.prepData
c<-muniTheta.prepData_filtered

new_df <- muniTheta.prepData %>%
  select(historical_precip,
         historical_temp, lat,cattleSlaughter_farmGatePrice_2017,distance,zbar_2017_muni)

# drop spatial feature
#new_df <-
#  new_df %>%
#  sf::st_drop_geometry()
new_df <- new_df %>%
  mutate(`I(historical_temp^2)` = historical_temp^2)
new_df <- new_df %>%
  mutate(`I(lat^2)` = lat^2)
new_df <- cbind(1, new_df)
#new_df <- as.data.frame(new_df)
new_df <- new_df %>%
  select(1, historical_precip,historical_temp,I.historical_temp.2.,lat,I.lat.2.,cattleSlaughter_farmGatePrice_2017,distance,zbar_2017_muni)

new_df <- new_df %>%
  mutate(X1=X1) %>%
  mutate(historical_precip = (historical_precip-mean(muniTheta.prepData_filtered$historical_precip))/sd(muniTheta.prepData_filtered$historical_precip)) %>%
  mutate(historical_temp = (historical_temp-mean(muniTheta.prepData_filtered$historical_temp))/sd(muniTheta.prepData_filtered$historical_temp)) %>%
  mutate(I.historical_temp.2. = (I.historical_temp.2.-mean(muniTheta.prepData_filtered$historical_temp^2))/sd(muniTheta.prepData_filtered$historical_temp^2)) %>%
  mutate(lat = (lat-mean(muniTheta.prepData_filtered$lat))/sd(muniTheta.prepData_filtered$lat)) %>%
  mutate(I.lat.2. =(I.lat.2.-mean(muniTheta.prepData_filtered$lat^2))/sd(muniTheta.prepData_filtered$lat^2)) %>%
  mutate(cattleSlaughter_farmGatePrice_2017=(cattleSlaughter_farmGatePrice_2017-mean(muniTheta.prepData_filtered$cattleSlaughter_farmGatePrice_2017))/sd(muniTheta.prepData_filtered$cattleSlaughter_farmGatePrice_2017))%>%
  mutate(distance = (distance-mean(muniTheta.prepData_filtered$distance))/sd(muniTheta.prepData_filtered$distance)) 





st_write(new_df, "data/HMC_norm/data_theta.geojson", driver = "GeoJSON",delete_dsn = TRUE)
