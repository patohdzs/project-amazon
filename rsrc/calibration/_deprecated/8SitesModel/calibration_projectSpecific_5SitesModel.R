
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (5 SITES MODEL)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# GROUNDHOG (REPRODUCIBILITY SOLUTION TO HANDLING DIFFERENT VERSIONS OF R AND ITS PACKAGES)

# check if groundhog is installed and load it
if ("groundhog" %in% installed.packages()) {
  library("groundhog")
} else {
  install.packages("groundhog")
  library("groundhog")
}

# define date of reference to load all packages
groundhog.date <- "2022-04-01"

# guarantee version 1.5 of groundhog is being used
groundhog::meta.groundhog(date = "2022-04-01")


# HERE
groundhog::groundhog.library("here", groundhog.date) # load package here


# TICTOC
groundhog::groundhog.library("tictoc", groundhog.date) # load package tictoc


# DECLARE LOCATION OF CURRENT SCRIPT TO SET UP PROJECT ROOT CORRECTLY
here::i_am("code/projectSpecific/5SitesModel/calibration_projectSpecific_5SitesModel.R", uuid = "0fef2915-0a76-4b7f-9a2c-c0b9830c058b")


# START TIME
tictoc::tic(msg = "calibration_projectSpecific_5SitesModel script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary
groundhog::groundhog.library("sf", groundhog.date)  # manipulate spatial data (vector format)
groundhog::groundhog.library("terra", groundhog.date)  # to calculate probability transition matrix


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster.1000Sites <- terra::rast(list.files(here::here("data/projectSpecific/1000SitesModel/aux_tifs"),
                                           pattern = "raster_",
                                           full.names = T))


# CALIBRATION GLOBAL MODEL
load(here::here("data/projectSpecific/globalModel/calibration_globalModel.Rdata"))

# clean enviroment
rm(matrixTransition.2prices, matrixTransition.2prices.alt, matrixTransition.4prices, matrixTransition.5prices)


# MUNI LEVEL SPATIAL SAMPLE
load(here::here("data/projectSpecific/prepData/sampleMuniSpatial_prepData.Rdata"))





# INITIAL CONDITIONS Z -------------------------------------------------------------------------------------------------------------------------------

# AGGREGATE FROM 1000 SITES TO 5 SITES
# transform shares to areas
raster.1000Sites$amazonBiomeArea_ha_5Sites <- raster.1000Sites$share_amazonBiome*raster.1000Sites$pixelArea_ha
raster.1000Sites$forestArea_1985_ha_5Sites <- raster.1000Sites$share_forest_1985*raster.1000Sites$pixelArea_ha
raster.1000Sites$agriculturalUseArea_1985_ha_5Sites <- raster.1000Sites$share_agriculturalUse_1985*raster.1000Sites$pixelArea_ha
raster.1000Sites$otherArea_1985_ha_5Sites <- raster.1000Sites$share_other_1985*raster.1000Sites$pixelArea_ha
raster.1000Sites$forestArea_2017_ha_5Sites <- raster.1000Sites$share_forest_2017*raster.1000Sites$pixelArea_ha
raster.1000Sites$agriculturalUseArea_2017_ha_5Sites <- raster.1000Sites$share_agriculturalUse_2017*raster.1000Sites$pixelArea_ha
raster.1000Sites$otherArea_2017_ha_5Sites <- raster.1000Sites$share_other_2017*raster.1000Sites$pixelArea_ha

# select area variables
raster.1000Sites <- terra::subset(raster.1000Sites, c(9:15))

# aggregate in two steps because using fact =12 generated a significant distortion while fact=4 followed by fact=3 was very close (in terms of total area compared to the 1000 sites)
raster.1000Sites.2 <- terra::aggregate(raster.1000Sites, fact = 4, fun = sum, na.rm = T)
raster.5Sites <- terra::aggregate(raster.1000Sites.2, fact = 5, fun = sum, na.rm = T)

# extract variables as polygons, transform to sf, and project data for faster spatial manipulation
calibration.5SitesModel <- terra::as.polygons(raster.5Sites, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)

# add id variable
calibration.5SitesModel$id <- 1:length(calibration.5SitesModel$geometry)

# transform share aggregate in area (ha)
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(zbar_1985_5Sites = agriculturalUseArea_1985_ha_5Sites + forestArea_1985_ha_5Sites,
                zbar_2017_5Sites = agriculturalUseArea_2017_ha_5Sites + forestArea_2017_ha_5Sites) %>%
  dplyr::select(id, amazonBiomeArea_ha_5Sites,
                forestArea_1985_ha_5Sites, otherArea_1985_ha_5Sites,
                z_1985_5Sites = agriculturalUseArea_1985_ha_5Sites, zbar_1985_5Sites,
                forestArea_2017_ha_5Sites, otherArea_2017_ha_5Sites,
                z_2017_5Sites = agriculturalUseArea_2017_ha_5Sites, zbar_2017_5Sites)



# PARAMETER GAMMA ------------------------------------------------------------------------------------------------------------------------------------

# # DATA INPUT
# # load pixel sample with biomass data
# load(here::here("data/projectSpecific/prepData/pixelBiomass_prepData.Rdata"))
#
# # select minicells of primary forest with co2 information and transform to spatial points
# pixelBiomass.prepData <-
#   pixelBiomass.prepData %>%
#   dplyr::filter(mapbiomas_classAgg == "primaryForest") %>%
#   sf::st_transform(sf::st_crs(calibration.5SitesModel)) %>%
#   dplyr::mutate(co2e_ha = (agb/2)*(44/12)) %>%
#   dplyr::select(co2e_ha)
#
# # match minicells with sites
# aux.gamma <- sf::st_join(calibration.5SitesModel, pixelBiomass.prepData)
#
# # calculate average carbon density on primary forest areas by site
# aux.gamma <-
#   aux.gamma %>%
#   sf::st_drop_geometry() %>%
#   dplyr::group_by(id) %>%
#   summarise(gamma_5Sites = mean(co2e_ha, na.rm = T))
#
#
# # add gamma_5Sites to spatial variables
# calibration.5SitesModel <- left_join(calibration.5SitesModel, aux.gamma)
#
# # clean environment
# rm(pixelBiomass.prepData, aux.gamma, raster.variables)
#
# # add alternative value for gamma with a 10.1% correction factor for small trees and lianas (see Malhi et al (2006))
# calibration.5SitesModel <-
#   calibration.5SitesModel %>%
#   dplyr::mutate(gamma_alt_5Sites = gamma_5Sites*1.101)



# DATA INPUT (2017)
# load pixel sample with biomass data
load(here::here("data/projectSpecific/prepData/pixelBiomass2017_prepData.Rdata"))

# select minicells of primary forest with co2 information and transform to spatial points
pixelBiomass2017.prepData <-
  pixelBiomass2017.prepData %>%
  dplyr::filter(mapbiomas_classAgg == "primaryForest") %>%
  sf::st_transform(sf::st_crs(calibration.5SitesModel)) %>%
  dplyr::mutate(co2e_ha_2017 = (agb_2017/2)*(44/12),
                co2eSD_ha_2017 = (agbSD_2017/2)*(44/12)) %>%
  dplyr::select(lon, lat, co2e_ha_2017, co2eSD_ha_2017) %>%
  dplyr::filter(co2e_ha_2017 > 0, co2eSD_ha_2017 > 0, !is.na(co2e_ha_2017), !is.na(co2eSD_ha_2017))

# cattle value per ha
reg.carbonUncertainty <-
  pixelBiomass2017.prepData %>%
  lm(formula = log(co2eSD_ha_2017) ~ log(co2e_ha_2017) + lon*lat + I(lon^2) + I(lat^2))

# regression results
summary(reg.carbonUncertainty)

# extract fitted values
pixelBiomass2017.prepData <-
  pixelBiomass2017.prepData %>%
  dplyr::mutate(co2eSD_ha_2017_fitted = predict(reg.carbonUncertainty, .),
                co2eSD_ha_2017_fitted = exp(co2eSD_ha_2017_fitted))

# match minicells with sites
aux.gamma2017 <- sf::st_join(calibration.5SitesModel, pixelBiomass2017.prepData)

# calculate average carbon density on primary forest areas by site
aux.gamma2017 <- aux.gamma2017 %>% sf::st_drop_geometry() %>% dplyr::group_by(id) %>% summarise(gamma_5Sites = mean(co2e_ha_2017, na.rm = T),
                                                                                                co2e_ha_2017 = mean(co2e_ha_2017, na.rm = T),
                                                                                                gammaSD_5Sites = mean(co2eSD_ha_2017_fitted, na.rm = T),
                                                                                                co2eSD_ha_2017 = mean(co2eSD_ha_2017, na.rm = T),
                                                                                                lon = mean(lon, na.rm = T),
                                                                                                lat = mean(lat, na.rm = T))

# add gamma_5Sites to spatial variables
calibration.5SitesModel <- left_join(calibration.5SitesModel, aux.gamma2017)

# clean environment
rm(pixelBiomass2017.prepData, aux.gamma2017)

calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(gammaSD_5Sites = predict(reg.carbonUncertainty, .),
                gammaSD_5Sites = exp(gammaSD_5Sites))

# clean environment
rm(reg.carbonUncertainty)







# PARAMETERS A AND B ---------------------------------------------------------------------------------------------------------------------------------

# estimate of a same as in the global model >
# estimate of b given a such that b_i/a is the average carbon in primary forest
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(a_5Sites = calibration.globalModel$a_global,
                b_5Sites = a_5Sites*gamma_5Sites*zbar_2017_5Sites)





# PARAMETER K ----------------------------------------------------------------------------------------------------------------------------------------

# estimate of k same as in the global model
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(k_5Sites = calibration.globalModel$k_global)





# PARAMETER ZETA ----------------------------------------------------------------------------------------------------------------------------------------

# estimate of zeta same as in the global model
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(zeta_5Sites = calibration.globalModel$zeta_global,
                zeta_alt_5Sites = calibration.globalModel$zeta_alt_global)




# PARAMETER RHO ----------------------------------------------------------------------------------------------------------------------------------------

# estimate of rho same as in the global model
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(rho_5Sites = calibration.globalModel$rho_global)





# PARAMETER R ----------------------------------------------------------------------------------------------------------------------------------------

# estimate of r same as in the global model
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(r_5Sites = calibration.globalModel$r_global)





# INITIAL CONDITIONS X -------------------------------------------------------------------------------------------------------------------------------

# x_2017_5Sites estimated as in the old way of global model, just considering the stock of carbon stored in forest areas assuming that all forests are primary
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(x_2017_5Sites = gamma_5Sites*forestArea_2017_ha_5Sites)





# INITIAL CONDITIONS C -------------------------------------------------------------------------------------------------------------------------------

# estimates of C same as in the global model
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(C_2017_5Sites = calibration.globalModel$C_2017_global,
                C0bern_2017_5Sites = calibration.globalModel$C0bern_2017_global,
                C1bern_2017_5Sites = calibration.globalModel$C1bern_2017_global,
                C2bern_2017_5Sites = calibration.globalModel$C2bern_2017_global,
                C3bern_2017_5Sites = calibration.globalModel$C3bern_2017_global)




# PARAMETER U ----------------------------------------------------------------------------------------------------------------------------------------

# estimates of u same as in the global model
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(u_5Sites = calibration.globalModel$u_global,
                u_bern_5Sites = calibration.globalModel$u_bern_global)





# PARAMETER THETA ------------------------------------------------------------------------------------------------------------------------------------

# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/projectSpecific/prepData/muniTheta_prepData.Rdata")

# load cattle price series
load("data/projectSpecific/prepData/seriesPriceCattle_prepData.Rdata")


# DATA MANIPULATION

# EXTRACT AVERAGE 2017 PRICE (use real prices because it is normalized to 2017 )
aux.price.2017 <-
  seriesPriceCattle.prepData %>%
  dplyr::mutate(price_real_mon_cattle = price_real_mon_cattle/3.1966) %>%
  dplyr::filter(year == 2017) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(mean_price_2017 = mean(price_real_mon_cattle)) %>%
  dplyr::pull(mean_price_2017)

# EXTRACT AVERAGE 2006 PRICE (use real prices because it is normalized to 2006)
aux.price.2006 <-
  seriesPriceCattle.prepData %>%
  dplyr::mutate(price_real_mon_cattle = price_real2006_mon_cattle/2.1761) %>%
  dplyr::filter(year == 2006) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(mean_price_2006 = mean(price_real_mon_cattle)) %>%
  dplyr::pull(mean_price_2006)


# REGRESSION - CATTLE VALUE (2017)

# cattle value per ha
reg.cattleValueperHa <-
  muniTheta.prepData %>%
  lm(formula = cattleSlaughter_value_ha ~ pastureArea_value + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
       lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = pastureArea_value)

# regression results
summary(reg.cattleValueperHa)

# extract fitted values
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(cattleSlaughter_value_ha_fitted = predict(reg.cattleValueperHa, .))

# extract minimum positive fitted value
aux.min.positive.cattleSlaughter.value.ha.fitted <-
  muniTheta.prepData %>%
  dplyr::filter(cattleSlaughter_value_ha_fitted > 0) %>%
  dplyr::pull(cattleSlaughter_value_ha_fitted) %>%
  min()

# winsorize cattleSlaughter_value_ha_fitted
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(d_theta_winsorized = dplyr::if_else(cattleSlaughter_value_ha_fitted <= 0,
                                                    1,
                                                    0),
                cattleSlaughter_value_ha_fitted = dplyr::if_else(cattleSlaughter_value_ha_fitted <= 0,
                                                                 aux.min.positive.cattleSlaughter.value.ha.fitted,
                                                                 cattleSlaughter_value_ha_fitted))
# clean environment
rm(reg.cattleValueperHa, aux.min.positive.cattleSlaughter.value.ha.fitted)


# REGRESSION - CATTLE VALUE (2006)

# cattle value per ha
reg.cattleValueperHa <-
  muniTheta.prepData %>%
  lm(formula = cattleSlaughter2006_value_ha ~ pastureArea2006_value + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
       lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = pastureArea2006_value)

# regression results
summary(reg.cattleValueperHa)

# extract fitted values
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(cattleSlaughter2006_value_ha_fitted = predict(reg.cattleValueperHa, .))

# extract minimum positive fitted value
aux.min.positive.cattleSlaughter.value.ha.fitted <-
  muniTheta.prepData %>%
  dplyr::filter(cattleSlaughter2006_value_ha_fitted > 0) %>%
  dplyr::pull(cattleSlaughter2006_value_ha_fitted) %>%
  min()

# winsorize cattleSlaughter_value_ha_fitted
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(d_theta2006_winsorized = dplyr::if_else(cattleSlaughter2006_value_ha_fitted <= 0,
                                                        1,
                                                        0),
                cattleSlaughter2006_value_ha_fitted = dplyr::if_else(cattleSlaughter2006_value_ha_fitted <= 0,
                                                                     aux.min.positive.cattleSlaughter.value.ha.fitted,
                                                                     cattleSlaughter2006_value_ha_fitted))
# clean environment
rm(reg.cattleValueperHa, aux.min.positive.cattleSlaughter.value.ha.fitted)


# REGRESSION - FARM GATE PRICE (2017)

# cattle value per ha
reg.cattleFarmGatePricePerArroba <-
  muniTheta.prepData %>%
  lm(formula = cattleSlaughter_farmGatePrice ~ cattleSlaughter_head + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
       lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = cattleSlaughter_head)

# regression results
summary(reg.cattleFarmGatePricePerArroba)

# extract fitted values
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(cattleSlaughter_farmGatePrice_fitted = predict(reg.cattleFarmGatePricePerArroba, .))

# clean environment
rm(reg.cattleFarmGatePricePerArroba)



# REGRESSION - FARM GATE PRICE (2006)

# cattle value per ha
reg.cattleFarmGatePricePerArroba <-
  muniTheta.prepData %>%
  lm(formula = cattleSlaughter2006_farmGatePrice ~ cattleSlaughter2006_head + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
       lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = cattleSlaughter2006_head)

# regression results
summary(reg.cattleFarmGatePricePerArroba)

# extract fitted values
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(cattleSlaughter2006_farmGatePrice_fitted = predict(reg.cattleFarmGatePricePerArroba, .))

# clean environment
rm(reg.cattleFarmGatePricePerArroba)


# change crs to match with calibration.5SitesModel and select variables
muniTheta.prepData <-
  muniTheta.prepData %>%
  sf::st_transform(crs = sf::st_crs(calibration.5SitesModel)) %>%
  dplyr::select(muni_code, muni_area, cattleSlaughter_value_ha_fitted, cattleSlaughter_farmGatePrice_fitted, pastureArea_value, d_theta_winsorized,
                cattleSlaughter2006_value_ha_fitted, cattleSlaughter2006_farmGatePrice_fitted, pastureArea2006_value, d_theta2006_winsorized)

# match munis with sites
aux.theta <- sf::st_intersection(calibration.5SitesModel, muniTheta.prepData)

# calculate muni areas inside each site
aux.theta$muni_site_area <-
  sf::st_area(aux.theta) %>%
  units::set_units(ha) %>%
  unclass()

# calculate cattleSlaughter_value_ha_fitted and pastureArea_value by site (for each muni adjust the value by the share of the muni area inside the site)
aux.theta.2017 <-
  aux.theta %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(cattleSlaughter_value_ha_fitted = weighted.mean(cattleSlaughter_value_ha_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter_value_ha_scenario2 = weighted.mean((74.7/15)*cattleSlaughter_farmGatePrice_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter_value_ha_scenario3 = weighted.mean((140.2/15)*cattleSlaughter_farmGatePrice_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter_value_ha_scenario4 = weighted.mean((201.7/15)*cattleSlaughter_farmGatePrice_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter_value_ha_scenario5 = weighted.mean((221.4/15)*cattleSlaughter_farmGatePrice_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   pastureArea_value = sum(pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   d_theta_winsorized = min(d_theta_winsorized, na.rm = T))

# add cattleSlaughter_value_ha_fitted and pastureArea_value to spatial variables
calibration.5SitesModel <- left_join(calibration.5SitesModel, aux.theta.2017)

aux.theta.2006 <-
  aux.theta %>%
  sf::st_drop_geometry() %>%
  dplyr::filter(!is.na(pastureArea2006_value)) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(cattleSlaughter2006_value_ha_fitted = weighted.mean(cattleSlaughter2006_value_ha_fitted, w = pastureArea2006_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter2006_value_ha_scenario2 = weighted.mean((74.7/15)*cattleSlaughter2006_farmGatePrice_fitted, w = pastureArea2006_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter2006_value_ha_scenario3 = weighted.mean((140.2/15)*cattleSlaughter2006_farmGatePrice_fitted, w = pastureArea2006_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter2006_value_ha_scenario4 = weighted.mean((201.7/15)*cattleSlaughter2006_farmGatePrice_fitted, w = pastureArea2006_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter2006_value_ha_scenario5 = weighted.mean((221.4/15)*cattleSlaughter2006_farmGatePrice_fitted, w = pastureArea2006_value*(muni_site_area/muni_area), na.rm = T),
                   pastureArea2006_value = sum(pastureArea2006_value*(muni_site_area/muni_area), na.rm = T),
                   d_theta2006_winsorized = min(d_theta2006_winsorized, na.rm = T))

# add cattleSlaughter_value_ha_fitted and pastureArea_value to spatial variables
calibration.5SitesModel <- left_join(calibration.5SitesModel, aux.theta.2006)

# clean environment
rm(muniTheta.prepData, aux.theta)

# calculate theta_5Sites
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(theta_5Sites = cattleSlaughter_value_ha_fitted/(aux.price.2017),
                theta_scenario2_5Sites = cattleSlaughter_value_ha_scenario2/(aux.price.2017),
                theta_scenario3_5Sites = cattleSlaughter_value_ha_scenario3/(aux.price.2017),
                theta_scenario4_5Sites = cattleSlaughter_value_ha_scenario4/(aux.price.2017),
                theta_scenario5_5Sites = cattleSlaughter_value_ha_scenario5/(aux.price.2017)) %>%
  dplyr::mutate(theta2006_5Sites = cattleSlaughter2006_value_ha_fitted/(aux.price.2006),
                theta2006_scenario2_5Sites = cattleSlaughter2006_value_ha_scenario2/(aux.price.2006),
                theta2006_scenario3_5Sites = cattleSlaughter2006_value_ha_scenario3/(aux.price.2006),
                theta2006_scenario4_5Sites = cattleSlaughter2006_value_ha_scenario4/(aux.price.2006),
                theta2006_scenario5_5Sites = cattleSlaughter2006_value_ha_scenario5/(aux.price.2006)) %>%
  dplyr::select(-cattleSlaughter_value_ha_fitted, -cattleSlaughter_value_ha_scenario2, -cattleSlaughter_value_ha_scenario3,
                -cattleSlaughter_value_ha_scenario4, -cattleSlaughter_value_ha_scenario5,
                -cattleSlaughter2006_value_ha_fitted, -cattleSlaughter2006_value_ha_scenario2, -cattleSlaughter2006_value_ha_scenario3,
                -cattleSlaughter2006_value_ha_scenario4, -cattleSlaughter2006_value_ha_scenario5)


# identify adjacent neighbors
aux.neighbors <- sf::st_is_within_distance(calibration.5SitesModel, calibration.5SitesModel, dist = 100, remove_self = TRUE)

# impute values for missing thetas based on the average of adjacent neighbors
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(d_imputation_theta_scenario2_5Sites = dplyr::if_else(is.na(theta_scenario2_5Sites), 1, 0),
                theta_scenario2_5Sites = dplyr::if_else(is.na(theta_scenario2_5Sites),
                                                           apply(aux.neighbors, 1, function(i){mean(.$theta_scenario2_5Sites[i], na.rm = TRUE)}),
                                                           theta_scenario2_5Sites),
                d_imputation_theta_scenario3_5Sites = dplyr::if_else(is.na(theta_scenario3_5Sites), 1, 0),
                theta_scenario3_5Sites = dplyr::if_else(is.na(theta_scenario3_5Sites),
                                                           apply(aux.neighbors, 1, function(i){mean(.$theta_scenario3_5Sites[i], na.rm = TRUE)}),
                                                           theta_scenario3_5Sites),
                d_imputation_theta_scenario4_5Sites = dplyr::if_else(is.na(theta_scenario4_5Sites), 1, 0),
                theta_scenario4_5Sites = dplyr::if_else(is.na(theta_scenario4_5Sites),
                                                           apply(aux.neighbors, 1, function(i){mean(.$theta_scenario4_5Sites[i], na.rm = TRUE)}),
                                                           theta_scenario4_5Sites),
                d_imputation_theta_scenario5_5Sites = dplyr::if_else(is.na(theta_scenario5_5Sites), 1, 0),
                theta_scenario5_5Sites = dplyr::if_else(is.na(theta_scenario5_5Sites),
                                                           apply(aux.neighbors, 1, function(i){mean(.$theta_scenario5_5Sites[i], na.rm = TRUE)}),
                                                           theta_scenario5_5Sites))

# THETA (SOYBEAN)

# DATA INPUT
# load soybean potential yield raster
raster.soybean <- terra::rast(here::here("data/raw2clean/potentialYield_faogaez/output/clean_potentialYield.tif"))



# STORE PARAMETER VALUES
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::bind_cols(dplyr::tibble(theta_soybean_5Sites = terra::extract(raster.soybean, terra::vect(calibration.5SitesModel %>% sf::st_transform(4326)), mean, na.rm = T)$soyb200b_yld))


# clean environment
rm(raster.soybean)





# EMPLOYMENT LOSS ------------------------------------------------------------------------------------------------------------------------------------

# DATA INPUT
# load variables at the muni level to calibrate theta
load(here::here("data/projectSpecific/prepData/muniAgCensusCattleRaising_prepData.Rdata"))



# DATA MANIPULATION

# REGRESSION - CATTLE VALUE

# cattle value per ha
reg.workersperHa <-
  muniAgCensusCattleRaising.prepData %>%
  lm(formula = workers_ha ~ pastureArea_value + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
       lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = pastureArea_value)

# regression results
summary(reg.workersperHa)

# extract fitted values
muniAgCensusCattleRaising.prepData <-
  muniAgCensusCattleRaising.prepData %>%
  dplyr::mutate(workers_ha_fitted = predict(reg.workersperHa, .))

# extract minimum positive fitted value
aux.min.positive.workers.ha.fitted <-
  muniAgCensusCattleRaising.prepData %>%
  dplyr::filter(workers_ha_fitted > 0) %>%
  dplyr::pull(workers_ha_fitted) %>%
  min()

# windsorize workers_ha_fitted
muniAgCensusCattleRaising.prepData <-
  muniAgCensusCattleRaising.prepData %>%
  dplyr::mutate(workers_ha_fitted = dplyr::if_else(workers_ha_fitted <= 0,
                                                   aux.min.positive.workers.ha.fitted,
                                                   workers_ha_fitted))

# clean environment
rm(reg.workersperHa, aux.min.positive.workers.ha.fitted)


# change crs to match with calibration.5SitesModel and select variables
muniAgCensusCattleRaising.prepData <-
  muniAgCensusCattleRaising.prepData %>%
  sf::st_transform(crs = sf::st_crs(calibration.5SitesModel)) %>%
  dplyr::select(muni_code, muni_area, workers_ha_fitted, pastureArea_value)

# match munis with sites
aux.worker <- sf::st_intersection(calibration.5SitesModel, muniAgCensusCattleRaising.prepData)

# calculate muni areas inside each site
aux.worker$muni_site_area <-
  sf::st_area(aux.worker) %>%
  units::set_units(ha) %>%
  unclass()

# calculate workers_ha_fitted and pastureArea_value by site (for each muni adjust the value by the share of the muni area inside the site)
aux.worker <-
  aux.worker %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(workers_ha_fitted = weighted.mean(workers_ha_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T))

# add workers_ha_fitted and pastureArea_value to spatial variables
calibration.5SitesModel <- left_join(calibration.5SitesModel, aux.worker)

# clean environment
rm(muniAgCensusCattleRaising.prepData, aux.worker)

# calculate theta_5Sites
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(workerPerHa_5Sites = workers_ha_fitted) %>%
  dplyr::select(-workers_ha_fitted)




# PRICE INITIAL CONDITION ----------------------------------------------------------------------------------------------------------------------------

# estimates of p_2017 same as in the global model
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::mutate(p_2017_5Sites = calibration.globalModel$p_2017_global,
                p_2017_alt_5Sites = calibration.globalModel$p_2017_alt_global,
                pSoybean_2017_5Sites = calibration.globalModel$pSoybean_2017_global,
                pSoybean_2017_alt_5Sites = calibration.globalModel$pSoybean_2017_alt_global)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# REMOVE NAs
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::filter(!is.na(gamma_5Sites))


# ORDER VARIABLES
calibration.5SitesModel <-
  calibration.5SitesModel %>%
  dplyr::select(id, ends_with("ha_5Sites"),
                z_2017_5Sites, zbar_2017_5Sites, x_2017_5Sites, b_5Sites, theta_5Sites, theta2006_5Sites, gamma_5Sites, gammaSD_5Sites,
                d_theta_winsorized, d_theta2006_winsorized,
                pastureArea_value, pastureArea2006_value, ends_with("_5Sites"))


# POST-TREATMENT OVERVIEW
# summary(calibration.5SitesModel)
# View(calibration.5SitesModel)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(calibration.5SitesModel,
     file = here::here("data/projectSpecific/5SitesModel",
                       "calibration_5SitesModel.Rdata"))

# remove spatial feature
calibration.5SitesModel <- calibration.5SitesModel %>% sf::st_drop_geometry()

readr::write_csv(calibration.5SitesModel,
                 file = here::here("data/projectSpecific/5SitesModel", "calibration_5SitesModel.csv"))


# CLEAN TEMP DIR
terra::tmpFiles(current = T, remove = T)
gc()



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/5SitesModel")






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------