
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (2 SITES MODEL)
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
here::i_am("code/projectSpecific/2SitesModel/calibration_projectSpecific_2SitesModel.R", uuid = "0fef2915-0a76-4b7f-9a2c-c0b9830c058b")


# START TIME
tictoc::tic(msg = "calibration_projectSpecific_2SitesModel script", log = T)


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

# AGGREGATE FROM 1000 SITES TO 2 SITES
# transform shares to areas
raster.1000Sites$amazonBiomeArea_ha_2Sites <- raster.1000Sites$share_amazonBiome*raster.1000Sites$pixelArea_ha
raster.1000Sites$forestArea_1985_ha_2Sites <- raster.1000Sites$share_forest_1985*raster.1000Sites$pixelArea_ha
raster.1000Sites$agriculturalUseArea_1985_ha_2Sites <- raster.1000Sites$share_agriculturalUse_1985*raster.1000Sites$pixelArea_ha
raster.1000Sites$otherArea_1985_ha_2Sites <- raster.1000Sites$share_other_1985*raster.1000Sites$pixelArea_ha
raster.1000Sites$forestArea_2017_ha_2Sites <- raster.1000Sites$share_forest_2017*raster.1000Sites$pixelArea_ha
raster.1000Sites$agriculturalUseArea_2017_ha_2Sites <- raster.1000Sites$share_agriculturalUse_2017*raster.1000Sites$pixelArea_ha
raster.1000Sites$otherArea_2017_ha_2Sites <- raster.1000Sites$share_other_2017*raster.1000Sites$pixelArea_ha

# select area variables
raster.1000Sites <- terra::subset(raster.1000Sites, c(9:15))

# aggregate in two steps because using fact =12 generated a significant distortion while fact=4 followed by fact=3 was very close (in terms of total area compared to the 1000 sites)
raster.1000Sites.2 <- terra::aggregate(raster.1000Sites, fact = 5, fun = sum, na.rm = T)
raster.2Sites <- terra::aggregate(raster.1000Sites.2, fact = 5, fun = sum, na.rm = T)

# extract variables as polygons, transform to sf, and project data for faster spatial manipulation
calibration.2SitesModel <- terra::as.polygons(raster.2Sites, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)

# add id variable
calibration.2SitesModel$id <- 1:length(calibration.2SitesModel$geometry)

# transform share aggregate in area (ha)
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(zbar_1985_2Sites = agriculturalUseArea_1985_ha_2Sites + forestArea_1985_ha_2Sites,
                zbar_2017_2Sites = agriculturalUseArea_2017_ha_2Sites + forestArea_2017_ha_2Sites) %>%
  dplyr::select(id, amazonBiomeArea_ha_2Sites,
                forestArea_1985_ha_2Sites, otherArea_1985_ha_2Sites,
                z_1985_2Sites = agriculturalUseArea_1985_ha_2Sites, zbar_1985_2Sites,
                forestArea_2017_ha_2Sites, otherArea_2017_ha_2Sites,
                z_2017_2Sites = agriculturalUseArea_2017_ha_2Sites, zbar_2017_2Sites)

# combine sites 2,3, and 4 to have a similar zbar to site 1
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(id = dplyr::if_else(id %in% c(1,4), 1, 2)) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(across(where(is.numeric), .fns = sum))





# PARAMETER GAMMA ------------------------------------------------------------------------------------------------------------------------------------

# # DATA INPUT
# # load pixel sample with biomass data
# load(here::here("data/projectSpecific/prepData/pixelBiomass_prepData.Rdata"))
#
# # select minicells of primary forest with co2 information and transform to spatial points
# pixelBiomass.prepData <-
#   pixelBiomass.prepData %>%
#   dplyr::filter(mapbiomas_classAgg == "primaryForest") %>%
#   sf::st_transform(sf::st_crs(calibration.2SitesModel)) %>%
#   dplyr::mutate(co2e_ha = (agb/2)*(44/12)) %>%
#   dplyr::select(co2e_ha)
#
# # match minicells with sites
# aux.gamma <- sf::st_join(calibration.2SitesModel, pixelBiomass.prepData)
#
# # calculate average carbon density on primary forest areas by site
# aux.gamma <-
#   aux.gamma %>%
#   sf::st_drop_geometry() %>%
#   dplyr::group_by(id) %>%
#   summarise(gamma_2Sites = mean(co2e_ha, na.rm = T))
#
#
# # add gamma_2Sites to spatial variables
# calibration.2SitesModel <- left_join(calibration.2SitesModel, aux.gamma)
#
# # clean environment
# rm(pixelBiomass.prepData, aux.gamma, raster.variables)
#
# # add alternative value for gamma with a 10.1% correction factor for small trees and lianas (see Malhi et al (2006))
# calibration.2SitesModel <-
#   calibration.2SitesModel %>%
#   dplyr::mutate(gamma_alt_2Sites = gamma_2Sites*1.101)



# DATA INPUT (2017)
# load pixel sample with biomass data
load(here::here("data/projectSpecific/prepData/pixelBiomass2017_prepData.Rdata"))

# select minicells of primary forest with co2 information and transform to spatial points
pixelBiomass2017.prepData <-
  pixelBiomass2017.prepData %>%
  dplyr::filter(mapbiomas_classAgg == "primaryForest") %>%
  sf::st_transform(sf::st_crs(calibration.2SitesModel)) %>%
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
aux.gamma2017 <- sf::st_join(calibration.2SitesModel, pixelBiomass2017.prepData)

# calculate average carbon density on primary forest areas by site
aux.gamma2017 <- aux.gamma2017 %>% sf::st_drop_geometry() %>% dplyr::group_by(id) %>% summarise(gamma_2Sites = mean(co2e_ha_2017, na.rm = T),
                                                                                                co2e_ha_2017 = mean(co2e_ha_2017, na.rm = T),
                                                                                                gammaSD_2Sites = mean(co2eSD_ha_2017_fitted, na.rm = T),
                                                                                                co2eSD_ha_2017 = mean(co2eSD_ha_2017, na.rm = T),
                                                                                                lon = mean(lon, na.rm = T),
                                                                                                lat = mean(lat, na.rm = T))

# add gamma_2Sites to spatial variables
calibration.2SitesModel <- left_join(calibration.2SitesModel, aux.gamma2017)

# clean environment
rm(pixelBiomass2017.prepData, aux.gamma2017)

calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(gammaSD_2Sites = predict(reg.carbonUncertainty, .),
                gammaSD_2Sites = exp(gammaSD_2Sites))

# clean environment
rm(reg.carbonUncertainty)







# PARAMETERS A AND B ---------------------------------------------------------------------------------------------------------------------------------

# estimate of a same as in the global model >
# estimate of b given a such that b_i/a is the average carbon in primary forest
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(a_2Sites = calibration.globalModel$a_global,
                b_2Sites = a_2Sites*gamma_2Sites*zbar_2017_2Sites)





# PARAMETER K ----------------------------------------------------------------------------------------------------------------------------------------

# estimate of k same as in the global model
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(k_2Sites = calibration.globalModel$k_global)





# PARAMETER ZETA ----------------------------------------------------------------------------------------------------------------------------------------

# estimate of zeta same as in the global model
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(zeta_2Sites = calibration.globalModel$zeta_global,
                zeta_alt_2Sites = calibration.globalModel$zeta_alt_global)




# PARAMETER RHO ----------------------------------------------------------------------------------------------------------------------------------------

# estimate of rho same as in the global model
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(rho_2Sites = calibration.globalModel$rho_global)





# PARAMETER R ----------------------------------------------------------------------------------------------------------------------------------------

# estimate of r same as in the global model
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(r_2Sites = calibration.globalModel$r_global)





# INITIAL CONDITIONS X -------------------------------------------------------------------------------------------------------------------------------

# x_2017_2Sites estimated as in the old way of global model, just considering the stock of carbon stored in forest areas assuming that all forests are primary
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(x_2017_2Sites = gamma_2Sites*forestArea_2017_ha_2Sites)





# INITIAL CONDITIONS C -------------------------------------------------------------------------------------------------------------------------------

# estimates of C same as in the global model
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(C_2017_2Sites = calibration.globalModel$C_2017_global,
                C0bern_2017_2Sites = calibration.globalModel$C0bern_2017_global,
                C1bern_2017_2Sites = calibration.globalModel$C1bern_2017_global,
                C2bern_2017_2Sites = calibration.globalModel$C2bern_2017_global,
                C3bern_2017_2Sites = calibration.globalModel$C3bern_2017_global)




# PARAMETER U ----------------------------------------------------------------------------------------------------------------------------------------

# estimates of u same as in the global model
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(u_2Sites = calibration.globalModel$u_global,
                u_bern_2Sites = calibration.globalModel$u_bern_global)





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
  dplyr::mutate(cattleSlaughter_value_ha_fitted = predict(reg.cattleValueperHa, .),
                cattleSlaughter_value_ha_resid = cattleSlaughter_value_ha-cattleSlaughter_value_ha_fitted)

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


# change crs to match with calibration.2SitesModel and select variables
muniTheta.prepData <-
  muniTheta.prepData %>%
  sf::st_transform(crs = sf::st_crs(calibration.2SitesModel)) %>%
  dplyr::select(muni_code, muni_area, cattleSlaughter_value_ha_fitted, cattleSlaughter_farmGatePrice_fitted, cattleSlaughter_value_ha_resid, pastureArea_value, d_theta_winsorized)

# match munis with sites
aux.theta <- sf::st_intersection(calibration.2SitesModel, muniTheta.prepData)

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
                   cattleSlaughter_value_ha_resid = weighted.mean(cattleSlaughter_value_ha_resid, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter_value_ha_scenario2 = weighted.mean((74.7/15)*cattleSlaughter_farmGatePrice_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter_value_ha_scenario3 = weighted.mean((140.2/15)*cattleSlaughter_farmGatePrice_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter_value_ha_scenario4 = weighted.mean((201.7/15)*cattleSlaughter_farmGatePrice_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter_value_ha_scenario5 = weighted.mean((221.4/15)*cattleSlaughter_farmGatePrice_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   pastureArea_value = sum(pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   d_theta_winsorized = min(d_theta_winsorized, na.rm = T))

# add cattleSlaughter_value_ha_fitted and pastureArea_value to spatial variables
calibration.2SitesModel <- left_join(calibration.2SitesModel, aux.theta.2017)

# clean environment
rm(muniTheta.prepData, aux.theta)

# calculate theta_2Sites
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(theta_2Sites = cattleSlaughter_value_ha_fitted/(aux.price.2017),
                theta_resid_2Sites = cattleSlaughter_value_ha_resid/(aux.price.2017),
                theta_scenario2_2Sites = cattleSlaughter_value_ha_scenario2/(aux.price.2017),
                theta_scenario3_2Sites = cattleSlaughter_value_ha_scenario3/(aux.price.2017),
                theta_scenario4_2Sites = cattleSlaughter_value_ha_scenario4/(aux.price.2017),
                theta_scenario5_2Sites = cattleSlaughter_value_ha_scenario5/(aux.price.2017)) %>%
  dplyr::select(-cattleSlaughter_value_ha_fitted, -cattleSlaughter_value_ha_scenario2, -cattleSlaughter_value_ha_scenario3,
                -cattleSlaughter_value_ha_scenario4, -cattleSlaughter_value_ha_scenario5)


# identify adjacent neighbors
aux.neighbors <- sf::st_is_within_distance(calibration.2SitesModel, calibration.2SitesModel, dist = 100, remove_self = TRUE)

# impute values for missing thetas based on the average of adjacent neighbors
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(d_imputation_theta_scenario2_2Sites = dplyr::if_else(is.na(theta_scenario2_2Sites), 1, 0),
                theta_scenario2_2Sites = dplyr::if_else(is.na(theta_scenario2_2Sites),
                                                           apply(aux.neighbors, 1, function(i){mean(.$theta_scenario2_2Sites[i], na.rm = TRUE)}),
                                                           theta_scenario2_2Sites),
                d_imputation_theta_scenario3_2Sites = dplyr::if_else(is.na(theta_scenario3_2Sites), 1, 0),
                theta_scenario3_2Sites = dplyr::if_else(is.na(theta_scenario3_2Sites),
                                                           apply(aux.neighbors, 1, function(i){mean(.$theta_scenario3_2Sites[i], na.rm = TRUE)}),
                                                           theta_scenario3_2Sites),
                d_imputation_theta_scenario4_2Sites = dplyr::if_else(is.na(theta_scenario4_2Sites), 1, 0),
                theta_scenario4_2Sites = dplyr::if_else(is.na(theta_scenario4_2Sites),
                                                           apply(aux.neighbors, 1, function(i){mean(.$theta_scenario4_2Sites[i], na.rm = TRUE)}),
                                                           theta_scenario4_2Sites),
                d_imputation_theta_scenario5_2Sites = dplyr::if_else(is.na(theta_scenario5_2Sites), 1, 0),
                theta_scenario5_2Sites = dplyr::if_else(is.na(theta_scenario5_2Sites),
                                                           apply(aux.neighbors, 1, function(i){mean(.$theta_scenario5_2Sites[i], na.rm = TRUE)}),
                                                           theta_scenario5_2Sites))

# THETA (SOYBEAN)

# DATA INPUT
# load soybean potential yield raster
raster.soybean <- terra::rast(here::here("data/raw2clean/potentialYield_faogaez/output/clean_potentialYield.tif"))



# STORE PARAMETER VALUES
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::bind_cols(dplyr::tibble(theta_soybean_2Sites = terra::extract(raster.soybean, terra::vect(calibration.2SitesModel %>% sf::st_transform(4326)), mean, na.rm = T)$soyb200b_yld))


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


# change crs to match with calibration.2SitesModel and select variables
muniAgCensusCattleRaising.prepData <-
  muniAgCensusCattleRaising.prepData %>%
  sf::st_transform(crs = sf::st_crs(calibration.2SitesModel)) %>%
  dplyr::select(muni_code, muni_area, workers_ha_fitted, pastureArea_value)

# match munis with sites
aux.worker <- sf::st_intersection(calibration.2SitesModel, muniAgCensusCattleRaising.prepData)

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
calibration.2SitesModel <- left_join(calibration.2SitesModel, aux.worker)

# clean environment
rm(muniAgCensusCattleRaising.prepData, aux.worker)

# calculate theta_2Sites
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(workerPerHa_2Sites = workers_ha_fitted) %>%
  dplyr::select(-workers_ha_fitted)




# PRICE INITIAL CONDITION ----------------------------------------------------------------------------------------------------------------------------

# estimates of p_2017 same as in the global model
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::mutate(p_2017_2Sites = calibration.globalModel$p_2017_global,
                p_2017_alt_2Sites = calibration.globalModel$p_2017_alt_global,
                pSoybean_2017_2Sites = calibration.globalModel$pSoybean_2017_global,
                pSoybean_2017_alt_2Sites = calibration.globalModel$pSoybean_2017_alt_global)




# ALTERNATIVE CALCULATION Z, X, C SERIES -------------------------------------------------------------------------------------------------------------

# clean environment
rm(list = ls(pattern = "aux"))
rm(list = ls(pattern = "raster"))
rm(sampleMuniSpatial.prepData, seriesPriceCattle.prepData, calibration.globalModel)

# load pixel sample with mapbiomas categories data
load(here::here("data/projectSpecific/prepData/pixelCategories_prepData.Rdata"))

# calculate z series
aux.zSeries <-
  pixelCategories.prepData %>%
  dplyr::slice_sample(prop = 0.01) %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  sf::st_join(sf::st_transform(calibration.2SitesModel, sf::st_crs(4326))) %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(z = dplyr::if_else(mapbiomas_classAgg %in% c("pasture", "otherCrops","soybean"), 1, 0)) %>% # create z variable (dummy =1 if agricultural use, =0 if forest)
  dplyr::select(id, year, z, amazonBiomeArea_ha_2Sites) %>% # select variables of interest
  dplyr::group_by(id, year) %>%
  dplyr::summarise(z = (sum(z)/n())*mean(amazonBiomeArea_ha_2Sites))

rm(pixelCategories.prepData)


aux.zSeries.site1 <- aux.zSeries %>% dplyr::filter(id == 1)
aux.zSeries.site2 <- aux.zSeries %>% dplyr::filter(id == 2)

aux.zSeries.site1 <-
  aux.zSeries.site1 %>%
  dplyr::bind_rows(tibble(id = 1,year = 1975:1984, z = NA)) %>%
  dplyr::arrange(year) %>%
  dplyr::mutate(z_predicted = predict.lm(object = lm(data = ., z ~ year),
                                         newdata = data.frame(year = 1975:2019))) %>%
  dplyr::mutate(z_mixed = dplyr::case_when(year == 1975 ~ 0,
                                           year %in% 1976:1984 ~ z_predicted,
                                           TRUE ~ z)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(z_dot = z_mixed - lag(z_mixed)) # calculate forest annual changes (deforestation if positive, regeneration if negative)

aux.zSeries.site2 <-
  aux.zSeries.site2 %>%
  dplyr::bind_rows(tibble(id = 2, year = 1971:1984, z = NA)) %>%
  dplyr::arrange(year) %>%
  dplyr::mutate(z_predicted = predict.lm(object = lm(data = ., z ~ year),
                                         newdata = data.frame(year = 1971:2019))) %>%
  dplyr::mutate(z_mixed = dplyr::case_when(year == 1971 ~ 0,
                                           year %in% 1972:1984 ~ z_predicted,
                                           TRUE ~ z)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(z_dot = z_mixed - lag(z_mixed)) # calculate forest annual changes (deforestation if positive, regeneration if negative)


# generate variables to store values starting at b/a and 0
aux.zSeries.site1$x_t <- calibration.2SitesModel %>% dplyr::filter(id == 1) %>% dplyr::mutate(x_t = a_2Sites*gamma_2Sites*zbar_2017_2Sites) %>% pull(x_t)/calibration.2SitesModel[1,]$a_2Sites
aux.zSeries.site1$e_t <- 0

# calculate e_t and x_t
for (y in seq_along(1975:2017)) {

  aux.zSeries.site1$x_t[y+1] <- aux.zSeries.site1$x_t[y] - aux.zSeries.site1$z_dot[y+1]*calibration.2SitesModel[1,]$gamma_2Sites - aux.zSeries.site1$x_t[y]*calibration.2SitesModel[1,]$a_2Sites + calibration.2SitesModel[1,]$b_2Sites*(1-aux.zSeries.site1$z_mixed[y+1]/calibration.2SitesModel[1,]$zbar_2017_2Sites)
  aux.zSeries.site1$e_t[y+1] <- -(aux.zSeries.site1$x_t[y+1]-aux.zSeries.site1$x_t[y]) + aux.zSeries.site1$z_mixed[y+1]*calibration.2SitesModel[1,]$k_2Sites

}

aux.zSeries.site2$x_t <- calibration.2SitesModel %>% dplyr::filter(id == 2) %>% dplyr::mutate(x_t = a_2Sites*gamma_2Sites*zbar_2017_2Sites) %>% pull(x_t)/calibration.2SitesModel[2,]$a_2Sites
aux.zSeries.site2$e_t <- 0

# calculate e_t and x_t
for (y in seq_along(1971:2017)) {

  aux.zSeries.site2$x_t[y+1] <- aux.zSeries.site2$x_t[y] - aux.zSeries.site2$z_dot[y+1]*calibration.2SitesModel[2,]$gamma_2Sites - aux.zSeries.site2$x_t[y]*calibration.2SitesModel[2,]$a_2Sites + calibration.2SitesModel[2,]$b_2Sites*(1-aux.zSeries.site2$z_mixed[y+1]/calibration.2SitesModel[2,]$zbar_2017_2Sites)
  aux.zSeries.site2$e_t[y+1] <- -(aux.zSeries.site2$x_t[y+1]-aux.zSeries.site2$x_t[y]) + aux.zSeries.site2$z_mixed[y+1]*calibration.2SitesModel[2,]$k_2Sites

}

# SIMPLIFIED CARBON ACCUMULATION
aux.zSeries.site1$C_t <- 0

for (y in seq_along(1975:2017)) {

  aux.zSeries.site1$C_t[y+1] <-  aux.zSeries.site1$e_t[y+1] +aux.zSeries.site1$C_t[y] -0.025*aux.zSeries.site1$C_t[y]

}

aux.zSeries.site2$C_t <- 0

for (y in seq_along(1971:2017)) {

  aux.zSeries.site2$C_t[y+1] <-  aux.zSeries.site2$e_t[y+1] +aux.zSeries.site2$C_t[y] -0.025*aux.zSeries.site2$C_t[y]

}

# EXPORT TABLE WITH Z, X AND C VALUES FOR DIFFERENT YEARS (1973-2017)
variablesEvolution.2SitesModel <-
  aux.zSeries.site1 %>%
  dplyr::bind_rows(aux.zSeries.site2) %>%
  dplyr::filter(year <= 2017) %>%
  dplyr::select(id, year, `z (ha)` = z_mixed, `x (CO2e Mg)` = x_t, `C (CO2e Mg)` = C_t) %>%
  dplyr::mutate(data_type = dplyr::if_else(year < 1985, "linear extrapolation", "observed"))

readr::write_csv(variablesEvolution.2SitesModel,
                 file = here::here("data/projectSpecific/2SitesModel", "variablesEvolution_2SitesModel.csv"))




# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# REMOVE NAs
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::filter(!is.na(gamma_2Sites))


# ORDER VARIABLES
calibration.2SitesModel <-
  calibration.2SitesModel %>%
  dplyr::select(id, ends_with("ha_2Sites"),
                z_2017_2Sites, zbar_2017_2Sites, x_2017_2Sites, b_2Sites, theta_2Sites, gamma_2Sites, gammaSD_2Sites,
                d_theta_winsorized, theta_resid_2Sites,
                pastureArea_value, ends_with("_2Sites"))


# POST-TREATMENT OVERVIEW
# summary(calibration.2SitesModel)
# View(calibration.2SitesModel)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(calibration.2SitesModel,
     file = here::here("data/projectSpecific/2SitesModel",
                       "calibration_2SitesModel.Rdata"))

# remove spatial feature
calibration.2SitesModel <- calibration.2SitesModel %>% sf::st_drop_geometry()

readr::write_csv(calibration.2SitesModel,
                 file = here::here("data/projectSpecific/2SitesModel", "calibration_2SitesModel.csv"))


# CLEAN TEMP DIR
terra::tmpFiles(current = T, remove = T)
gc()



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/2SitesModel")






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------