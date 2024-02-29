
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (1059 SITES MODEL)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("code/setup.R")


# START TIMER
tictoc::tic(msg = "calibration_1059SitesModel.R script", log = T)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster.variables <- terra::rast(list.files(here::here("data/calibration/1055SitesModel/aux_tifs"),
                                           pattern = "raster_",
                                           full.names = T))


# MUNI LEVEL SPATIAL SAMPLE
load(here::here("data/calibration/prepData/sampleMuniSpatial_prepData.Rdata"))





# INITIAL CONDITIONS Z -------------------------------------------------------------------------------------------------------------------------------

# MAPBIOMAS VARIABLES + AMAZON BIOME + PIXEL AREA (Z_1059Sites CONSTRUCTION)
# extract variables as polygons, transform to sf, and project data for faster spatial manipulation
calibration.1059SitesModel <- terra::as.polygons(raster.variables, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)

# remove sites with less than 1% of its are intersecting with the amazon biome
calibration.1059SitesModel <-
  calibration.1059SitesModel %>%
  dplyr::filter(share_amazonBiome >= 0.001)

# add id variable
calibration.1059SitesModel$id <- 1:nrow(calibration.1059SitesModel)

# transform share variables in area (ha)
calibration.1059SitesModel <-
  calibration.1059SitesModel %>%
  dplyr::mutate(amazonBiomeArea_ha_1059Sites = share_amazonBiome*pixelArea_ha,
                forestArea_2017_ha_1059Sites = share_forest_2017*pixelArea_ha,
                z_2017_1059Sites = share_agriculturalUse_2017*pixelArea_ha,
                otherArea_2017_ha_1059Sites = share_other_2017*pixelArea_ha,
                zbar_2017_1059Sites = forestArea_2017_ha_1059Sites + z_2017_1059Sites,
                forestArea_1995_ha_1059Sites = share_forest_1995*pixelArea_ha,
                z_1995_1059Sites = share_agriculturalUse_1995*pixelArea_ha,
                otherArea_1995_ha_1059Sites = share_other_1995*pixelArea_ha,
                zbar_1995_1059Sites = forestArea_1995_ha_1059Sites + z_1995_1059Sites) %>%
  dplyr::select(id, siteArea_ha_1059Sites = pixelArea_ha, amazonBiomeArea_ha_1059Sites,
                forestArea_2017_ha_1059Sites, otherArea_2017_ha_1059Sites, z_2017_1059Sites, zbar_2017_1059Sites,
                forestArea_1995_ha_1059Sites, otherArea_1995_ha_1059Sites, z_1995_1059Sites, zbar_1995_1059Sites)



# PARAMETER GAMMA ------------------------------------------------------------------------------------------------------------------------------------

# DATA INPUT (2010)
# load pixel sample with biomass data
load(here::here("data/calibration/prepData/pixelBiomass2010_prepData.Rdata"))

# select minicells of primary forest with co2 information and transform to spatial points
pixelBiomass2010.prepData <-
  pixelBiomass2010.prepData %>%
  sf::st_transform(sf::st_crs(calibration.1059SitesModel)) %>%
  dplyr::mutate(co2e_ha_2010 = (agb_2010/2)*(44/12)) %>%
  dplyr::filter(co2e_ha_2010 > 0, !is.na(co2e_ha_2010))

# match minicells with sites
site.gamma2010 <-
  sf::st_join(calibration.1059SitesModel %>% dplyr::select(id),
              pixelBiomass2010.prepData %>% dplyr::select(co2e_ha_2010)) %>%
  sf::st_drop_geometry()

# calculate average carbon density on primary forest areas by site
aux.gamma2010 <-
  site.gamma2010 %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(gamma2010_1059Sites = mean(co2e_ha_2010, na.rm = T))

# add gamma_1059Sites to spatial variables
calibration.1059SitesModel <- dplyr::left_join(calibration.1059SitesModel, aux.gamma2010)


# clean environment
rm(pixelBiomass2010.prepData, aux.gamma2010)


# DATA INPUT (2017)
# load pixel sample with biomass data
load(here::here("data/calibration/prepData/pixelBiomass2017_prepData.Rdata"))

# select minicells of primary forest with co2 information and transform to spatial points
pixelBiomass2017.prepData <-
  pixelBiomass2017.prepData %>%
  sf::st_transform(sf::st_crs(calibration.1059SitesModel)) %>%
  dplyr::mutate(co2e_ha_2017 = (agb_2017/2)*(44/12)) %>%
  dplyr::filter(co2e_ha_2017 > 0, !is.na(co2e_ha_2017))

# match minicells with sites
site.gamma2017 <-
  sf::st_join(calibration.1059SitesModel %>% dplyr::select(id),
              pixelBiomass2017.prepData %>% dplyr::select(co2e_ha_2017)) %>%
  sf::st_drop_geometry()

# calculate average carbon density on primary forest areas by site
aux.gamma2017 <-
  site.gamma2017 %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(gamma2017_1059Sites = mean(co2e_ha_2017, na.rm = T))

# add gamma_1059Sites to spatial variables
calibration.1059SitesModel <- dplyr::left_join(calibration.1059SitesModel, aux.gamma2017)


# clean environment
rm(pixelBiomass2017.prepData, aux.gamma2017)


# DATA INPUT (2018)
# load pixel sample with biomass data
load(here::here("data/calibration/prepData/pixelBiomass2018_prepData.Rdata"))

# select minicells of primary forest with co2 information and transform to spatial points
pixelBiomass2018.prepData <-
  pixelBiomass2018.prepData %>%
  sf::st_transform(sf::st_crs(calibration.1059SitesModel)) %>%
  dplyr::mutate(co2e_ha_2018 = (agb_2018/2)*(44/12)) %>%
  dplyr::filter(co2e_ha_2018 > 0, !is.na(co2e_ha_2018))

# match minicells with sites
site.gamma2018 <-
  sf::st_join(calibration.1059SitesModel %>% dplyr::select(id),
              pixelBiomass2018.prepData %>% dplyr::select(co2e_ha_2018)) %>%
  sf::st_drop_geometry()

# calculate average carbon density on primary forest areas by site
aux.gamma2018 <-
  site.gamma2018 %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(gamma2018_1059Sites = mean(co2e_ha_2018, na.rm = T))

# add gamma_1059Sites to spatial variables
calibration.1059SitesModel <- dplyr::left_join(calibration.1059SitesModel, aux.gamma2018)


# clean environment
rm(pixelBiomass2018.prepData, aux.gamma2018)


# identify adjacent neighbors
aux.neighbors <- sf::st_is_within_distance(calibration.1059SitesModel, calibration.1059SitesModel, dist = 100, remove_self = TRUE)

# impute values for missing gammas based on the average of adjacent neighbors
calibration.1059SitesModel <-
  calibration.1059SitesModel %>%
  dplyr::mutate(gamma2010_1059Sites = dplyr::if_else(is.na(gamma2010_1059Sites),
                                                     apply(aux.neighbors, 1, function(i){mean(.$gamma2010_1059Sites[i], na.rm = TRUE)}),
                                                     gamma2010_1059Sites),
                gamma2017_1059Sites = dplyr::if_else(is.na(gamma2017_1059Sites),
                                                     apply(aux.neighbors, 1, function(i){mean(.$gamma2017_1059Sites[i], na.rm = TRUE)}),
                                                     gamma2017_1059Sites),
                gamma2018_1059Sites = dplyr::if_else(is.na(gamma2018_1059Sites),
                                                     apply(aux.neighbors, 1, function(i){mean(.$gamma2018_1059Sites[i], na.rm = TRUE)}),
                                                     gamma2018_1059Sites))


# set baseline gamma and gammaSD as the data from 2017 and calculate alternative  gamma gammaSD based on the mean and sd of gamma2010,gamma2017, and gamma2018
calibration.1059SitesModel <-
  calibration.1059SitesModel %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(gamma_1059Sites = rowMeans(across(c("gamma2010_1059Sites", "gamma2017_1059Sites", "gamma2018_1059Sites"))),
                gammaSD_1059Sites = apply(across(c("gamma2010_1059Sites", "gamma2017_1059Sites", "gamma2018_1059Sites")), 1, sd)) %>%
  dplyr::ungroup()





# PARAMETER ALPHA ------------------------------------------------------------------------------------------------------------------------------------

# estimate of alpha same as in the global model
calibration.1059SitesModel <-
  calibration.1059SitesModel %>%
  dplyr::mutate(alpha_1059Sites = 1 - (1-0.99)^(1/100))





# PARAMETER KAPPA ----------------------------------------------------------------------------------------------------------------------------------------

# DATA INPUT
# load pixel sample with biomass data
load(here::here("data/calibration/prepData/stateEmission_prepData.Rdata"))


# DATA MANIPULATION
# calculate average of net emission factor from agricultural use across years and states
avg.netEmissionFactor <-
  stateEmission.prepData %>%
  dplyr::ungroup() %>%
  dplyr::summarise(netEmissionFactor_co2e_ha = sum(netEmission_co2e)/sum(agriculturalUse_area)) %>%
  pull(netEmissionFactor_co2e_ha)


# estimate of kappa same as in the global model
calibration.1059SitesModel <-
  calibration.1059SitesModel %>%
  dplyr::mutate(kappa_1059Sites = avg.netEmissionFactor)

# clean environment
rm(avg.netEmissionFactor)





# PARAMETER ZETA ----------------------------------------------------------------------------------------------------------------------------------------

# zeta is calibrated such that the marginal cost of changing land use (zeta*forestToPastureTransitionArea) matches the forest to pasture transition cost >
# estimated by Araujo, Costa and Sant'Anna (2022), reported on column 4 of table 4 (on the right). We transform to dollars using an FX rate of 4.14 (December 2019) >
# as in the paper: 1614.54/4.14 = 390 USD/ha. The paper also calculates that 6.5% of the forest area was converted to pasture between 2008 and 2017. >
aux.transitionCost <- 1614.54/4.14

# The forest area in 2008 represented 72% of the Legal Amazon area, which covers 501,506,775 ha, so the transition in hectares is >
#  0.065*0.72*501,506,775 = 23,470,517, resulting in an annual average of 2,347,052 ha.
aux.transitionArea <- (0.065*0.72*501506775)/(2017-2008+1)

zeta <- aux.transitionCost/aux.transitionArea
zeta_alt <- 483/aux.transitionArea # Alternative value based on a quote from (https://www.otempo.com.br/brasil/investigacoes-revelam-quadrilhas-e-ganho-milionario-por-tras-do-desmate-1.2229571)

# estimate of zeta same as in the global model
calibration.1059SitesModel <-
  calibration.1059SitesModel %>%
  dplyr::mutate(zeta_1059Sites = zeta,
                zeta_alt_1059Sites = zeta_alt)





# INITIAL CONDITIONS X -------------------------------------------------------------------------------------------------------------------------------

# x_2017_1059Sites estimated as in the old way of global model, just considering the stock of carbon stored in forest areas assuming that all forests are primary
calibration.1059SitesModel <-
  calibration.1059SitesModel %>%
  dplyr::mutate(x_2017_1059Sites = gamma_1059Sites*(zbar_2017_1059Sites-z_2017_1059Sites))





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
reg.cattleValueperHa.2017 <-
  muniTheta.prepData %>%
  lm(formula = cattleSlaughter_valuePerHa_2017 ~ pasture_area_2017 + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
       lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = pasture_area_2017)

# regression results
summary(reg.cattleValueperHa.2017)

# extract fitted values
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(cattleSlaughter_valuePerHa_fitted_2017 = predict(reg.cattleValueperHa.2017, .))

# extract minimum positive fitted value
aux.min.positive.cattleSlaughter.value.ha.fitted.2017 <-
  muniTheta.prepData %>%
  dplyr::filter(cattleSlaughter_valuePerHa_fitted_2017 > 0) %>%
  dplyr::pull(cattleSlaughter_valuePerHa_fitted_2017) %>%
  min()

# winsorize cattleSlaughter_valuePerHa_fitted
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(d_theta_winsorized_2017 = dplyr::if_else(cattleSlaughter_valuePerHa_fitted_2017 <= 0,
                                                         1,
                                                         0),
                cattleSlaughter_valuePerHa_fitted_2017 = dplyr::if_else(cattleSlaughter_valuePerHa_fitted_2017 <= 0,
                                                                        aux.min.positive.cattleSlaughter.value.ha.fitted.2017,
                                                                        cattleSlaughter_valuePerHa_fitted_2017))
# clean environment
rm(reg.cattleValueperHa.2017, aux.min.positive.cattleSlaughter.value.ha.fitted.2017)

# match munis with sites
site.theta.2017 <- sf::st_intersection(calibration.1059SitesModel %>% dplyr::select(id),
                                       muniTheta.prepData       %>% dplyr::select(muni_code, muni_area, cattleSlaughter_valuePerHa_fitted_2017,
                                                                                  pasture_area_2017, d_theta_winsorized_2017))

# calculate muni areas inside each site
site.theta.2017$muni_site_area <-
  sf::st_area(site.theta.2017) %>%
  units::set_units(ha) %>%
  unclass()

# drop spatial feature
site.theta.2017 <-
  site.theta.2017 %>%
  sf::st_drop_geometry()

# calculate cattleSlaughter_valuePerHa_fitted and pastureArea_value by site (for each muni adjust the value by the share of the muni area inside the site)
aux.theta.2017 <-
  site.theta.2017 %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(theta2017_1059Sites = weighted.mean(cattleSlaughter_valuePerHa_fitted_2017/aux.price.2017, w = pasture_area_2017*(muni_site_area/muni_area), na.rm = T),
                   pasture_area_2017 = sum(pasture_area_2017*(muni_site_area/muni_area), na.rm = T),
                   d_theta_winsorized_2017 = min(d_theta_winsorized_2017, na.rm = T))

# add cattleSlaughter_valuePerHa_fitted and pastureArea_value to spatial variables
calibration.1059SitesModel <- dplyr::left_join(calibration.1059SitesModel, aux.theta.2017)

# clean environment
rm(aux.theta.2017)



# EXTRACT AVERAGE 2006 PRICE (use real prices because it is normalized to 2017 )
aux.price.2006 <-
  seriesPriceCattle.prepData %>%
  dplyr::filter(year == 2006) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(mean_price_2006 = mean(price_real_mon_cattle)/3.192) %>% # BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
  dplyr::pull(mean_price_2006)


# REGRESSION - CATTLE VALUE (2006)

# cattle value per ha
reg.cattleValueperHa.2006 <-
  muniTheta.prepData %>%
  lm(formula = cattleSlaughter_valuePerHa_2006 ~ pasture_area_2006 + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
       lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = pasture_area_2006)

# regression results
summary(reg.cattleValueperHa.2006)

# extract fitted values
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(cattleSlaughter_valuePerHa_fitted_2006 = predict(reg.cattleValueperHa.2006, .))

# extract minimum positive fitted value
aux.min.positive.cattleSlaughter.value.ha.fitted.2006 <-
  muniTheta.prepData %>%
  dplyr::filter(cattleSlaughter_valuePerHa_fitted_2006 > 0) %>%
  dplyr::pull(cattleSlaughter_valuePerHa_fitted_2006) %>%
  min()

# winsorize cattleSlaughter_valuePerHa_fitted
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(d_theta_winsorized_2006 = dplyr::if_else(cattleSlaughter_valuePerHa_fitted_2006 <= 0,
                                                         1,
                                                         0),
                cattleSlaughter_valuePerHa_fitted_2006 = dplyr::if_else(cattleSlaughter_valuePerHa_fitted_2006 <= 0,
                                                                        aux.min.positive.cattleSlaughter.value.ha.fitted.2006,
                                                                        cattleSlaughter_valuePerHa_fitted_2006))
# clean environment
rm(reg.cattleValueperHa.2006, aux.min.positive.cattleSlaughter.value.ha.fitted.2006)

# match munis with sites
site.theta.2006 <- sf::st_intersection(calibration.1059SitesModel %>% dplyr::select(id),
                                       muniTheta.prepData       %>% dplyr::select(muni_code, muni_area, cattleSlaughter_valuePerHa_fitted_2006,
                                                                                  pasture_area_2006, d_theta_winsorized_2006))

# calculate muni areas inside each site
site.theta.2006$muni_site_area <-
  sf::st_area(site.theta.2006) %>%
  units::set_units(ha) %>%
  unclass()

# drop spatial feature
site.theta.2006 <-
  site.theta.2006 %>%
  sf::st_drop_geometry()

# calculate cattleSlaughter_valuePerHa_fitted and pastureArea_value by site (for each muni adjust the value by the share of the muni area inside the site)
aux.theta.2006 <-
  site.theta.2006 %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(theta2006_1059Sites = weighted.mean(cattleSlaughter_valuePerHa_fitted_2006/aux.price.2006, w = pasture_area_2006*(muni_site_area/muni_area), na.rm = T),
                   pasture_area_2006 = sum(pasture_area_2006*(muni_site_area/muni_area), na.rm = T),
                   d_theta_winsorized_2006 = min(d_theta_winsorized_2006, na.rm = T))

# add cattleSlaughter_valuePerHa_fitted and pastureArea_value to spatial variables
calibration.1059SitesModel <- dplyr::left_join(calibration.1059SitesModel, aux.theta.2006)

# clean environment
rm(aux.theta.2006)


# calculate average and SD theta using the values of 2006 and 2017
calibration.1059SitesModel <-
  calibration.1059SitesModel %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(theta_1059Sites = rowMeans(across(starts_with("theta20")), na.rm = T),
                thetaSD_1059Sites = apply(across(starts_with("theta20")), 1, sd)) %>%
  dplyr::ungroup()




# PRICE INITIAL CONDITION AND TRANSITIONS ------------------------------------------------------------------------------------------------------------

seriesPriceCattle.prepData <-
  seriesPriceCattle.prepData %>%
  dplyr::mutate(price_real_mon_cattle = price_real_mon_cattle/3.1966) %>%  # change from BRL to USD (commercial exchange rate - selling - average - monthly - january/2017 - ipeadata)
  dplyr::filter(year >= 1995, year <= 2017)  # select time period


# 2 PRICES (LOW X HIGH)
seriesPriceCattle.prepData <-
  seriesPriceCattle.prepData %>%
  dplyr::mutate(price_low = quantile(price_real_mon_cattle, 0.25), # define low price value as the 33th percentile
                price_high = quantile(price_real_mon_cattle, 0.75), # define high price value as the 66th percentile
                price_median = quantile(price_real_mon_cattle, 0.5),
                price_mean = mean(price_real_mon_cattle))


# create discretized version of the price series (high and low values only - using 25th and 75th percentiles as cut-offs)
seriesPriceCattle.prepData$d_high <- as.numeric(NA) # initialize dummy indicating if the price is high or low
seriesPriceCattle.prepData[1, "d_high"] <- 1 # initial value set to 1 because the first price is the highest of the series

for (i in 2:nrow(seriesPriceCattle.prepData)) {

  # a change from high to low only occurs if price reaches a value below the low
  if (seriesPriceCattle.prepData[i-1, "d_high"] == 1 & seriesPriceCattle.prepData[i, "price_real_mon_cattle"] > seriesPriceCattle.prepData[i, "price_low"]) {
    seriesPriceCattle.prepData[i, "d_high"] <- 1
  } else if (seriesPriceCattle.prepData[i-1, "d_high"] == 1 & seriesPriceCattle.prepData[i, "price_real_mon_cattle"] < seriesPriceCattle.prepData[i, "price_low"]) {
    seriesPriceCattle.prepData[i, "d_high"] <- 0
    # a change from low to high only occurs if price reaches a value above the high
  } else if (seriesPriceCattle.prepData[i-1, "d_high"] == 0 & seriesPriceCattle.prepData[i, "price_real_mon_cattle"] > seriesPriceCattle.prepData[i, "price_high"]) {
    seriesPriceCattle.prepData[i, "d_high"] <- 1
  } else if (seriesPriceCattle.prepData[i-1, "d_high"] == 0 & seriesPriceCattle.prepData[i, "price_real_mon_cattle"] < seriesPriceCattle.prepData[i, "price_high"]) {
    seriesPriceCattle.prepData[i, "d_high"] <- 0
  }
}

# construct discretized version of prices (2 states: high x low)
seriesPriceCattle.prepData <-
  seriesPriceCattle.prepData %>%
  dplyr::mutate(discrete_2prices = if_else(d_high == 1, price_high, price_low)) %>%
  dplyr::select(date, year, month, price_real_mon_cattle, discrete_2prices, price_high, price_median, price_mean, price_low)


# calculate probability transition matrix (give the same results as manually computing the number of consecutive prices at the same level divided by the ocurrence of that price level)
matrixTransition.2prices <- markovchain::markovchainFit(seriesPriceCattle.prepData$discrete_2prices)$estimate@transitionMatrix


# STORE PARAMETER VALUES
calibration.1059SitesModel <-
  calibration.1059SitesModel %>%
  dplyr::mutate(p_2017_1059Sites = max(seriesPriceCattle.prepData$price_high))





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# ORDER VARIABLES
calibration.1059SitesModel <-
  calibration.1059SitesModel %>%
  dplyr::select(id, z_2017_1059Sites, zbar_2017_1059Sites, x_2017_1059Sites, gamma_1059Sites, theta_1059Sites,
                d_theta_winsorized_2017, d_theta_winsorized_2006, pasture_area_2017, ends_with("_1059Sites"))


# POST-TREATMENT OVERVIEW
# summary(calibration.1059SitesModel)
# View(calibration.1059SitesModel)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(calibration.1059SitesModel,
     file = here::here("data/calibration/_deprecated",
                       "calibration_1059SitesModel.Rdata"))

# remove spatial feature
calibration.1059SitesModel <- calibration.1059SitesModel %>% sf::st_drop_geometry()

readr::write_csv(calibration.1059SitesModel,
                 file = here::here("data/calibration/_deprecated", "calibration_1059SitesModel.csv"))


# CLEAN TEMP DIR
terra::tmpFiles(current = T, remove = T)
gc()



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("code/calibration")






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------