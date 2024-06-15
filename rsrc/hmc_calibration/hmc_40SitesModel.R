
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (40 SITES MODEL)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -




# START TIMER
tictoc::tic(msg = "calibration_40SitesModel.R script", log = TRUE)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tmpdir = "data/_temp",
                      timer  = T)





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster_40Sites <- terra::rast(list.files("data/calibration/1043SitesModel/aux_tifs",
                                           pattern = "raster_",
                                           full.names = T))



# MUNI LEVEL SPATIAL SAMPLE
load("data/prepData/sampleMuniSpatial_prepData.Rdata")





# INITIAL CONDITIONS Z -------------------------------------------------------------------------------------------------------------------------------

# AGGREGATE FROM 1000 SITES TO 100 SITES
# transform shares to areas
raster_40Sites$amazonBiomeArea_ha_40Sites <- raster_40Sites$share_amazonBiome*raster_40Sites$pixelArea_ha
raster_40Sites$forestArea_1995_ha_40Sites <- raster_40Sites$share_forest_1995*raster_40Sites$pixelArea_ha
raster_40Sites$agriculturalUseArea_1995_ha_40Sites <- raster_40Sites$share_agriculturalUse_1995*raster_40Sites$pixelArea_ha
raster_40Sites$otherArea_1995_ha_40Sites <- raster_40Sites$share_other_1995*raster_40Sites$pixelArea_ha
raster_40Sites$forestArea_2017_ha_40Sites <- raster_40Sites$share_forest_2017*raster_40Sites$pixelArea_ha
raster_40Sites$agriculturalUseArea_2017_ha_40Sites <- raster_40Sites$share_agriculturalUse_2017*raster_40Sites$pixelArea_ha
raster_40Sites$otherArea_2017_ha_40Sites <- raster_40Sites$share_other_2017*raster_40Sites$pixelArea_ha
raster_40Sites$forestArea_2008_ha_40Sites <- raster_40Sites$share_forest_2008*raster_40Sites$pixelArea_ha
raster_40Sites$agriculturalUseArea_2008_ha_40Sites <- raster_40Sites$share_agriculturalUse_2008*raster_40Sites$pixelArea_ha
raster_40Sites$otherArea_2008_ha_40Sites <- raster_40Sites$share_other_2008*raster_40Sites$pixelArea_ha

# select area variables
# select area variables
raster_40Sites <- terra::subset(raster_40Sites,
                                c("amazonBiomeArea_ha_40Sites", "pixelArea_ha",
                                  "forestArea_1995_ha_40Sites", "agriculturalUseArea_1995_ha_40Sites", "otherArea_1995_ha_40Sites",
                                  "forestArea_2017_ha_40Sites", "agriculturalUseArea_2017_ha_40Sites", "otherArea_2017_ha_40Sites",
                                  "forestArea_2008_ha_40Sites", "agriculturalUseArea_2008_ha_40Sites", "otherArea_2008_ha_40Sites"))
# aggregate from 1000 sites to 100
raster_40Sites <- terra::aggregate(raster_40Sites, fact = 6, fun = sum, na.rm = T)

# extract variables as polygons, transform to sf, and project data for faster spatial manipulation
calibration_40SitesModel <- terra::as.polygons(raster_40Sites, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)

# transform share aggregate in area (ha)
calibration_40SitesModel <-
  calibration_40SitesModel %>%
  dplyr::mutate(zbar_1995_40Sites = agriculturalUseArea_1995_ha_40Sites + forestArea_1995_ha_40Sites,
                zbar_2017_40Sites = agriculturalUseArea_2017_ha_40Sites + forestArea_2017_ha_40Sites,
                zbar_2008_40Sites = agriculturalUseArea_2008_ha_40Sites + forestArea_2008_ha_40Sites) %>%
  dplyr::select(siteArea_ha_40Sites = pixelArea_ha, amazonBiomeArea_ha_40Sites,
                forestArea_1995_ha_40Sites, otherArea_1995_ha_40Sites,
                z_1995_40Sites = agriculturalUseArea_1995_ha_40Sites, zbar_1995_40Sites,
                forestArea_2017_ha_40Sites, otherArea_2017_ha_40Sites,
                z_2017_40Sites = agriculturalUseArea_2017_ha_40Sites, zbar_2017_40Sites,
                forestArea_2008_ha_40Sites, otherArea_2008_ha_40Sites,
                z_2008_40Sites = agriculturalUseArea_2008_ha_40Sites, zbar_2008_40Sites)

# remove sites with less than 1% of its are intersecting with the amazon biome
calibration_40SitesModel <-
  calibration_40SitesModel %>%
  dplyr::filter(amazonBiomeArea_ha_40Sites/siteArea_ha_40Sites >= 0.03)

# add id variable
calibration_40SitesModel$id <- 1:nrow(calibration_40SitesModel)



# PARAMETER GAMMA ------------------------------------------------------------------------------------------------------------------------------------


# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/prepData/muniTheta_prepData.Rdata")

muniTheta_prepData<-muniTheta_prepData %>%
  dplyr::mutate(co2e_ha_2017 = (agb_2017/2)*(44/12))




# DATA INPUT (2017)
# load pixel sample with biomass data
# Load pixelBiomass2017_prepData.Rdata




reg_gamma_2017 <-
  muniTheta_prepData  %>%
  lm(formula = log(co2e_ha_2017)  ~ log(historical_precip) + log(historical_temp) +log(lat)+log(lon), na.action = na.exclude)

summary(reg_gamma_2017)


muniTheta_prepData  <-   muniTheta_prepData %>%
  dplyr::mutate(co2e_ha_2017_fitted = exp(predict(reg_gamma_2017, .)))



# match minicells with Sites
site_gamma2017 <-
  sf::st_intersection(calibration_40SitesModel %>% dplyr::select(id),
              muniTheta_prepData %>% dplyr::select(muni_code, muni_area, co2e_ha_2017,co2e_ha_2017_fitted,historical_precip,historical_temp,lat,lon)) 
# sf::st_drop_geometry()

site_gamma2017$muni_site_area <-
  sf::st_area(site_gamma2017) %>%
  units::set_units(ha) %>%
  unclass()

# drop spatial feature
site_gamma2017 <-
  site_gamma2017 %>%
  sf::st_drop_geometry()

# calculate average carbon density on primary forest areas by site
aux_gamma2017 <-
  site_gamma2017 %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(gamma2017_40Sites = weighted.mean(co2e_ha_2017_fitted, w=muni_site_area, na.rm = T))




# add gamma_40Sites to spatial variables
calibration_40SitesModel <- dplyr::left_join(calibration_40SitesModel, aux_gamma2017)


# clean environment
rm(aux_gamma2017)


# identify adjacent neighbors
aux_neighbors <- sf::st_is_within_distance(calibration_40SitesModel, calibration_40SitesModel, dist = 100, remove_self = TRUE)

# impute values for missing gammas based on the average of adjacent neighbors
calibration_40SitesModel <-
  calibration_40SitesModel %>%
  dplyr::mutate(gamma2017_40Sites = dplyr::if_else(is.na(gamma2017_40Sites),
                                                   apply(aux_neighbors, 1, function(i){mean(.$gamma2017_40Sites[i], na.rm = TRUE)}),
                                                   gamma2017_40Sites))


# set baseline gamma and gammaSD as the data from 2017 and calculate alternative  gamma gammaSD based on the mean and sd of gamma2010,gamma2017, and gamma2018
calibration_40SitesModel <-
  calibration_40SitesModel %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(gamma_40Sites = rowMeans(across(c("gamma2017_40Sites")))) %>%
  dplyr::ungroup()






# PARAMETER ALPHA ------------------------------------------------------------------------------------------------------------------------------------

# estimate of alpha same as in the global model
calibration_40SitesModel <-
  calibration_40SitesModel %>%
  dplyr::mutate(alpha_40Sites = 1 - (1-0.99)^(1/100))





# PARAMETER KAPPA ----------------------------------------------------------------------------------------------------------------------------------------

# DATA INPUT
# load pixel sample with biomass data
load("data/prepData/stateEmission_prepData.Rdata")


# DATA MANIPULATION
# calculate average of net emission factor from agricultural use across years and states
avg_netEmissionFactor <-
  stateEmission_prepData %>%
  dplyr::ungroup() %>%
  dplyr::summarise(netEmissionFactor_co2e_ha = sum(netEmission_co2e)/sum(agriculturalUse_area)) %>%
  pull(netEmissionFactor_co2e_ha)


# estimate of kappa same as in the global model
calibration_40SitesModel <-
  calibration_40SitesModel %>%
  dplyr::mutate(kappa_40Sites = avg_netEmissionFactor)

# clean environment
rm(avg_netEmissionFactor)





# PARAMETER ZETA ----------------------------------------------------------------------------------------------------------------------------------------

# zeta is calibrated such that the marginal cost of changing land use (zeta*forestToPastureTransitionArea) matches the forest to pasture transition cost >
# estimated by Araujo, Costa and Sant'Anna (2022), reported on column 4 of table 4 (on the right). We transform to dollars using an FX rate of 4.14 (December 2019) >
# as in the paper: 1614.54/4.14 = 390 USD/ha. The paper also calculates that 6.5% of the forest area was converted to pasture between 2008 and 2017. >
aux_transitionCost <- 1614.54/4.14

# The forest area in 2008 represented 72% of the Legal Amazon area, which covers 501,506,775 ha, so the transition in hectares is >
#  0.065*0.72*501,506,775 = 23,470,517, resulting in an annual average of 2,347,052 ha.
aux_transitionArea <- (0.065*0.72*501506775)/(2017-2008+1)

zeta <- aux_transitionCost/aux_transitionArea
zeta_alt <- 483/aux_transitionArea # Alternative value based on a quote from (https://www.otempo.com.br/brasil/investigacoes-revelam-quadrilhas-e-ganho-milionario-por-tras-do-desmate-1.2229571)

# estimate of zeta same as in the global model
calibration_40SitesModel <-
  calibration_40SitesModel %>%
  dplyr::mutate(zeta_40Sites = zeta,
                zeta_alt_40Sites = zeta_alt)





# INITIAL CONDITIONS X -------------------------------------------------------------------------------------------------------------------------------

# x_2017_40Sites estimated as in the old way of global model, just considering the stock of carbon stored in forest areas assuming that all forests are primary
calibration_40SitesModel <-
  calibration_40SitesModel %>%
  dplyr::mutate(x_2017_40Sites = gamma_40Sites*(zbar_2017_40Sites-z_2017_40Sites),
                x_1995_40Sites = gamma_40Sites*(zbar_1995_40Sites-z_1995_40Sites),
                x_2008_40Sites = gamma_40Sites*(zbar_2008_40Sites-z_2008_40Sites))





# PARAMETER THETA ------------------------------------------------------------------------------------------------------------------------------------

distance_data <-
  read_excel("data/raw/ipea/distance_to_capital/ipeadata[21-08-2023-01-28].xls")

distance_data$muni_code <- as.numeric(distance_data$muni_code)
# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/prepData/muniTheta_prepData.Rdata")

# load cattle price series
load("data/prepData/seriesPriceCattle_prepData.Rdata")


# DATA MANIPULATION

# EXTRACT AVERAGE 2017 PRICE (use real prices because it is normalized to 2017 )
aux_price_2017 <-
  seriesPriceCattle_prepData %>%
  dplyr::filter(year == 2017) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(mean_price_2017 = mean(price_real_mon_cattle)/3.192) %>% # BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
  dplyr::pull(mean_price_2017)



# Remove rows from attribute data
muniTheta_prepData_data <- as.data.frame(muniTheta_prepData)  # Convert to regular dataframe
muniTheta_prepData_data <- muniTheta_prepData_data[-c(142, 106, 112), ]

# Remove geometries
geo_backup <- st_geometry(muniTheta_prepData)
geo_backup <- geo_backup[-c(142, 106, 112)]


predicted_values <-
  read_excel("data/raw/ipea/farm_gate_price/farm_gate_price.xlsx")

# Combine back into an sf object
muniTheta_prepData <- st_sf(muniTheta_prepData_data, geometry = geo_backup)

# 2. Merging the cleaned muniTheta_prepData with my_data

# Convert to non-spatial dataframe for the merge
muniTheta_no_geo <- as.data.frame(muniTheta_prepData)

# Perform the merge
merged_data <- left_join(muniTheta_no_geo, distance_data, by = "muni_code")

# Reattach the geometry
merged_data_sf <- st_sf(merged_data, geometry = geo_backup)

muniTheta_prepData<-merged_data_sf

merged_data <- muniTheta_prepData %>%
  left_join(predicted_values, by = "muni_code") %>%
  mutate(cattleSlaughter_farmGatePrice_2017 = ifelse(is.na(cattleSlaughter_farmGatePrice_2017),
                                                     average_weighted_price,
                                                     cattleSlaughter_farmGatePrice_2017))


muniTheta_prepData<-merged_data


muniTheta_prepData<- muniTheta_prepData %>%
  dplyr::filter(!is.na(distance))

muniTheta_prepData_filtered <- muniTheta_prepData %>%
  dplyr::filter(cattleSlaughter_valuePerHa_2017 > 0) # Exclude zeros

# REGRESSION - CATTLE VALUE (2017)

# cattle value per ha


reg_cattleValueperHa_2017 <-
  muniTheta_prepData_filtered  %>%
  lm(formula = log(cattleSlaughter_valuePerHa_2017) ~  historical_precip+ historical_temp + I(historical_temp^2)
     + lat+I(lat^2)+distance+log(cattleSlaughter_farmGatePrice_2017), na.action = na.exclude, weights = pasture_area_2017)

# regression results
summary(reg_cattleValueperHa_2017)

# extract fitted values
muniTheta_prepData <-
  muniTheta_prepData %>%
  dplyr::mutate(cattleSlaughter_valuePerHa_fitted_2017 = exp(predict(reg_cattleValueperHa_2017, .)))

# extract minimum positive fitted value
aux_min_positive_cattleSlaughter_value_ha_fitted_2017 <-
  muniTheta_prepData %>%
  dplyr::filter(cattleSlaughter_valuePerHa_fitted_2017 > 0) %>%
  dplyr::pull(cattleSlaughter_valuePerHa_fitted_2017) %>%
  min()

# winsorize cattleSlaughter_valuePerHa_fitted
muniTheta_prepData <-
  muniTheta_prepData %>%
  dplyr::mutate(d_theta_winsorized_2017 = dplyr::if_else(cattleSlaughter_valuePerHa_fitted_2017 <= 0,
                                                         1,
                                                         0))
# clean environment
rm(reg_cattleValueperHa_2017, aux_min_positive_cattleSlaughter_value_ha_fitted_2017)

# match munis with sites
site_theta_2017 <- sf::st_intersection(calibration_40SitesModel %>% dplyr::select(id),
                                       muniTheta_prepData       %>% dplyr::select(muni_code, muni_area, cattleSlaughter_valuePerHa_fitted_2017,
                                                                                  pasture_area_2017, d_theta_winsorized_2017))



# calculate muni areas inside each site
site_theta_2017$muni_site_area <-
  sf::st_area(site_theta_2017) %>%
  units::set_units(ha) %>%
  unclass()

# drop spatial feature
site_theta_2017 <-
  site_theta_2017 %>%
  sf::st_drop_geometry()

# calculate cattleSlaughter_valuePerHa_fitted and pastureArea_value by site (for each muni adjust the value by the share of the muni area inside the site)
aux_theta_2017 <-
  site_theta_2017 %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(theta2017_40Sites = weighted.mean(cattleSlaughter_valuePerHa_fitted_2017/aux_price_2017, w = muni_site_area, na.rm = T),
                   pasture_area_2017 = sum(pasture_area_2017*(muni_site_area/muni_area), na.rm = T),
                   d_theta_winsorized_2017 = min(d_theta_winsorized_2017, na.rm = T))

# add cattleSlaughter_valuePerHa_fitted and pastureArea_value to spatial variables
calibration_40SitesModel <- dplyr::left_join(calibration_40SitesModel, aux_theta_2017)

# clean environment
rm(aux_theta_2017)




# calculate average and SD theta using the values of 2006 and 2017
calibration_40SitesModel <-
  calibration_40SitesModel %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(theta_40Sites = rowMeans(across(starts_with("theta20")), na.rm = T)) %>%
  dplyr::ungroup()




# PRICE INITIAL CONDITION AND TRANSITIONS ------------------------------------------------------------------------------------------------------------

seriesPriceCattle_prepData <-
  seriesPriceCattle_prepData %>%
  dplyr::mutate(price_real_mon_cattle = price_real_mon_cattle/3.1966) %>%  # change from BRL to USD (commercial exchange rate - selling - average - monthly - january/2017 - ipeadata)
  dplyr::filter(year >= 1995, year <= 2017)  # select time period


# 2 PRICES (LOW X HIGH)
seriesPriceCattle_prepData <-
  seriesPriceCattle_prepData %>%
  dplyr::mutate(price_low = quantile(price_real_mon_cattle, 0.25), # define low price value as the 33th percentile
                price_high = quantile(price_real_mon_cattle, 0.75), # define high price value as the 66th percentile
                price_median = quantile(price_real_mon_cattle, 0.5),
                price_mean = mean(price_real_mon_cattle))


# create discretized version of the price series (high and low values only - using 25th and 75th percentiles as cut-offs)
seriesPriceCattle_prepData$d_high <- as.numeric(NA) # initialize dummy indicating if the price is high or low
seriesPriceCattle_prepData[1, "d_high"] <- 1 # initial value set to 1 because the first price is the highest of the series

for (i in 2:nrow(seriesPriceCattle_prepData)) {

  # a change from high to low only occurs if price reaches a value below the low
  if (seriesPriceCattle_prepData[i-1, "d_high"] == 1 & seriesPriceCattle_prepData[i, "price_real_mon_cattle"] > seriesPriceCattle_prepData[i, "price_low"]) {
    seriesPriceCattle_prepData[i, "d_high"] <- 1
  } else if (seriesPriceCattle_prepData[i-1, "d_high"] == 1 & seriesPriceCattle_prepData[i, "price_real_mon_cattle"] < seriesPriceCattle_prepData[i, "price_low"]) {
    seriesPriceCattle_prepData[i, "d_high"] <- 0
    # a change from low to high only occurs if price reaches a value above the high
  } else if (seriesPriceCattle_prepData[i-1, "d_high"] == 0 & seriesPriceCattle_prepData[i, "price_real_mon_cattle"] > seriesPriceCattle_prepData[i, "price_high"]) {
    seriesPriceCattle_prepData[i, "d_high"] <- 1
  } else if (seriesPriceCattle_prepData[i-1, "d_high"] == 0 & seriesPriceCattle_prepData[i, "price_real_mon_cattle"] < seriesPriceCattle_prepData[i, "price_high"]) {
    seriesPriceCattle_prepData[i, "d_high"] <- 0
  }
}

# construct discretized version of prices (2 states: high x low)
seriesPriceCattle_prepData <-
  seriesPriceCattle_prepData %>%
  dplyr::mutate(discrete_2prices = if_else(d_high == 1, price_high, price_low)) %>%
  dplyr::select(date, year, month, price_real_mon_cattle, discrete_2prices, price_high, price_median, price_mean, price_low)


# calculate probability transition matrix (give the same results as manually computing the number of consecutive prices at the same level divided by the ocurrence of that price level)
matrixTransition_2prices <- markovchain::markovchainFit(seriesPriceCattle_prepData$discrete_2prices)$estimate@transitionMatrix


# STORE PARAMETER VALUES
calibration_40SitesModel <-
  calibration_40SitesModel %>%
  dplyr::mutate(p_2017_40Sites = max(seriesPriceCattle_prepData$price_high))





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# ORDER VARIABLES
calibration_40SitesModel <-
  calibration_40SitesModel %>%
  dplyr::select(id, z_2017_40Sites, zbar_2017_40Sites, x_2017_40Sites, gamma_40Sites, theta_40Sites,
                d_theta_winsorized_2017,pasture_area_2017, ends_with("_40Sites"))






# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(calibration_40SitesModel,
     file = "data/calibration/hmc/hmc_40SitesModel.Rdata")

# remove spatial feature
calibration_40SitesModel <- calibration_40SitesModel %>% sf::st_drop_geometry()

readr::write_csv(calibration_40SitesModel,
                 file = "data/calibration/hmc/hmc_40SitesModel.csv")


# CLEAN TEMP DIR
terra::tmpFiles(current = T, remove = T)
gc()



# END TIMER
tictoc::toc(log = TRUE)



# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
