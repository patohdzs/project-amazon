
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (78 SITES MODEL)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "calibration_78SitesModel.R script", log = T)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster.78Sites <- terra::rast(list.files(here::here("data/calibration/1055SitesModel/aux_tifs"),
                                          pattern = "raster_",
                                          full.names = T))


# MUNI LEVEL SPATIAL SAMPLE
load(here::here("data/calibration/prepData/sampleMuniSpatial_prepData.Rdata"))





# INITIAL CONDITIONS Z -------------------------------------------------------------------------------------------------------------------------------

# AGGREGATE FROM 1000 SITES TO 43 SITES
# transform shares to areas
raster.78Sites$amazonBiomeArea_ha_78Sites <- raster.78Sites$share_amazonBiome*raster.78Sites$pixelArea_ha
raster.78Sites$forestArea_1995_ha_78Sites <- raster.78Sites$share_forest_1995*raster.78Sites$pixelArea_ha
raster.78Sites$agriculturalUseArea_1995_ha_78Sites <- raster.78Sites$share_agriculturalUse_1995*raster.78Sites$pixelArea_ha
raster.78Sites$otherArea_1995_ha_78Sites <- raster.78Sites$share_other_1995*raster.78Sites$pixelArea_ha
raster.78Sites$forestArea_2017_ha_78Sites <- raster.78Sites$share_forest_2017*raster.78Sites$pixelArea_ha
raster.78Sites$agriculturalUseArea_2017_ha_78Sites <- raster.78Sites$share_agriculturalUse_2017*raster.78Sites$pixelArea_ha
raster.78Sites$otherArea_2017_ha_78Sites <- raster.78Sites$share_other_2017*raster.78Sites$pixelArea_ha
raster.78Sites$forestArea_2008_ha_78Sites <- raster.78Sites$share_forest_2008*raster.78Sites$pixelArea_ha
raster.78Sites$agriculturalUseArea_2008_ha_78Sites <- raster.78Sites$share_agriculturalUse_2008*raster.78Sites$pixelArea_ha
raster.78Sites$otherArea_2008_ha_78Sites <- raster.78Sites$share_other_2008*raster.78Sites$pixelArea_ha
# select area variables
raster.78Sites <- terra::subset(raster.78Sites,
                                c("amazonBiomeArea_ha_78Sites", "pixelArea_ha",
                                  "forestArea_1995_ha_78Sites", "agriculturalUseArea_1995_ha_78Sites", "otherArea_1995_ha_78Sites",
                                  "forestArea_2017_ha_78Sites", "agriculturalUseArea_2017_ha_78Sites", "otherArea_2017_ha_78Sites",
                                  "forestArea_2008_ha_78Sites", "agriculturalUseArea_2008_ha_78Sites", "otherArea_2008_ha_78Sites"))
# aggregate from 1000 sites to 100
# aggregate from 1000 sites to 43
raster.78Sites <- terra::aggregate(raster.78Sites, fact = 4, fun = sum, na.rm = T)

# extract variables as polygons, transform to sf, and project data for faster spatial manipulation
calibration.78SitesModel <- terra::as.polygons(raster.78Sites, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)

# transform share aggregate in area (ha)
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::mutate(zbar_1995_78Sites = agriculturalUseArea_1995_ha_78Sites + forestArea_1995_ha_78Sites,
                zbar_2017_78Sites = agriculturalUseArea_2017_ha_78Sites + forestArea_2017_ha_78Sites,
                zbar_2008_78Sites = agriculturalUseArea_2008_ha_78Sites + forestArea_2008_ha_78Sites) %>%
  dplyr::select(amazonBiomeArea_ha_78Sites, siteArea_ha_78Sites = pixelArea_ha,
                forestArea_1995_ha_78Sites,
                z_1995_78Sites = agriculturalUseArea_1995_ha_78Sites, zbar_1995_78Sites,
                forestArea_2017_ha_78Sites, otherArea_2017_ha_78Sites,
                z_2017_78Sites = agriculturalUseArea_2017_ha_78Sites, zbar_2017_78Sites,
                forestArea_2008_ha_78Sites, otherArea_2008_ha_78Sites,
                z_2008_78Sites = agriculturalUseArea_2008_ha_78Sites, zbar_2008_78Sites)


# remove sites with less than 2% of its are intersecting with the amazon biome
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::filter(amazonBiomeArea_ha_78Sites/siteArea_ha_78Sites >= 0.03)

# add id variable
calibration.78SitesModel$id <- 1:nrow(calibration.78SitesModel)




# PARAMETER GAMMA ------------------------------------------------------------------------------------------------------------------------------------



# DATA INPUT (2017)
# load pixel sample with biomass data
# Load pixelBiomass2017_prepData.Rdata
load(here::here("data/calibration/prepData/pixelBiomass2017_prepData.Rdata"))


# select minicells of primary forest with co2 information and transform to spatial points
pixelBiomass2017.prepData <-
  pixelBiomass2017.prepData %>%
  sf::st_transform(sf::st_crs(calibration.78SitesModel)) %>%
  dplyr::mutate(co2e_ha_2017 = (agb_2017/2)*(44/12)) %>%
  dplyr::filter(co2e_ha_2017 > 0, !is.na(co2e_ha_2017))

# match minicells with Sites
site.gamma2017 <-
  sf::st_join(calibration.78SitesModel %>% dplyr::select(id),
              pixelBiomass2017.prepData %>% dplyr::select(co2e_ha_2017)) %>%
  sf::st_drop_geometry()

# calculate average carbon density on primary forest areas by site
aux.gamma2017 <-
  site.gamma2017 %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(gamma2017_78Sites = mean(co2e_ha_2017, na.rm = T))

# add gamma_78Sites to spatial variables
calibration.78SitesModel <- dplyr::left_join(calibration.78SitesModel, aux.gamma2017)


# clean environment
rm(pixelBiomass2017.prepData, aux.gamma2017)


# identify adjacent neighbors
aux.neighbors <- sf::st_is_within_distance(calibration.78SitesModel, calibration.78SitesModel, dist = 100, remove_self = TRUE)

# impute values for missing gammas based on the average of adjacent neighbors
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::mutate(gamma2017_78Sites = dplyr::if_else(is.na(gamma2017_78Sites),
                                                   apply(aux.neighbors, 1, function(i){mean(.$gamma2017_78Sites[i], na.rm = TRUE)}),
                                                   gamma2017_78Sites))


# set baseline gamma and gammaSD as the data from 2017 and calculate alternative  gamma gammaSD based on the mean and sd of gamma2010,gamma2017, and gamma2018
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(gamma_78Sites = rowMeans(across(c("gamma2017_78Sites")))) %>%
  #gammaSD_10Sites = apply(across(c("gamma2010_10Sites", "gamma2017_10Sites", "gamma2018_10Sites")), 1, sd)) %>%
  dplyr::ungroup()








# PARAMETER ALPHA ------------------------------------------------------------------------------------------------------------------------------------

# estimate of alpha same as in the global model
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::mutate(alpha_78Sites = 1 - (1-0.99)^(1/100))





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
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::mutate(kappa_78Sites = avg.netEmissionFactor)

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
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::mutate(zeta_78Sites = zeta,
                zeta_alt_78Sites = zeta_alt)





# INITIAL CONDITIONS X -------------------------------------------------------------------------------------------------------------------------------

# x_2017_78Sites estimated as in the old way of global model, just considering the stock of carbon stored in forest areas assuming that all forests are primary
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::mutate(x_2017_78Sites = gamma_78Sites*(zbar_2017_78Sites-z_2017_78Sites),
                x_1995_78Sites = gamma_78Sites*(zbar_1995_78Sites-z_1995_78Sites),
                x_2008_78Sites = gamma_78Sites*(zbar_2008_78Sites-z_2008_78Sites))





# PARAMETER THETA ------------------------------------------------------------------------------------------------------------------------------------

distance_data <-
  read_excel("data/calibration/ipeadata[21-08-2023-01-28].xls")

distance_data$muni_code <- as.numeric(distance_data$muni_code)
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



# Remove rows from attribute data
muniTheta.prepData_data <- as.data.frame(muniTheta.prepData)  # Convert to regular dataframe
muniTheta.prepData_data <- muniTheta.prepData_data[-c(142, 106, 112), ]

# Remove geometries
geo_backup <- st_geometry(muniTheta.prepData)
geo_backup <- geo_backup[-c(142, 106, 112)]


predicted_values <-
  read_excel("data/calibration/farm_gate_price.xlsx")

# Combine back into an sf object
muniTheta.prepData <- st_sf(muniTheta.prepData_data, geometry = geo_backup)

# 2. Merging the cleaned muniTheta.prepData with my_data

# Convert to non-spatial dataframe for the merge
muniTheta_no_geo <- as.data.frame(muniTheta.prepData)

# Perform the merge
merged_data <- left_join(muniTheta_no_geo, distance_data, by = "muni_code")

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

# REGRESSION - CATTLE VALUE (2017)

# cattle value per ha


reg.cattleValueperHa.2017 <-
  muniTheta.prepData_filtered  %>%
  lm(formula = log(cattleSlaughter_valuePerHa_2017) ~  historical_precip+ historical_temp + I(historical_temp^2)
     + lat+I(lat^2)+distance+cattleSlaughter_farmGatePrice_2017, na.action = na.exclude, weights = pasture_area_2017)

# regression results
summary(reg.cattleValueperHa.2017)

# extract fitted values
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(cattleSlaughter_valuePerHa_fitted_2017 = exp(predict(reg.cattleValueperHa.2017, .)))

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
                                                         0))
# clean environment
rm(reg.cattleValueperHa.2017, aux.min.positive.cattleSlaughter.value.ha.fitted.2017)

# match munis with sites
site.theta.2017 <- sf::st_intersection(calibration.78SitesModel %>% dplyr::select(id),
                                       muniTheta.prepData       %>% dplyr::select(muni_code, muni_area, cattleSlaughter_valuePerHa_fitted_2017,
                                                                                  pasture_area_2017, d_theta_winsorized_2017,zbar_2017_muni))



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
  dplyr::filter(!is.na(zbar_2017_muni)) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(theta2017_78Sites = weighted.mean(cattleSlaughter_valuePerHa_fitted_2017/aux.price.2017, w = muni_site_area, na.rm = T),
                   pasture_area_2017 = sum(pasture_area_2017*(muni_site_area/muni_area), na.rm = T),
                   d_theta_winsorized_2017 = min(d_theta_winsorized_2017, na.rm = T))

# add cattleSlaughter_valuePerHa_fitted and pastureArea_value to spatial variables
calibration.78SitesModel <- dplyr::left_join(calibration.78SitesModel, aux.theta.2017)

# clean environment
rm(aux.theta.2017)




# calculate average and SD theta using the values of 2006 and 2017
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(theta_78Sites = rowMeans(across(starts_with("theta20")), na.rm = T)) %>%
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
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::mutate(p_2017_78Sites = max(seriesPriceCattle.prepData$price_high))





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# ORDER VARIABLES
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::select(id, z_2017_78Sites, zbar_2017_78Sites, x_2017_78Sites, gamma_78Sites, theta_78Sites,
                d_theta_winsorized_2017,pasture_area_2017, ends_with("_78Sites"))



# POST-TREATMENT OVERVIEW
# summary(calibration.78SitesModel)
# View(calibration.78SitesModel)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(calibration.78SitesModel,
     file = paste(getwd(), "data/calibration/78SitesModel", "calibration_78SitesModel.Rdata", sep = "/"))

# remove spatial feature
calibration.78SitesModel <- calibration.78SitesModel %>% sf::st_drop_geometry()

# Save calibration.24SitesModel as CSV
readr::write_csv(calibration.78SitesModel,
                 file = paste(getwd(), "data/calibration/78SitesModel", "calibration_78SitesModel.csv", sep = "/"))


# CLEAN TEMP DIR
terra::tmpFiles(current = T, remove = T)
gc()



# END TIMER
tictoc::toc(log = T)



# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
