
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (GLOBAL MODEL)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

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



# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")



# START TIMER
tictoc::tic(msg = "calibration_globalModel.R script", log = T)

# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = paste(getwd(), "data", "_temp", sep = "/"))



# INITIAL CONDITIONS Z -------------------------------------------------------------------------------------------------------------------------------

# load pixel sample with mapbiomas categories data
load(here::here("data/calibration/prepData/pixelCategories_prepData.Rdata"))

# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster.variables <-  terra::rast(list.files(paste(getwd(), "data/calibration/1043SitesModel/aux_tifs", sep = "/"),
                                         pattern = "raster_",
                                         full.names = T))


# MAPBIOMAS VARIABLES + AMAZON BIOME + PIXEL AREA (Z_1000Sites CONSTRUCTION)
# extract variables as polygons, transform to sf, and project data for faster spatial manipulation
calibration.1000SitesModel <- terra::as.polygons(raster.variables, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)

# add id variable
calibration.1000SitesModel$id <- 1:length(calibration.1000SitesModel$geometry)


# transform share variables in area (ha)
calibration.globalModel <-
  calibration.1000SitesModel %>%
  sf::st_drop_geometry() %>%
  dplyr::summarise(amazonBiomeArea_ha_global = sum(share_amazonBiome*pixelArea_ha),
                forestArea_2017_ha_global = sum(share_forest_2017*pixelArea_ha),
                z_2017_global = sum(share_agriculturalUse_2017*pixelArea_ha),
                otherArea_2017_ha_global = sum(share_other_2017*pixelArea_ha),
                zbar_2017_global = sum(forestArea_2017_ha_global + z_2017_global),
                forestArea_2008_ha_global = sum(share_forest_2008*pixelArea_ha),
                z_2008_global = sum(share_agriculturalUse_2008*pixelArea_ha),
                otherArea_2008_ha_global = sum(share_other_2008*pixelArea_ha),
                zbar_2008_global = sum(forestArea_2008_ha_global + z_2008_global),
                forestArea_1995_ha_global = sum(share_forest_1995*pixelArea_ha),
                z_1995_global = sum(share_agriculturalUse_1995*pixelArea_ha),
                otherArea_1995_ha_global = sum(share_other_1995*pixelArea_ha),
                zbar_1995_global = sum(forestArea_1995_ha_global + z_1995_global)
                )





# PARAMETER GAMMA ------------------------------------------------------------------------------------------------------------------------------------




# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/calibration/prepData/muniTheta_prepData.Rdata")

muniTheta.prepData<-muniTheta.prepData %>%
  dplyr::mutate(co2e_ha_2017 = (agb_2017/2)*(44/12))



# DATA INPUT (2017)
# load pixel sample with biomass data
# Load pixelBiomass2017_prepData.Rdata




reg.gamma.2017 <-
  muniTheta.prepData  %>%
  lm(formula = log(co2e_ha_2017)  ~ log(historical_precip) + log(historical_temp) +log(lat)+log(lon), na.action = na.exclude)

summary(reg.gamma.2017)


muniTheta.prepData  <-   muniTheta.prepData %>%
  dplyr::mutate(co2e_ha_2017_fitted = exp(predict(reg.gamma.2017, .)))


aux.gamma <-
  muniTheta.prepData %>%
  sf::st_drop_geometry() %>%
  dplyr::select(muni_code, muni_area, biomeAmazon_share, co2e_ha_2017_fitted) %>%
  dplyr::ungroup() %>%
  dplyr::summarise(gamma    = weighted.mean(co2e_ha_2017_fitted, w = biomeAmazon_share*muni_area))


# STORE PARAMETER VALUES
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::bind_cols(dplyr::tibble(gamma2017_global = aux.gamma$gamma))





# set baseline gamma and gammaSD as the data from 2017 and calculate alternative  gamma gammaSD based on the mean and sd of gamma2010,gamma2017, and gamma2018
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::mutate(gamma_global = rowMeans(across(c("gamma2017_global")))) %>%
  dplyr::ungroup()




# PARAMETER A ----------------------------------------------------------------------------------------------------------------------------------------

# DATA MANIPULATION
# estimate of a to make convergence (0.99*b/a) happens in 100 years (time based on Heinrich et al (2021))
a <- 1 - (1-0.99)^(1/100)


# STORE PARAMETER VALUES
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::bind_cols(dplyr::tibble(a_global = a))

# clean environment
rm(a)





# PARAMETER KAPPA ----------------------------------------------------------------------------------------------------------------------------------------

# DATA INPUT
# load pixel sample with biomass data
load(paste(getwd(), "data/calibration/prepData/stateEmission_prepData.Rdata", sep = "/"))


# DATA MANIPULATION
# calculate average of net emission factor from agricultural use across years and states
avg.netEmissionFactor <-
  stateEmission.prepData %>%
  dplyr::ungroup() %>%
  dplyr::summarise(netEmissionFactor_co2e_ha = sum(netEmission_co2e)/sum(agriculturalUse_area)) %>%
  pull(netEmissionFactor_co2e_ha)


# STORE PARAMETER VALUES
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::bind_cols(dplyr::tibble(kappa_global = avg.netEmissionFactor))


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


# STORE PARAMETER VALUES
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::bind_cols(dplyr::tibble(zeta_global = zeta,
                                 zeta_alt_global = zeta_alt))


# clean environment
rm(zeta, zeta_alt, aux.transitionCost, aux.transitionArea)





# INITIAL CONDITION X -------------------------------------------------------------------------------------------------------------------------

# x_2017_1000Sites estimated as in the old way of global model, just considering the stock of carbon stored in forest areas assuming that all forests are primary
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::mutate(x_2017_global = gamma_global*(zbar_2017_global-z_2017_global),
                x_1995_global = gamma_global*(zbar_1995_global-z_1995_global),
                x_2008_global = gamma_global*(zbar_2008_global-z_2008_global))





# PARAMETER THETA ------------------------------------------------------------------------------------------------------------------------------------

# DATA INPUT

distance_data <-
  read_excel("data/raw2clean/distance_to_capital/ipeadata[21-08-2023-01-28].xls")
distance_data$muni_code <- as.numeric(distance_data$muni_code)


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
  read_excel("data/raw2clean/farm_gate_price/farm_gate_price.xlsx")



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
     + lat+I(lat^2)+distance+log(cattleSlaughter_farmGatePrice_2017), na.action = na.exclude, weights = pasture_area_2017)

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



# POSSIBLE VALUES OF THETA (2017)
aux.theta <-
  muniTheta.prepData %>%
  sf::st_drop_geometry() %>%
  dplyr::select(muni_code, muni_area, biomeAmazon_share, cattleSlaughter_valuePerHa_fitted_2017, pasture_area_2017) %>%
  dplyr::ungroup() %>%
  dplyr::summarise(cattleSlaughter_valuePerHa_fitted    = weighted.mean(cattleSlaughter_valuePerHa_fitted_2017, w = biomeAmazon_share*muni_area),
                   theta = cattleSlaughter_valuePerHa_fitted/(aux.price.2017))


# STORE PARAMETER VALUES
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::bind_cols(dplyr::tibble(theta2017_global = aux.theta$theta))



# calculate average and SD theta using the values of 2006 and 2017
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::mutate(theta_global = rowMeans(across(starts_with("theta20"))))



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
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::bind_cols(dplyr::tibble(p_2017_global = max(seriesPriceCattle.prepData$price_high)))




# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(calibration.globalModel$gamma_global) <- "average CO2e density (Mg/ha) on primary forests"
sjlabelled::set_label(calibration.globalModel$a_global) <- "estimate of a to make convergence (0.99*b/a) happens in 100 years (time based on Heinrich et al (2021))"
sjlabelled::set_label(calibration.globalModel$kappa_global) <- "average net emission per hectare of agricultural use across years and states"
sjlabelled::set_label(calibration.globalModel$zbar_2017_global) <- "maximum value of z, sum of areas of forest and agricultural use in 2017 in the Amazon Biome (hectares)"
sjlabelled::set_label(calibration.globalModel$z_2017_global) <- "area of agricultural use in 2017 in the Amazon Biome (hectares)"
sjlabelled::set_label(calibration.globalModel$zeta_global) <- "cost of transition from forest to pasture (per hectare) using estimate from Araujo, Costa and Sant'Anna (2022)"
sjlabelled::set_label(calibration.globalModel$zeta_alt_global) <- "cost of transition from forest to pasture (per hectare) using estimate from Araujo, Costa and Sant'Anna (2022) and (https://www.otempo.com.br/brasil/investigacoes-revelam-quadrilhas-e-ganho-milionario-por-tras-do-desmate-1.2229571)"
sjlabelled::set_label(calibration.globalModel$x_2017_global) <- "total stock of carbon stored in the Amazon forest (Mg CO2e) (using carbon accumulation equation with gamma_global)"
sjlabelled::set_label(calibration.globalModel$theta_global) <- "cattle productivity using average predicted value of cattle sold for slaughter per hectare of pasture area adjusted by p_2017_global"
sjlabelled::set_label(calibration.globalModel$p_2017_global) <- "cattle price in 2017 (high level)"




# POST-TREATMENT OVERVIEW
# summary(calibration.globalModel)
# View(calibration.globalModel)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------



# export csv files
readr::write_csv(calibration.globalModel,
                 file = paste(getwd(), "data/hmc", "hmc_globalModel.csv", sep = "/"))



# END TIMER
tictoc::toc(log = T)






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
