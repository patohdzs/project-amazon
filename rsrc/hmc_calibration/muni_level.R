
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








# PARAMETER GAMMA ------------------------------------------------------------------------------------------------------------------------------------


# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/calibration/prepData/muniTheta_prepData_new.Rdata")

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
  dplyr::mutate(gamma = exp(predict(reg.gamma.2017, .)))




# INITIAL CONDITIONS X -------------------------------------------------------------------------------------------------------------------------------

# x_2017_78Sites estimated as in the old way of global model, just considering the stock of carbon stored in forest areas assuming that all forests are primary
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(x_2017_muni = gamma*(zbar_2017_muni-z_2017_muni),
                x_1995_muni = gamma*(zbar_1995_muni-z_1995_muni),
                x_2008_muni = gamma*(zbar_2008_muni-z_2008_muni))






# PARAMETER THETA ------------------------------------------------------------------------------------------------------------------------------------

distance_data <-
  read_excel("data/calibration/ipeadata[21-08-2023-01-28].xls")

distance_data$muni_code <- as.numeric(distance_data$muni_code)
# DATA INPUT


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
     + lat+I(lat^2)+distance+log(cattleSlaughter_farmGatePrice_2017), na.action = na.exclude, weights = pasture_area_2017)

# regression results
summary(reg.cattleValueperHa.2017)

# extract fitted values
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(theta = exp(predict(reg.cattleValueperHa.2017, .))/aux.price.2017)


muniTheta.prepData <-
  muniTheta.prepData %>%
    dplyr::filter(zbar_1995_muni >0) %>%
    dplyr::filter(z_1995_muni >0) %>%
    dplyr::filter(!is.na(zbar_1995_muni)) %>% 
    dplyr::filter(!is.na(z_1995_muni))  


# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(muniTheta.prepData,
     file = here::here("data/hmc",
                       "muni_level.Rdata"))

# remove spatial feature
muniTheta.prepData <- muniTheta.prepData %>% sf::st_drop_geometry()

readr::write_csv(muniTheta.prepData,
                 file = here::here("data/hmc", "muni_level.csv"))


# CLEAN TEMP DIR
terra::tmpFiles(current = T, remove = T)
gc()



# END TIMER
tictoc::toc(log = T)



# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
