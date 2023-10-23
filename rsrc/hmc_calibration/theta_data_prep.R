library(tidyverse)
library(sf)
library(readxl)

# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/calibration/prepData/muniTheta_prepData.Rdata")

# load cattle price series
load("data/calibration/prepData/seriesPriceCattle_prepData.Rdata")

# DATA MANIPULATION
# EXTRACT AVERAGE 2017 PRICE (use real prices because it is normalized to 2017 )
# BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
aux_price_2017 <-
  seriesPriceCattle.prepData %>%
  filter(year == 2017) %>%
  group_by(year) %>%
  summarise(mean_price_2017 = mean(price_real_mon_cattle) / 3.192) %>%
  pull(mean_price_2017)


# REGRESSION - CATTLE VALUE (2017)
distance_data <-
  read_excel("data/calibration/ipeadata[21-08-2023-01-28].xls")

distance_data$muni_code <- as.numeric(distance_data$muni_code)


# Convert to regular dataframe
muni_theta_prepData_data <-
  as.data.frame(muniTheta.prepData)

# Remove rows from attribute data
muni_theta_prepData_data <-
  muni_theta_prepData_data[-c(142, 106, 112),]

# Remove geometries
geo_backup <- st_geometry(muniTheta.prepData)[-c(142, 106, 112)]


predicted_values <-
  read_excel("data/calibration/farm_gate_price.xlsx")

# Combine back into an sf object
muni_theta_prep_data <-
  st_sf(muni_theta_prepData_data, geometry = geo_backup)


# Convert to non-spatial dataframe for the merge
muni_theta_no_geo <- as.data.frame(muni_theta_prep_data)

# Perform the merge
merged_data <-
  left_join(muni_theta_no_geo, distance_data, by = "muni_code")

# Reattach the geometry
merged_data_sf <- st_sf(merged_data, geometry = geo_backup)

muni_theta_prep_data <- merged_data_sf

merged_data <- muni_theta_prep_data %>%
  left_join(predicted_values, by = "muni_code") %>%
  mutate(
    cattleSlaughter_farmGatePrice_2017 = ifelse(
      is.na(cattleSlaughter_farmGatePrice_2017),
      average_weighted_price,
      cattleSlaughter_farmGatePrice_2017
    )
  )

muni_theta_prep_data <- merged_data

muni_theta_prep_data <- muni_theta_prep_data %>%
  filter(!is.na(distance))

# Select columns
df <- muni_theta_prep_data %>%
  select(
    historical_precip,
    historical_temp,
    lat,
    cattleSlaughter_farmGatePrice_2017,
    distance,
    zbar_2017_muni,
    cattleSlaughter_valuePerHa_2017,
    pasture_area_2017
  )

# Add column of 1's, add polynomial terms, and scale
df <- cbind(1, df)
df <- df %>%
  mutate(
    I.historical_temp.2. = historical_temp ^ 2,
    I.lat.2. = lat ^ 2,
    log_cattleSlaughter_valuePerHa_2017 = log(cattleSlaughter_valuePerHa_2017),
    weights = pasture_area_2017
  ) %>%
  select(
    1,
    historical_precip,
    historical_temp,
    I.historical_temp.2.,
    lat,
    I.lat.2.,
    cattleSlaughter_farmGatePrice_2017,
    distance,
    zbar_2017_muni,
    log_cattleSlaughter_valuePerHa_2017,
    weights
  ) %>%
  mutate(X1 = X1) %>%
  mutate(
    historical_precip = scale(historical_precip),
    historical_temp = scale(historical_temp),
    I.historical_temp.2. = scale(I.historical_temp.2.),
    lat = scale(lat),
    I.lat.2. = scale(I.lat.2.),
    cattleSlaughter_farmGatePrice_2017 = scale(cattleSlaughter_farmGatePrice_2017),
    distance = scale(distance)
  )

# Write output
st_write(df,
         "data/hmc/data_theta.geojson",
         driver = "GeoJSON",
         delete_dsn = TRUE)
