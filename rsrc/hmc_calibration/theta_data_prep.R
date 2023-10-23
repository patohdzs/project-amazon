library(tidyverse)
library(sf)
library(readxl)

# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/calibration/prepData/muniTheta_prepData.Rdata")

# load cattle price series
load("data/calibration/prepData/seriesPriceCattle_prepData.Rdata")

# load predicted values
predicted_values <-
  read_excel("data/calibration/farm_gate_price.xlsx")

# load distance data
distance_data <-
  read_excel("data/calibration/ipeadata[21-08-2023-01-28].xls") %>%
  mutate(muni_code = as.numeric(muni_code))


# DATA MANIPULATION

# Drop outliers and separate geometry from df
muni_theta_prep_data <-
  as.data.frame(muniTheta.prepData)[-c(142, 106, 112),]

geo <- st_geometry(muniTheta.prepData)[-c(142, 106, 112)]

# Merge with distance data and farm price data
muni_theta_prep_data <- muni_theta_prep_data %>%
  left_join(distance_data, by = "muni_code") %>%
  left_join(predicted_values, by = "muni_code") %>%
  st_sf(geometry = geo) %>%
  mutate(
    cattle_price_2017 = ifelse(
      is.na(cattleSlaughter_farmGatePrice_2017),
      average_weighted_price,
      cattleSlaughter_farmGatePrice_2017
    )
  ) %>%
  filter(!is.na(distance))

# Add vector of ones, add polynomial terms, and scale
df <- muni_theta_prep_data %>%
  mutate(
    X1 = 1,
    historical_temp_sq = historical_temp ^ 2,
    lat_sq = lat ^ 2,
    log_cattleSlaughter_valuePerHa_2017 = log(cattleSlaughter_valuePerHa_2017),
    weights = pasture_area_2017
  ) %>%
  mutate(across(
    c(
      historical_precip,
      historical_temp,
      historical_temp_sq,
      lat,
      lat_sq,
      cattle_price_2017,
      distance
    ),
    scale
  )) %>%
  select(
    X1,
    historical_precip,
    historical_temp,
    historical_temp_sq,
    lat,
    lat_sq,
    cattle_price_2017,
    distance,
    zbar_2017_muni,
    log_cattleSlaughter_valuePerHa_2017,
    weights
  )

# Write output
st_write(df,
         "data/hmc/data_theta.geojson",
         driver = "GeoJSON",
         delete_dsn = TRUE)
