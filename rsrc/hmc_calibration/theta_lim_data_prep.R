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
aux.price.2017 <-
  seriesPriceCattle.prepData %>%
  dplyr::filter(year == 2017) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(mean_price_2017 = mean(price_real_mon_cattle) / 3.192) %>% # BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
  dplyr::pull(mean_price_2017)



# REGRESSION - CATTLE VALUE (2017)
distance_data <-
  read_excel("data/calibration/ipeadata[21-08-2023-01-28].xls")

distance_data$muni_code <- as.numeric(distance_data$muni_code)

# Remove rows from attribute data

muniTheta.prepData_data <-
  as.data.frame(muniTheta.prepData)  # Convert to regular dataframe
muniTheta.prepData_data <-
  muniTheta.prepData_data[-c(142, 106, 112),]

# Remove geometries
geo_backup <- st_geometry(muniTheta.prepData)
geo_backup <- geo_backup[-c(142, 106, 112)]


predicted_values <-
  read_excel("data/calibration/farm_gate_price.xlsx")

# Combine back into an sf object
muniTheta.prepData <-
  st_sf(muniTheta.prepData_data, geometry = geo_backup)


# Convert to non-spatial dataframe for the merge
muniTheta_no_geo <- as.data.frame(muniTheta.prepData)

# Perform the merge
merged_data <-
  left_join(muniTheta_no_geo, distance_data, by = "muni_code")

# Reattach the geometry
merged_data_sf <- st_sf(merged_data, geometry = geo_backup)


muniTheta.prepData <- merged_data_sf



merged_data <- muniTheta.prepData %>%
  left_join(predicted_values, by = "muni_code") %>%
  mutate(
    cattleSlaughter_farmGatePrice_2017 = ifelse(
      is.na(cattleSlaughter_farmGatePrice_2017),
      average_weighted_price,
      cattleSlaughter_farmGatePrice_2017
    )
  )


muniTheta.prepData <- merged_data


muniTheta.prepData <- muniTheta.prepData %>%
  filter(!is.na(distance))

# Exclude zeros
muniTheta.prepData_filtered <- muniTheta.prepData %>%
  filter(cattleSlaughter_valuePerHa_2017 > 0)





# Create the weight matrix
weights <- muniTheta.prepData_filtered$pasture_area_2017
W <- diag(sqrt(weights))

# Select columns
df <- muniTheta.prepData %>%
  select(
    historical_precip,
    historical_temp,
    lat,
    cattleSlaughter_farmGatePrice_2017,
    distance,
    zbar_2017_muni,
    cattleSlaughter_valuePerHa_2017
  )

# drop spatial feature
df <- cbind(1, df)
df <- df %>%
  mutate(
    I.historical_temp.2. = historical_temp ^ 2,
    I.lat.2. = lat ^ 2,
    log_cattleSlaughter_valuePerHa_2017 = log(cattleSlaughter_valuePerHa_2017)
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
    log_cattleSlaughter_valuePerHa_2017
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




st_write(df,
         "data/hmc/data_theta.geojson",
         driver = "GeoJSON",
         delete_dsn = TRUE)
