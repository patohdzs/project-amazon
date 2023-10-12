library(tidyverse)

# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/calibration/muniTheta_prepData.Rdata")

# load cattle price series
load("data/calibration/seriesPriceCattle_prepData.Rdata")


# DATA MANIPULATION

# EXTRACT AVERAGE 2017 PRICE (use real prices because it is normalized to 2017 )
aux.price.2017 <-
  seriesPriceCattle.prepData %>%
  dplyr::filter(year == 2017) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(mean_price_2017 = mean(price_real_mon_cattle)/3.192) %>% # BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
  dplyr::pull(mean_price_2017)


# REGRESSION - CATTLE VALUE (2017)


my_data <- read_excel("data/calibration/ipeadata[21-08-2023-01-28].xls")
my_data$muni_code <- as.numeric(my_data$muni_code)

# Remove rows from attribute data
a<-muniTheta.prepData
muniTheta.prepData_data <- as.data.frame(muniTheta.prepData)  # Convert to regular dataframe
muniTheta.prepData_data <- muniTheta.prepData_data[-c(142, 106, 112), ]

# Remove geometries
geo_backup <- st_geometry(muniTheta.prepData)
geo_backup <- geo_backup[-c(142, 106, 112)]


predicted_values <- read_excel("data/calibration/farm_gate_price.xlsx")

# Combine back into an sf object
muniTheta.prepData <- st_sf(muniTheta.prepData_data, geometry = geo_backup)

# 2. Merging the cleaned muniTheta.prepData with my_data

# Convert to non-spatial dataframe for the merge
muniTheta_no_geo <- as.data.frame(muniTheta.prepData)

# Perform the merge
merged_data <- left_join(muniTheta_no_geo, my_data, by = "muni_code")

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
#muniTheta.prepData_filtered <- na.omit(muniTheta.prepData_filtered)


b<-muniTheta.prepData
c<-muniTheta.prepData_filtered

new_df <- muniTheta.prepData %>%
  select(historical_precip,
         historical_temp, lat,cattleSlaughter_farmGatePrice_2017,distance,zbar_2017_muni)

# drop spatial feature
new_df <- new_df %>%
  mutate(`I(historical_temp^2)` = historical_temp^2)
new_df <- new_df %>%
  mutate(`I(lat^2)` = lat^2)
new_df <- cbind(1, new_df)
new_df <- new_df %>%
  select(1, historical_precip,historical_temp,I.historical_temp.2.,lat,I.lat.2.,cattleSlaughter_farmGatePrice_2017,distance,zbar_2017_muni)

new_df <- new_df %>%
  mutate(X1=X1) %>%
  mutate(historical_precip = historical_precip/sqrt(mean(muniTheta.prepData_filtered$historical_precip^2))) %>%
  mutate(historical_temp = historical_temp/sqrt(mean(muniTheta.prepData_filtered$historical_temp^2))) %>%
  mutate(I.historical_temp.2. = I.historical_temp.2./sqrt(mean(muniTheta.prepData_filtered$historical_temp^4))) %>%
  mutate(lat = lat/sqrt(mean(muniTheta.prepData_filtered$lat^2))) %>%
  mutate(I.lat.2. =I.lat.2./sqrt(mean(muniTheta.prepData_filtered$lat^4))) %>%
  #mutate(cattleSlaughter_farmGatePrice_2017=cattleSlaughter_farmGatePrice_2017/35.75924071280666)%>%
  mutate(cattleSlaughter_farmGatePrice_2017=cattleSlaughter_farmGatePrice_2017/sqrt(mean(muniTheta.prepData_filtered$cattleSlaughter_farmGatePrice_2017^2)))%>%
  mutate(distance = distance/sqrt(mean(muniTheta.prepData_filtered$distance^2)))




st_write(new_df, "data/hmc/data_theta.geojson", driver = "GeoJSON",delete_dsn = TRUE)
