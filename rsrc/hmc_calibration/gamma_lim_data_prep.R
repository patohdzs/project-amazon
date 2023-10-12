library(tidyverse)
library(sf)

# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/calibration/muniTheta_prepData_gamma.Rdata")

df <- muniTheta.prepData %>%
  mutate(co2e_ha_2017 = (agb_2017 / 2) * (44 / 12)) %>%
  select(historical_precip, historical_temp, lat, lon, co2e_ha_2017) %>%
  mutate(

    log_historical_precip = log(historical_precip),
    log_historical_temp = log(historical_temp),
    log_lat = log(lat),
    log_lon = log(lon),
    log_co2e_ha_2017 = log(co2e_ha_2017)
  ) %>%
  mutate(
    log_historical_precip = scale(log_historical_precip),
    log_historical_temp = scale(log_historical_temp),
    log_lat = scale(log_lat),
    log_lon = scale(log_lon)
  )



reg_gamma_2017 <-
  df  %>%
  lm(
    formula = log_co2e_ha_2017  ~ log_historical_precip + log_historical_temp + log_lat + log_lon,
    na.action = na.exclude
  )

reg_summary <- summary(reg_gamma_2017)
regressor_df <- as.data.frame(reg_gamma_2017$model[-1])

df <- cbind(1, df)
df <- df %>%
  select(1, log_historical_precip,log_historical_temp,log_lat, log_lon, log_co2e_ha_2017)


st_write(df, "data/hmc/data_gamma.geojson", driver = "GeoJSON", delete_dsn = TRUE)
