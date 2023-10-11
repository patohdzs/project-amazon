library(tidyverse)
library(sf)

# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/calibration/muniTheta_prepData_gamma.Rdata")

muniTheta.prepData<-muniTheta.prepData %>%
  dplyr::mutate(co2e_ha_2017 = (agb_2017/2)*(44/12))


reg.gamma.2017 <-
  muniTheta.prepData  %>%
  lm(formula = log(co2e_ha_2017)  ~ log(historical_precip) + log(historical_temp) +log(lat)+log(lon), na.action = na.exclude)

reg_summary <- summary(reg.gamma.2017)
regressor_df <- as.data.frame(reg.gamma.2017$model[-1])


new_df <- muniTheta.prepData %>%
  select(historical_precip, historical_temp, lat, lon, co2e_ha_2017) %>%
  mutate(log_historical_precip = log(historical_precip), log_co2e_ha_2017 = log(co2e_ha_2017)) %>%
  mutate(log_historical_temp = log(historical_temp)) %>%
  mutate(log_lat = log(lat)) %>%
  mutate(log_lon = log(lon))  %>%
  mutate(log_historical_precip = (log_historical_precip-mean(regressor_df$`log(historical_precip)`))/ sd(regressor_df$`log(historical_precip)`)) %>%
  mutate(log_historical_temp = (log_historical_temp-mean(regressor_df$`log(historical_temp)`))/ sd(regressor_df$`log(historical_temp)`)) %>%
  mutate(log_lat = (log_lat-mean(regressor_df$`log(lat)`))/ sd(regressor_df$`log(lat)`)) %>%
  mutate(log_lon = (log_lon-mean(regressor_df$`log(lon)`))/ sd(regressor_df$`log(lon)`))


new_df <- cbind(1, new_df)
new_df <- new_df %>%
  select(1, log_historical_precip,log_historical_temp,log_lat, log_lon, log_co2e_ha_2017)


st_write(new_df, "data/hmc_new/data_gamma.geojson", driver = "GeoJSON", delete_dsn = TRUE)
