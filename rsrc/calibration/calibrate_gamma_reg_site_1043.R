library(sf)
library(tictoc)
library(tidyverse)

tictoc::tic(msg = "calibrate_gamma_reg.R script", log = TRUE)

# Load variables at the muni level to calibrate gamma
load("data/calibration/gamma_calibration_1043_sites.Rdata")

# Convert biomass into CO2e, add column of ones, take logs, and scale
df <- calib_df %>%
  mutate(
    X1 = 1,
    log_hist_precip = log(hist_precip),
    log_hist_temp = log(hist_temp),
    latlon = lat * lon,
    log_co2e_ha_2017 = log(co2e)
  ) %>%
  mutate(
    log_hist_precip = scale(log_hist_precip),
    log_hist_temp = scale(log_hist_temp),
    lat = scale(lat),
    lon = scale(lon),
    latlon = scale(latlon),
  ) %>%
  select(
    X1,
    log_hist_precip,
    log_hist_temp,
    lat,
    lon,
    latlon,
    log_co2e_ha_2017,
    id
  )

model_3 <- lm(
  formula = log_co2e_ha_2017 ~
    log_hist_precip +
    log_hist_temp +
    lat +
    lon +
    latlon,
  data = df,
  na.action = na.exclude
)

summary(model_3)


# Output municipality-level regression data
st_write(df,
  "data/calibration/hmc/gamma_reg_site_1043.geojson",
  driver = "GeoJSON",
  delete_dsn = TRUE
)


# END TIMER
tictoc::toc(log = TRUE)
