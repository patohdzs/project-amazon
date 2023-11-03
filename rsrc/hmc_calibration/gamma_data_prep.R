library(tidyverse)
library(sf)

# DATA INPUT
load("data/calibration/prepData/muniTheta_prepData_gamma.Rdata")

# Convert biomass into CO2e, add column of ones, take logs, and scale
df <- muniTheta.prepData %>%
  mutate(co2e_ha_2017 = (agb_2017 / 2) * (44 / 12)) %>%
  mutate(
    X1 = 1,
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
  ) %>%
  select(X1,
         log_historical_precip,
         log_historical_temp,
         log_lat,
         log_lon,
         log_co2e_ha_2017)

st_write(df,
         "data/hmc/muni_data_gamma.geojson",
         driver = "GeoJSON",
         delete_dsn = TRUE)


for (n in list(10, 24, 40, 78)) {
  id_df <- st_read(sprintf("data/hmc/id_%d.geojson", n))

  site_gamma2017 <- df %>%
    st_intersection(id_df)

  site_gamma2017$muni_site_area <-
    st_area(site_gamma2017) %>%
    units::set_units(ha) %>%
    unclass()

  st_write(
    site_gamma2017,
    sprintf("data/hmc/site_%d_data_gamma.geojson", n),
    driver = "GeoJSON",
    delete_dsn = TRUE
  )

}
