library(tidyverse)
library(sf)
library(readxl)

# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/calibration/prepData/muniTheta_prepData.Rdata")

# load cattle price series
load("data/calibration/prepData/seriesPriceCattle_prepData.Rdata")

# load farm gate price data
fgp_data <-
  read_excel("data/calibration/farm_gate_price.xlsx")

# load distance data
distance_data <-
  read_excel("data/calibration/ipeadata[21-08-2023-01-28].xls") %>%
  mutate(muni_code = as.numeric(muni_code))


# DATA MANIPULATION

# Drop outliers and separate geometry from df
muni_theta_prep_data <-
  as.data.frame(muniTheta.prepData)[-c(142, 106, 112), ]

geo <- st_geometry(muniTheta.prepData)[-c(142, 106, 112)]

# Merge with distance data and farm price data
muni_theta_prep_data <- muni_theta_prep_data %>%
  left_join(distance_data, by = "muni_code") %>%
  left_join(fgp_data, by = "muni_code") %>%
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
    historical_temp_sq = historical_temp^2,
    lat_sq = lat^2,
    log_cattleSlaughter_valuePerHa_2017 = log(cattleSlaughter_valuePerHa_2017),
    log_cattle_price = log(cattle_price_2017),
    weights = pasture_area_2017
  ) %>%
  mutate(across(
    c(
      historical_precip,
      historical_temp,
      historical_temp_sq,
      lat,
      lat_sq,
      log_cattle_price,
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
    log_cattle_price,
    distance,
    zbar_2017_muni,
    log_cattleSlaughter_valuePerHa_2017,
    weights
  )

# Output municipality-level regression data
st_write(df,
  "data/hmc/muni_data_theta.geojson",
  driver = "GeoJSON",
  delete_dsn = TRUE
)

# Output site-level regression data
for (n in list(10, 24, 40, 78)) {
  id_df <- st_read(sprintf("data/hmc/id_%d.geojson", n))

  # Project data to site level
  site_theta2017 <- df %>%
    st_intersection(id_df)

  # Set area units to hectares
  site_theta2017$muni_site_area <-
    st_area(site_theta2017) %>%
    units::set_units(ha) %>%
    unclass()

  # Write to file
  st_write(
    site_theta2017,
    sprintf("data/hmc/site_%d_data_theta.geojson", n),
    driver = "GeoJSON",
    delete_dsn = TRUE
  )
}
