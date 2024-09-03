library(stargazer)
library(ggplot2)
library(terra)
library(tidyverse)
library(sf)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

# Historical temperature
temp_rst <- rast("data/clean/temperature.tif")

# Historical percipitation
precip_rst <- rast("data/clean/precipitation.tif")

# Load pixel-level biomass data for 2017
load("data/processed/pixel_biomass_2017.Rdata")

# Load calibration data set
load("data/calibration/calibration_78_sites.Rdata")

# Convert and filter biomass -> CO2e
pixel_biomass_2017 <- pixel_biomass_2017 %>%
  st_transform(st_crs(calib_df)) %>%
  mutate(co2e = (agb_2017 / 2) * (44 / 12)) %>%
  filter(co2e > 0, !is.na(co2e))

# Match pixels with sites and calculate average CO2e by site
site_level_co2e <- pixel_biomass_2017 %>%
  st_join(calib_df) %>%
  select(id, co2e) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(co2e = mean(co2e, na.rm = TRUE))

# Add site-level CO2e to spatial variables
calib_df <- left_join(calib_df, site_level_co2e)

# Match spatial sample crs with raster
calib_df <- calib_df %>%
  st_transform(st_crs(precip_rst))

# Crop rasters
precip_rst <- crop(precip_rst, calib_df)
temp_rst <- crop(temp_rst, calib_df)

# Calculate average yearly precipitation and temperature
precip_rst <- mean(precip_rst)
temp_rst <- mean(temp_rst)

# Extract total precipitation data by muni
calib_df$hist_precip <- extract(precip_rst, vect(calib_df), fun = mean, na.rm = TRUE)[, 2]

# Extract average temperature data by muni
calib_df$hist_temp <- extract(temp_rst, vect(calib_df), fun = mean, na.rm = TRUE)[, 2]


# Add site centroid coordinates
centroids <- calib_df %>%
  st_centroid() %>%
  st_coordinates()

calib_df$lon <- centroids[, "X"]
calib_df$lat <- centroids[, "Y"]

# Regress log-gamma on geographic covariates
model_1 <- lm(
  formula = log(co2e) ~
    log(hist_precip) +
    log(hist_temp) +
    asinh(lat) +
    asinh(lon),
  data = calib_df,
  na.action = na.exclude
)

model_2 <- lm(
  formula = log(co2e) ~
    log(hist_precip) +
    log(hist_temp) +
    lat +
    lon,
  data = calib_df,
  na.action = na.exclude
)

model_3 <- lm(
  formula = log(co2e) ~
    log(hist_precip) +
    log(hist_temp) +
    lat +
    lon +
    lat * lon,
  data = calib_df,
  na.action = na.exclude
)
summary(model_3)

# Save regression tables
#stargazer(model_1, model_2, model_3, out = "plots/gamma_calib/1043_sites_reg_table.tex")
#stargazer(model_1, model_2, model_3, out = "plots/gamma_calib/1043_sites_reg_table.txt")


# Predict gammas
calib_df <- calib_df %>%
  mutate(muni_reg_gamma = gamma) %>%
  mutate(site_reg_gamma = exp(predict(model_3, .)))





save(calib_df, file = "data/calibration/gamma_calibration_78_sites.Rdata")


df <- calib_df %>%
  mutate(
    X1 = 1,
    log_hist_precip = log(hist_precip),
    log_hist_temp = log(hist_temp),
    latlon = lat * lon,
    log_co2e_ha_2017 = log(co2e)
  ) %>%
  select(
    X1,
    log_hist_precip,
    log_hist_temp,
    lat,
    lon,
    latlon,
    log_co2e_ha_2017
  ) 

# Output municipality-level regression data
st_write(df,
         "data/calibration/hmc/gamma_data_site_78.geojson",
         driver = "GeoJSON",
         delete_dsn = TRUE
)


