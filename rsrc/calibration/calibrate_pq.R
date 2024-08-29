# TODO: need to combine with calib script

library(AER)
library(ggplot2)
library(terra)
library(tidyverse)
library(sf)
library(conflicted)
library(stargazer)

conflicts_prefer(dplyr::filter)
conflicts_prefer(AER::tobit)
conflicts_prefer(terra::extract)

# Historical temperature
temp_rst <- rast("data/clean/temperature.tif")

# Historical percipitation
precip_rst <- rast("data/clean/precipitation.tif")

# Load calibration data set
load("data/calibration/calibration_78_sites.Rdata")

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

# Regress share of low quality pasture on geographic covariates
model_1 <- tobit(
  share_low_pq ~
    log(hist_precip) +
    log(hist_temp) +
    lat +
    lon +
    lat * lon, ,
  left = 0,
  right = 1,
  data = calib_df
)

# Get latent variable location and scale
mu <- predict(model_1, calib_df)
sigma <- model_1$scale

# Compute Pr(y>0)
prob_0 <- pnorm(mu / sigma)

# Compute conditional expectation
lambda <- function(x) dnorm(x) / pnorm(x)
exp_0 <- mu + sigma * lambda(mu / sigma)

# Compute unconditional expectation
calib_df$fit_share_low_pq <- prob_0 * exp_0

# Define the fill scale with consistent limits
fill_scale <- scale_fill_viridis_c(
  option = "plasma",
  na.value = "grey",
  limits = c(0, 100)
)

fig_1 <- calib_df %>%
  st_transform(4326) %>%
  ggplot() +
  geom_sf(aes(fill = fit_share_low_pq * 100), color = "white") +
  fill_scale +
  labs(
    fill = "Share of low quality pasture (%)",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/adj_costs_calib/78_sites_share_low_pq_fitted.pdf", plot = fig_1)


fig_2 <- calib_df %>%
  st_transform(4326) %>%
  ggplot() +
  geom_sf(aes(fill = share_low_pq * 100), color = "white") +
  fill_scale +
  labs(
    fill = "Share of low quality pasture (%)",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/adj_costs_calib/1043_sites_share_low_pq.pdf", plot = fig_2)

stargazer(model_1, out = "plots/adj_costs_calib/78_sites_low_pq_share_tobit.tex")
stargazer(model_1, out = "plots/adj_costs_calib/78_sites_low_pq_share_tobit.txt")

# Save calibration data CSV
calib_df %>%
  st_drop_geometry() %>%
  write_csv(file = "data/calibration/calibration_78_sites.csv")
