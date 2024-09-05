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

# Load raster files
pq_rst <- rast(
  list.files(
    "data/raw/mapbiomas/pasture_quality/",
    pattern = "2015",
    full.names = TRUE
  )
)

sec_veg_age_rst <- rast(
  list.files(
    "data/raw/mapbiomas/secondary_vegetation_age/",
    pattern = "2017",
    full.names = TRUE
  )
)

clean_mapbiomas <- rast("data/clean/land_use_cover_2000.tif")

# Historical temperature
temp_rst <- rast("data/clean/temperature.tif")

# Historical percipitation
precip_rst <- rast("data/clean/precipitation.tif")

# Load calibration data set
load("data/calibration/calibration_1043_sites.Rdata")

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
calib_df$hist_precip <- extract(
  precip_rst,
  vect(calib_df),
  fun = mean,
  na.rm = TRUE
)[, 2]

# Extract average temperature data by muni
calib_df$hist_temp <- extract(
  temp_rst,
  vect(calib_df),
  fun = mean,
  na.rm = TRUE
)[, 2]

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

ggsave(filename = "plots/adj_costs_calib/1043_sites_share_low_pq_fitted.pdf", plot = fig_1)


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

stargazer(model_1, out = "plots/adj_costs_calib/1043_sites_low_pq_share_tobit.tex")
stargazer(model_1, out = "plots/adj_costs_calib/1043_sites_low_pq_share_tobit.txt")

# Save calibration data CSV
calib_df %>%
  st_drop_geometry() %>%
  write_csv(file = "data/calibration/calibration_1043_sites.csv")

# Select Amazonia subset
sec_veg_age_rst <- sec_veg_age_rst |>
  crop(clean_mapbiomas) |>
  mask(clean_mapbiomas)

pq_rst <- pq_rst |>
  crop(clean_mapbiomas) |>
  mask(clean_mapbiomas)

no_pq <- (pq_rst == 0) & (sec_veg_age_rst == 2)

# Set areas with no two-year-old sec veg to 0
pq_rst[(sec_veg_age_rst != 2)] <- 0

# Compute pixel areas
pixel_areas <- cellSize(pq_rst, unit = "ha")

# Calculate the total area for each pasture quality
ref <- zonal(pixel_areas, pq_rst, sum, na.rm = TRUE) %>%
  slice(-1) %>%
  mutate(prop = area / sum(area, na.rm = TRUE))

print(ref)

no_pq_ref <- zonal(pixel_areas, no_pq, sum, na.rm = TRUE) %>%
  slice(-1) %>%
  mutate(prop = area / sum(area, na.rm = TRUE))

print(no_pq_ref)

no_pq_ref <- no_pq_ref[[1, "area"]]

# Get mean share of low pasture quality across sites
avg_share_low_pq <- mean(calib_df$fit_share_low_pq)

# Get site-specific MC's (up to factor of zeta)
r_low <- ref[1, 2] + avg_share_low_pq * no_pq_ref
r_high <- ref[2, 2] + ref[3, 2] + (1 - avg_share_low_pq) * no_pq_ref

mc_low_pq <- r_low * calib_df$fit_share_low_pq
mc_high_pq <- r_high * (1 - calib_df$fit_share_low_pq)

# Compute implied adjustment cost params (MC = MR condition)
zeta_2 <- 674 / mean(mc_low_pq)
zeta_3 <- 52 / mean(mc_high_pq)

print("Zeta (MC)")
print(zeta_2)
print(zeta_3)

# Deforestation cost calibration
# We transform to dollars using an FX rate of 4.14 (December 2019)
deforestation_cost <- 1614.54 / 4.14

# The forest area in 2008 represented 72% of the Legal Amazon area,
#   which covers 501,506,775 ha, so the transition in hectares is
#   0.065*0.72*501,506,775 = 23,470,517,
#   resulting in an annual average of 2,347,052 ha.
transition_area <- (0.065 * 0.72 * 501506775) / (2017 - 2008 + 1)
zeta_1 <- deforestation_cost / transition_area
print("Zeta (deforestation)")
print(zeta_1)

# Alternative value based on https://shorturl.at/jgp9i
zeta_alt <- 483 / transition_area
