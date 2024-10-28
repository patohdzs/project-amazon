library(sf)
library(terra)
library(stargazer)
library(ggplot2)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

# Load historical temperature and precipitation
temp_rst <- rast("data/clean/temperature.tif")

# Historical percipitation
precip_rst <- rast("data/clean/precipitation.tif")

# Load pixel-level biomass data for 2017
load("data/processed/tree_richness_2017.Rdata")

# Load calibration data set
load("data/calibration/calibration_1043_sites.Rdata")

# Match pixels with sites and calculate average tree richness by site
site_level_tree_richness <- tree_richness_2017 %>%
  st_transform(5880) %>%
  st_join(calib_df) %>%
  select(id, tree_richness_ha) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(
    mean_tree_richness_ha = mean(tree_richness_ha, na.rm = TRUE),
    tree_richness_ha = median(tree_richness_ha, na.rm = TRUE)
  )

# Add site-level tree richness to spatial variables
calib_df <- left_join(calib_df, site_level_tree_richness)

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
calib_df$precip <- extract(
  precip_rst, vect(calib_df),
  fun = mean, na.rm = TRUE
)[, 2]

# Extract average temperature data by muni
calib_df$temp <- extract(
  temp_rst, vect(calib_df),
  fun = mean, na.rm = TRUE
)[, 2]

# Add site centroid coordinates
centroids <- calib_df %>%
  st_centroid() %>%
  st_coordinates()

calib_df$lon <- centroids[, "X"]
calib_df$lat <- centroids[, "Y"]

# Regress log-gamma on geographic covariates
model_1 <- lm(
  formula = mean_tree_richness_ha ~
    lat + lon + (lat * lon) + log(precip) + log(temp),
  data = calib_df,
  na.action = na.exclude
)

model_2 <- lm(
  formula = log(mean_tree_richness_ha) ~
    lat + lon + (lat * lon) + log(precip) + log(temp),
  data = calib_df,
  na.action = na.exclude
)

# Save regression tables
stargazer(model_1, model_2, out = "plots/beta_calib/beta_reg_table.tex")
stargazer(model_1, model_2, out = "plots/beta_calib/beta_reg_table.txt")

# Predict betas
calib_df <- calib_df %>%
  mutate(site_reg_beta = exp(predict(model_2, .)))

fig <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = site_reg_beta)) +
  scale_fill_viridis_c(name = "Number of species / ha")

# Save plot to PDF
ggsave("plots/beta_calib/beta_fitted.pdf", width = 8, height = 6)


# Calibrate growth rate using Rozendaal et al. (2019)
calib_df <- calib_df %>%
  mutate(lambda = 1 - (1 - 0.90)^(1 / 32))
