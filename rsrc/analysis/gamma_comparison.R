library(stargazer)
library(ggplot2)
library(terra)
library(tidyverse)
library(sf)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

# Load pixel-level biomass data for 2017
load("data/processed/pixel_biomass_2017.Rdata")

# Load calibration data set
load("data/calibration/calibration_1043_sites.Rdata")

# Load carbon sequestration rate estimates
mean_seq_rate_rst <- rast("data/raw/global_forest_watch/sequestration_rate__mean__aboveground__full_extent__Mg_C_ha_yr.tif")

# Convert and filter biomass -> CO2e
pixel_biomass_2017 <- pixel_biomass_2017 %>%
  st_transform(st_crs(calib_df)) %>%
  mutate(co2e = (agb_2017 / 2) * (44 / 12)) %>%
  filter(co2e > 0, !is.na(co2e))

# Match pixels with sites
true_gamma <- pixel_biomass_2017 %>%
  st_join(calib_df) %>%
  select(id, co2e) %>%
  st_drop_geometry()

# Calculate average carbon density on primary forest areas by site
true_gamma <- true_gamma %>%
  group_by(id) %>%
  summarise(true_gamma = mean(co2e, na.rm = TRUE))

# Add true_gamma to spatial variables
calib_df <- left_join(calib_df, true_gamma)

# Convert estimates into gammas
alt_gamma_rst <- mean_seq_rate_rst * 148.62

# Extract new gamma
alt_gamma <- extract(
  alt_gamma_rst, calib_df,
  fun = mean,
  weights = TRUE,
  na.rm = TRUE
)
calib_df$alt_gamma <- alt_gamma[, -1]

# Transform CRS for better plots
calib_df <- calib_df %>% st_transform(crs = 4326)

fig_1 <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = true_gamma), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    fill = "Mg CO2e/ha",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/gamma_calib/true_gamma.pdf", plot = fig_1)

fig_2 <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = alt_gamma), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    fill = "Mg CO2e/ha",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/gamma_calib/alt_gamma.pdf", plot = fig_2)


fig_3 <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = gamma), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    fill = "Mg CO2e/ha",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/gamma_calib/gamma.pdf", plot = fig_3)

# Compare our estimates v.s and alt gammas
model_1 <- lm(true_gamma ~ 0 + gamma, data = calib_df)
model_2 <- lm(true_gamma ~ gamma, data = calib_df)
model_3 <- lm(true_gamma ~ 0 + alt_gamma, data = calib_df)
model_4 <- lm(true_gamma ~ alt_gamma, data = calib_df)
stargazer(model_1, model_2, model_3, model_4, out = "plots/gamma_calib/table_1.tex")
stargazer(model_1, model_2, model_3, model_4, out = "plots/gamma_calib/table_1.txt")
