library(stargazer)
library(ggplot2)
library(terra)
library(tidyverse)
library(sf)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

# Load calibration data set
load("data/calibration/gamma_calibration_1043_sites.Rdata")

# Load carbon sequestration rate estimates
mean_seq_rate_rst <- rast("data/raw/global_forest_watch/sequestration_rate__mean__aboveground__full_extent__Mg_C_ha_yr.tif")

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

# Compute the global min and max for the fill scale
global_min <- min(
  calib_df$co2e,
  calib_df$alt_gamma,
  calib_df$site_reg_gamma,
  na.rm = TRUE
)
global_max <- max(
  calib_df$co2e,
  calib_df$alt_gamma,
  calib_df$site_reg_gamma,
  na.rm = TRUE
)

# Define the fill scale with consistent limits
fill_scale <- scale_fill_viridis_c(
  option = "plasma",
  na.value = "grey",
  limits = c(global_min, global_max)
)


fig_1 <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = co2e), color = "white") +
  fill_scale +
  labs(
    fill = "Mg CO2e/ha",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/gamma_calib/true_gamma.pdf", plot = fig_1)

fig_2 <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = alt_gamma), color = "white") +
  fill_scale +
  labs(
    fill = "Mg CO2e/ha",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/gamma_calib/alt_gamma.pdf", plot = fig_2)


fig_3 <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = site_reg_gamma), color = "white") +
  fill_scale +
  labs(
    fill = "Mg CO2e/ha",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/gamma_calib/1043_sites_reg_gamma.pdf", plot = fig_3)

fig_4 <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = muni_reg_gamma), color = "white") +
  fill_scale +
  labs(
    fill = "Mg CO2e/ha",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/gamma_calib/muni_reg_gamma.pdf", plot = fig_3)

# Compare all estimates v.s data on CO2e
model_1 <- lm(co2e ~ 0 + muni_reg_gamma, data = calib_df)
model_2 <- lm(co2e ~ muni_reg_gamma, data = calib_df)
model_3 <- lm(co2e ~ 0 + site_reg_gamma, data = calib_df)
model_4 <- lm(co2e ~ site_reg_gamma, data = calib_df)
model_5 <- lm(co2e ~ 0 + alt_gamma, data = calib_df)
model_6 <- lm(co2e ~ alt_gamma, data = calib_df)

stargazer(model_1, model_2, model_3, model_4, model_5, model_6, out = "plots/gamma_calib/gamma_comparison.tex")
stargazer(model_1, model_2, model_3, model_4, model_5, model_6, out = "plots/gamma_calib/gamma_comparison.txt")

# Compare our estimates v.s Cook-patton
model_1 <- lm(alt_gamma ~ 0 + muni_reg_gamma, data = calib_df)
model_2 <- lm(alt_gamma ~ muni_reg_gamma, data = calib_df)
model_3 <- lm(log(alt_gamma) ~ log(muni_reg_gamma), data = calib_df)

stargazer(model_1, model_2, model_3, out = "plots/muni_reg_gamma_vs_cook_patton.tex")
stargazer(model_1, model_2, model_3, out = "plots/muni_reg_gamma_vs_cook_patton.txt")

# Compare our estimates v.s Cook-patton
model_1 <- lm(alt_gamma ~ 0 + site_reg_gamma, data = calib_df)
model_2 <- lm(alt_gamma ~ site_reg_gamma, data = calib_df)
model_3 <- lm(log(alt_gamma) ~ log(site_reg_gamma), data = calib_df)

stargazer(model_1, model_2, model_3, out = "plots/site_reg_gamma_vs_cook_patton.tex")
stargazer(model_1, model_2, model_3, out = "plots/site_reg_gamma_vs_cook_patton.txt")

# Scatterplot adjusting to stock
fig_5 <- calib_df %>%
  ggplot(aes(x = alt_gamma, y = site_reg_gamma)) +
  geom_point() +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    x = "Cook-Patton et al. (2020) gamma (Mg CO2e/ha)",
    y = "1043-sites regression gamma (Mg CO2e/ha)"
  )

ggsave(filename = "plots/gamma_calib/scatter_site_reg.pdf", plot = fig_5)

fig_6 <- calib_df %>%
  ggplot(aes(x = alt_gamma, y = muni_reg_gamma)) +
  geom_point() +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    x = "Cook-Patton et al. (2020) gamma (Mg CO2e/ha)",
    y = "Municipal regression gamma (Mg CO2e/ha)"
  )

ggsave(filename = "plots/gamma_calib/scatter_muni_reg.pdf", plot = fig_6)


# Find % of underestimates
pct_under <- calib_df %>%
  mutate(under_pred = gamma < alt_gamma)

# Print the result
print(mean(pct_under$under_pred, na.rm = TRUE) * 100)
