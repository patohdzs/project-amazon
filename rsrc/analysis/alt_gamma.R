library(stargazer)
library(ggplot2)
library(terra)
library(tidyverse)
library(sf)
library(conflicted)

conflicts_prefer(dplyr::filter())
conflicts_prefer(terra::extract)

# Load carbon sequestration rate estimates
mean_seq_rate_rst <- rast("data/raw/global_forest_watch/sequestration_rate__mean__aboveground__full_extent__Mg_C_ha_yr.tif")

# Load calibration data set
load("data/calibration/calibration_1043_sites.Rdata")

# Load municipal data
load("data/processed/muni_data.Rdata")

# Convert estimates into gammas
new_gamma_rst <- mean_seq_rate_rst * 148.62

grid <- calib_df

# Extract new gamma
new_gamma <- extract(
  new_gamma_rst, grid,
  fun = mean,
  weights = TRUE,
  na.rm = TRUE
)
grid$new_gamma <- new_gamma[, -1]

# Transform CRS for better plots
grid <- grid %>% st_transform(crs = 4326)

fig_2 <- grid %>%
  ggplot() +
  geom_sf(aes(fill = new_gamma), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    fill = "Mg CO2e/ha",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/alt_gamma.pdf", plot = fig_2)


fig_3 <- grid %>%
  ggplot() +
  geom_sf(aes(fill = gamma), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    fill = "Mg CO2e/ha",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/gamma.pdf", plot = fig_3)

# Compare old and new gammas
model_1 <- lm(gamma ~ 0 + new_gamma, data = grid)
model_2 <- lm(gamma ~ new_gamma, data = grid)
model_3 <- lm(log(gamma) ~ log(new_gamma), data = grid)
stargazer(model_1, model_2, model_3, out = "plots/table_1.tex")
stargazer(model_1, model_2, model_3, out = "plots/table_1.txt")


# Create a scatter plot with a regression line
fig_4 <- ggplot(grid, aes(x = new_gamma, y = gamma)) +
  geom_point() + # Scatter plot
  geom_smooth(method = "lm") + # Regression line
  labs(
    x = "Alternative gamma (Mg CO2/ha)",
    y = "Current gamma (Mg CO2/ha)"
  )

ggsave(filename = "plots/reg_levels.pdf", plot = fig_4)

fig_5 <- ggplot(grid, aes(x = log(new_gamma), y = log(gamma))) +
  geom_point() + # Scatter plot
  geom_smooth(method = "lm") + # Regression line
  labs(
    x = "log(Alternative gamma)",
    y = "log(Current gamma)"
  )

ggsave(filename = "plots/reg_logs.pdf", plot = fig_5)


# Scatterplot adjusting to stock
fig_6 <- grid %>%
  mutate(new_gamma = new_gamma) %>%
  ggplot(aes(x = new_gamma, y = gamma)) +
  geom_point() +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    x = "Alternative gamma (Mg CO2e/ha)",
    y = "Current gamma (Mg CO2e/ha)"
  )

ggsave(filename = "plots/scatter_levels.pdf", plot = fig_6)


# Find % of underestimates
pct_under <- grid %>%
  mutate(under_pred = gamma < new_gamma)

# Print the result
print(mean(pct_under$under_pred, na.rm = TRUE) * 100)

# Get alternative gammas at municipality level
new_gamma_muni <- extract(
  new_gamma_rst, muni_data,
  fun = mean,
  weights = TRUE,
  na.rm = TRUE
)

muni_data$new_gamma <- new_gamma_muni[, -1]

# Convert biomass to CO2 equivalent
muni_data <- muni_data %>%
  mutate(co2e_ha_2017 = (agb_2017 / 2) * (44 / 12))

# Regress log-gamma on geographic covariates
current_gamma_reg <- lm(
  formula = log(co2e_ha_2017) ~
    log(historical_precip) +
    log(historical_temp) +
    log(lat) +
    log(lon),
  data = muni_data,
  na.action = na.exclude
)


alt_gamma_reg <- lm(
  formula = log(new_gamma) ~
    log(historical_precip) +
    log(historical_temp) +
    log(lat) +
    log(lon),
  data = muni_data,
  na.action = na.exclude
)

stargazer(current_gamma_reg, alt_gamma_reg, out = "plots/table_2.tex")
stargazer(current_gamma_reg, alt_gamma_reg, out = "plots/table_2.txt")
