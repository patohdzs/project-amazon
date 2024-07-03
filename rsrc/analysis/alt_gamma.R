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

# Convert estimates into gammas
new_gamma_rst <- mean_seq_rate_rst * 148.62

# Load calibration data set
load("data/calibration/calibration_1043_sites.Rdata")
grid <- calib_df

# Extract new gamma
new_gamma <- extract(new_gamma_rst, grid, fun = mean, na.rm = TRUE)
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
model_1 <- lm(gamma ~ new_gamma, data = grid)
model_3 <- lm(log(gamma) ~ log(new_gamma), data = grid)
stargazer(model_1, model_3, out = "plots/table_1.tex")
stargazer(model_1, model_3, out = "plots/table_1.txt")


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
