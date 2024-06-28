library(stargazer)
library(ggplot2)
library(terra)
library(tidyverse)
library(sf)
library(conflicted)

conflicts_prefer(dplyr::filter())
conflicts_prefer(terra::extract)

# Load carbon sequestration rate estimates
ngs_gextent <- rast("data/raw/global_forest_watch/young_forest_sequestration_rate_Griscom_extent.tif")
ngs <- rast("data/raw/global_forest_watch/sequestration_rate__mean__aboveground__full_extent__Mg_C_ha_yr.tif")

# Load calibration data set
load("data/calibration/hmc/hmc_1043SitesModel.Rdata")
grid <- calibration.1043SitesModel

x <- extract(ngs_gextent, grid, fun = mean, na.rm = TRUE)
y <- extract(ngs, grid, fun = mean, na.rm = TRUE)
grid$new_gamma_gextent <- x[, -1]
grid$new_gamma <- y[, -1]


fig_1 <- grid %>%
  ggplot() +
  geom_sf(aes(fill = new_gamma_gextent), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    title = "Young Forest Sequestration Rate",
    subtitle = "Griscom et al. (2017) reforestable areas",
    fill = "Mg carbon/ha/yr",
    x = "Longitude",
    y = "Latitude"
  )


ggsave(filename = "plots/alt_gamma_fig_1.png", plot = fig_1)

fig_2 <- grid %>%
  ggplot() +
  geom_sf(aes(fill = new_gamma), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    title = "Young Forest Sequestration Rate",
    fill = "Mg carbon/ha/yr",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/alt_gamma_fig_2.png", plot = fig_2)


fig_3 <- grid %>%
  ggplot() +
  geom_sf(aes(fill = gamma_1043Sites), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    title = "Carbon Sequestration Rate",
    fill = "Mg CO2e/ha",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/gamma.png", plot = fig_3)

# Compare old and new gammas
model_1 <- lm(gamma_1043Sites ~ new_gamma, data = grid)
model_2 <- lm(log(gamma_1043Sites) ~ new_gamma, data = grid)
model_3 <- lm(log(gamma_1043Sites) ~ log(new_gamma), data = grid)
stargazer(model_1, model_2, model_3, out = "plots/table_1.tex")
stargazer(model_1, model_2, model_3, out = "plots/table_1.txt")


# Create a scatter plot with a regression line
fig_4 <- ggplot(grid, aes(x = new_gamma, y = gamma_1043Sites)) +
  geom_point() + # Scatter plot
  geom_smooth(method = "lm") + # Regression line
  labs(
    title = "Current gamma v.s new gamma",
    x = "New gamma (Mg carbon/ha/yr)",
    y = "Current gamma (Mg CO2/ha)"
  )

ggsave(filename = "plots/reg_levels.png", plot = fig_4)

fig_5 <- ggplot(grid, aes(x = log(new_gamma), y = log(gamma_1043Sites))) +
  geom_point() + # Scatter plot
  geom_smooth(method = "lm") + # Regression line
  labs(
    title = "Current gamma v.s new gamma",
    x = "log(New gamma)",
    y = "log(Current gamma)"
  )

ggsave(filename = "plots/reg_logs.png", plot = fig_5)


# Scatterplot adjusting to stock
fig_6 <- grid %>%
  mutate(new_gamma = new_gamma * 110.1) %>%
  ggplot(aes(x = new_gamma, y = gamma_1043Sites)) +
  geom_point() +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    title = "Current gamma v.s new gamma",
    x = "New gamma (Mg carbon/ha/yr)",
    y = "Current gamma (Mg CO2/ha)"
  )

ggsave(filename = "plots/scatter_levels.png", plot = fig_6)


# Find % of underestimates
pct_under <- grid %>%
  mutate(new_gamma = new_gamma * 110.1) %>%
  mutate(under_pred = gamma_1043Sites < new_gamma) %>%
  summarize(pct_under = mean(under_pred, na.rm = TRUE) * 100)

# Print the result
print(pct_under)
