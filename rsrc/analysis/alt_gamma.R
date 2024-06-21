library(stargazer)
library(ggplot2)
library(terra)
library(sf)

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
    fill = "Mg carbon/ha",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/gamma.png", plot = fig_3)

# Compare old and new gammas
model_1 <- lm(gamma_1043Sites ~ new_gamma, data = grid)
model_2 <- lm(scale(gamma_1043Sites) ~ scale(new_gamma), data = grid)
stargazer(model_1, model_2, out = "plots/table_1.tex")
stargazer(model_1, model_2, out = "plots/table_1.txt")
