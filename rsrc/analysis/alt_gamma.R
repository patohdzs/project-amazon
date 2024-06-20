library(ggplot2)
library(terra)
library(sf)

# Load carbon sequestration rate estimates
ngs_gextent <- rast("data/raw/global_forest_watch/young_forest_sequestration_rate_Griscom_extent.tif")
ngs <- rast("data/raw/global_forest_watch/sequestration_rate__mean__aboveground__full_extent__Mg_C_ha_yr.tif")

# Load amazon biome grid
biome_shares <- rast("data/processed/amazon_biome_1043_sites.tif")

grid <- biome_shares %>%
  as.polygons(dissolve = FALSE) %>%
  st_as_sf()

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
