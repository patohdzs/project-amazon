library(sf)
library(glue)
library(tictoc)
library(terra)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter)

# Load calibration data
load("data/calibration/calibration_1043_sites.Rdata")

# Load raster files
pq_rst <- rast(
  list.files(
    "data/raw/mapbiomas/pasture_quality/",
    pattern = "2014",
    full.names = TRUE
  )
)

sec_veg_age_rst <- rast(
  list.files(
    "data/raw/mapbiomas/secondary_vegetation_age/",
    pattern = "2016",
    full.names = TRUE
  )
)

# Select Amazonia subset
sec_veg_age_rst <- sec_veg_age_rst |>
  crop(pq_rst) |>
  mask(pq_rst)

# Set areas with no two-year-old sec veg to 0
pq_rst[(sec_veg_age_rst != 2)] <- 0

# Compute pixel areas
pixel_areas <- cellSize(pq_rst, unit = "ha")

# Calculate the total area for each pasture quality
reforested_areas <- zonal(pixel_areas, pq_rst, sum, na.rm = TRUE) %>%
  slice(-1) %>%
  mutate(prop = area / sum(area, na.rm = TRUE))

print(reforested_areas)

# Get site-specific MC's
MC <- reforested_areas[1, 2] * calib_df$share_low_pq

# Compute implied adjustment cost params (MC = MR condition)
zeta_1_mc <- 674 / reforested_areas[1, 2]
zeta_2_mc <- 472 / reforested_areas[2, 2]
zeta_3_mc <- 52 / reforested_areas[3, 2]
print("Zeta (MC)")
print(zeta_1_mc)
print(zeta_2_mc)
print(zeta_3_mc)


# Old calibration
# We transform to dollars using an FX rate of 4.14 (December 2019)
aux_transition_cost <- 1614.54 / 4.14

# The forest area in 2008 represented 72% of the Legal Amazon area,
#   which covers 501,506,775 ha, so the transition in hectares is
#   0.065*0.72*501,506,775 = 23,470,517,
#   resulting in an annual average of 2,347,052 ha.
aux_transition_area <- (0.065 * 0.72 * 501506775) / (2017 - 2008 + 1)
zeta_old <- aux_transition_cost / aux_transition_area
print("Cost one hectare (old)")
print(zeta_old / 2)

# Alternative value based on https://shorturl.at/jgp9i
zeta_alt <- 483 / aux_transition_area


# Plot share of type 1 pasture
fig_2 <- calib_df %>%
  st_transform(4326) %>%
  ggplot() +
  geom_sf(aes(fill = share_pq_1_pasture * 100), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    fill = "% of total pasture area",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/adj_costs_calib/share_pq_1_pasture.pdf", plot = fig_2)

# Plot share of type 2 pasture
fig_3 <- calib_df %>%
  st_transform(4326) %>%
  ggplot() +
  geom_sf(aes(fill = share_pq_2_pasture * 100), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    fill = "% of total pasture area",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/adj_costs_calib/share_pq_2_pasture.pdf", plot = fig_3)



# Plot share of type 3 pasture
fig_4 <- calib_df %>%
  st_transform(4326) %>%
  ggplot() +
  geom_sf(aes(fill = share_pq_3_pasture * 100), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    fill = "% of total pasture area",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/adj_costs_calib/share_pq_3_pasture.pdf", plot = fig_4)


calib_df <- calib_df %>%
  select(starts_with("share_pq")) %>%
  write_csv("data/calibration/pq_shares.csv")
