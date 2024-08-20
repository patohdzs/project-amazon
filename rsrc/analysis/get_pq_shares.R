library(sf)
library(glue)
library(tictoc)
library(terra)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter)

# Load raster data (amazon biome share, pixel areas, and land uses)
amazon_rsts <- rast(
  list.files(
    "data/processed/",
    pattern = ".tif",
    full.names = TRUE
  )
)

# Extract variables as polygons and project data
calib_df <- as.polygons(amazon_rsts, dissolve = FALSE) %>%
  st_as_sf() %>%
  st_transform(5880)

# Remove sites with less than 3% overlap with the amazon biome
calib_df <- calib_df %>%
  filter(share_amazon_biome >= 0.03)

# Add id variable
calib_df$id <- seq_len(nrow(calib_df))

# Transform share variables into area (ha)
calib_df <- calib_df %>%
  mutate(across(
    starts_with("share_"),
    .names = "area_{.col}",
    ~ . * pixel_area_ha
  )) %>%
  rename_with(
    ~ str_replace(., "share_", ""),
    starts_with("area_")
  )

# Compute total pasture area
calib_df <- calib_df %>%
  mutate(
    area_pasture_2017 = area_pasture_quality_1_2017 +
      area_pasture_quality_2_2017 +
      area_pasture_quality_3_2017
  )

# Compute pasture area as share of agricultural area (Z)
calib_df <- calib_df %>%
  mutate(
    share_pasture_z = area_pasture_2017 / area_agricultural_use_2017
  )

# Compute as share of total pasture area
calib_df <- calib_df %>%
  mutate(
    share_pq_1_pasture = area_pasture_quality_1_2017 / area_pasture_2017,
    share_pq_2_pasture = area_pasture_quality_2_2017 / area_pasture_2017,
    share_pq_3_pasture = area_pasture_quality_3_2017 / area_pasture_2017,
  )

# Compute as share of total agricultural area (Z)
calib_df <- calib_df %>%
  mutate(
    share_pq_1_z = area_pasture_quality_1_2017 / area_agricultural_use_2017,
    share_pq_2_z = area_pasture_quality_2_2017 / area_agricultural_use_2017,
    share_pq_3_z = area_pasture_quality_3_2017 / area_agricultural_use_2017,
  )

# Some checks
calib_df %>%
  pull(share_pasture_z) %>%
  summary()

calib_df %>%
  filter(share_pasture_z > 1) %>%
  nrow()

calib_df %>%
  select(
    share_pq_1_pasture,
    share_pq_2_pasture,
    share_pq_3_pasture
  ) %>%
  summary()

calib_df %>%
  select(
    share_pq_1_z,
    share_pq_2_z,
    share_pq_3_z
  ) %>%
  summary()

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

# Compute implied adjustment cost params
zeta_1 <- 674 * 2 / reforested_areas[1, 2]
zeta_2 <- 472 * 2 / reforested_areas[2, 2]
zeta_3 <- 52 * 2 / reforested_areas[3, 2]

print(zeta_1)
print(zeta_2)
print(zeta_3)


# Compute implied adjustment cost params (MC = MR condition)
zeta_1_mc <- 674 / reforested_areas[1, 2]
zeta_2_mc <- 472 / reforested_areas[2, 2]
zeta_3_mc <- 52 / reforested_areas[3, 2]

print(zeta_1_mc)
print(zeta_2_mc)
print(zeta_3_mc)

# Compute cost per hectare for each site
calib_df <- calib_df %>% mutate(
  cost_per_ha = (zeta_1_mc / 2) * share_pq_1_pasture^2 +
    (zeta_2_mc / 2) * share_pq_2_pasture^2 +
    (zeta_3_mc / 2) * share_pq_3_pasture^2
)

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



# Plot adjustment cost at average hectare
fig_1 <- calib_df %>%
  st_transform(4326) %>%
  ggplot() +
  geom_sf(aes(fill = cost_per_ha), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    fill = "USD",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/adj_costs_calib/adj_cost_at_avg.pdf", plot = fig_1)


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
