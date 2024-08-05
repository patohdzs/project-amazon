library(sf)
library(terra)
library(stargazer)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

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

# Get pixel areas
pixel_areas <- cellSize(pq_rst, unit = "ha")
pasture_quals <- zonal(pixel_areas, pq_rst, sum, na.rm = TRUE) %>%
  slice(-1) %>%
  mutate(prop = area / sum(area, na.rm = TRUE))

print(pasture_quals)

# Get pixel areas
pixel_areas <- cellSize(sec_veg_age_rst, unit = "ha")
area_by_age <- zonal(pixel_areas, sec_veg_age_rst, sum, na.rm = TRUE)
print(area_by_age)

# Resample pasture quality to match sec_veg resolution
pq_rst <- resample(pq_rst, sec_veg_age_rst, method = "near")

# Set areas with no two-year-old sec veg to 0
pq_rst[(sec_veg_age_rst != 2)] <- 0

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

# Multiply zetas by proportions squared
cost_one_hectare_mc <-
  (zeta_1_mc / 2) * pasture_quals[1, 3]^2 +
  (zeta_2_mc / 2) * pasture_quals[2, 3]^2 +
  (zeta_3_mc / 2) * pasture_quals[3, 3]^2

print("Cost one hectare (new, mc)")
print(cost_one_hectare_mc)

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
