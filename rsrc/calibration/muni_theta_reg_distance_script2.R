library(stargazer)
library(ggplot2)
library(terra)
library(tidyverse)
library(sf)
library(conflicted)
library(units)
conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)


# Load municipal data
load("data/processed/muni_data.Rdata")

load("data/processed/cattle_price_index.Rdata")

# Load calibration data set
load("data/calibration/calibration_78_sites.Rdata")

# load basin geo data
map_basin <- st_read(
  dsn = "data/raw/mapbiomas/basin",
  layer = "BASIN_LEVEL_1_PNRH"
)


calib_df <- calib_df %>%
  st_transform(st_crs(map_basin))

muni_data <- muni_data %>%
  st_transform(st_crs(map_basin))

centroids <- muni_data %>%
  st_centroid() %>%
  st_coordinates()

muni_data$lon <- centroids[, "X"]
muni_data$lat <- centroids[, "Y"]




map_basin_tf <- map_basin %>%
  select(FEATURE_ID)

muni_with_basin <- st_join(muni_data, map_basin_tf, join = st_nearest_feature)



# load basin geo data
map_basin_2 <- st_read(
  dsn = "data/raw/mapbiomas/basin",
  layer = "BASIN_LEVEL_2_PNRH"
)

map_basin_tf_2 <- map_basin_2 %>%
  select(FEATURE_ID)


muni_with_basin_full <- st_join(muni_data, map_basin_tf_2, join = st_nearest_feature)

muni_to_reassign <- muni_with_basin_full %>%
  filter(FEATURE_ID %in% c(11, 100, 112,154,256)) %>%
  select(-FEATURE_ID)


map_basin_tf_reassign <- map_basin_2 %>%
  select(FEATURE_ID) %>%
  filter(!FEATURE_ID %in% c(11, 100,107, 112,154,256))

muni_to_reassign <- st_join(muni_to_reassign, map_basin_tf_reassign, join = st_nearest_feature)

muni_with_basin_2 <- bind_rows(
  muni_with_basin_full %>% filter(!FEATURE_ID %in% c(11, 100, 112,154,256)),
  muni_to_reassign
)



muni_data <- muni_with_basin


# Remove outliers from municipal data
muni_data <- muni_data[-c(142, 106, 112), ]

# Filter out observations missing distance to capital
muni_data <- muni_data %>%
  filter(!is.na(distance))




muni_data_2 <- muni_with_basin_2


# Remove outliers from municipal data
muni_data_2 <- muni_data_2[-c(142, 106, 112), ]

# Filter out observations missing distance to capital
muni_data_2 <- muni_data_2 %>%
  filter(!is.na(distance))


# Theta regression
theta_reg <- muni_data %>%
  filter(slaughter_value_per_ha_2017 > 0) %>%
  lm(
    formula = log(slaughter_value_per_ha_2017) ~
      lat +
      I(lat^2) +
      historical_temp +
      I(historical_temp^2) +
      historical_precip +
      distance +
      log(farm_gate_price_2017),
    na.action = na.exclude,
    weights = pasture_area_2017
  )
summary(theta_reg)






# Load necessary library
library(lme4)
theta_reg_re <- muni_data %>%
  filter(slaughter_value_per_ha_2017 > 0) %>%
  lmer(
    formula = log(slaughter_value_per_ha_2017) ~
      lat +
      I(lat^2) +
      historical_temp +
      I(historical_temp^2) +
      historical_precip +
      distance +
      log(farm_gate_price_2017) +
      (1 | FEATURE_ID),  # Random intercept for FEATURE_ID
    na.action = na.exclude,
    weights = pasture_area_2017
  )

# Show the summary of the random effects model
summary(theta_reg_re)



theta_reg_re_2 <- muni_data_2 %>%
  filter(slaughter_value_per_ha_2017 > 0) %>%
  lmer(
    formula = log(slaughter_value_per_ha_2017) ~
      lat +
      I(lat^2) +
      historical_temp +
      I(historical_temp^2) +
      historical_precip +
      distance +
      log(farm_gate_price_2017) +
      (1 | FEATURE_ID),  # Random intercept for FEATURE_ID
    na.action = na.exclude,
    weights = pasture_area_2017
  )




# Load necessary libraries
library(broom)
library(stargazer)



# Use stargazer to create a combined summary and output it to .txt and .tex
stargazer(theta_reg,theta_reg_re_2,
          type = "text", 
          title = "Theta Regression and Random Effects Regression Results", 
          out = "theta_reg_results_combined.txt")

stargazer(theta_reg,theta_reg_re_2,
          type = "latex", 
          title = "Theta Regression and Random Effects Regression Results", 
          out = "theta_reg_results_combined.tex")









