# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (1043 SITES MODEL)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(sf)
library(units)
library(terra)
library(readxl)
library(tidyverse)
library(tictoc)

# Start timer
tic(msg = "calibrate_1043_sites_model.R script", log = TRUE)

# Load spatial municipal sample
load("data/processed/spatial_muni_sample.Rdata")

# Load municipal data
load("data/processed/muni_data.Rdata")

# Load state emissions data
load("data/processed/state_emissions.Rdata")

# Load cattle price index
load("data/processed/cattle_price_index.Rdata")

# Load distance to capital data
dist_to_capital <-
  read_excel("data/raw/ipea/distance_to_capital/ipeadata[21-08-2023-01-28].xls")

# Load raster data (amazon biome share, pixel areas, and land uses)
raster_variables <- rast(
  list.files(
    "data/processed/",
    pattern = "amazon_",
    full.names = TRUE
  )
)

# Extract variables as polygons and project data
calib_df <- as.polygons(raster_variables, dissolve = FALSE) %>%
  st_as_sf() %>%
  st_transform(5880)

# Remove sites with less than 3% overlap with the amazon biome
calib_df <- calib_df %>%
  filter(share_amazon_biome >= 0.03)

# Add id variable
calib_df$id <- seq_len(nrow(calib_df))

# Transform share variables into area (ha)
calib_df <- calib_df %>%
  rename(site_area_ha = pixel_area_ha) %>%
  mutate(
    amazon_biome_area_ha = share_amazon_biome * site_area_ha,
    forest_area_1995_ha = share_forest_1995 * site_area_ha,
    forest_area_2008_ha = share_forest_2008 * site_area_ha,
    forest_area_2017_ha = share_forest_2017 * site_area_ha,
    other_area_1995_ha = share_other_1995 * site_area_ha,
    other_area_2008_ha = share_other_2008 * site_area_ha,
    other_area_2017_ha = share_other_2017 * site_area_ha,
  ) %>%
  mutate(
    z_1995 = share_agricultural_use_1995 * site_area_ha,
    z_2008 = share_agricultural_use_2008 * site_area_ha,
    z_2017 = share_agricultural_use_2017 * site_area_ha,
    zbar_1995 = forest_area_1995_ha + z_1995,
    zbar_2008 = forest_area_2008_ha + z_2008,
    zbar_2017 = forest_area_2017_ha + z_2017,
  ) %>%
  select(id, contains("area"), contains("z"))

# Convert biomass to CO2 equivalent
muni_data <- muni_data %>%
  mutate(co2e_ha_2017 = (agb_2017 / 2) * (44 / 12))

# Regress log-gamma on geographic covariates
reg_gamma_2017 <- lm(
  formula = log(co2e_ha_2017) ~
    log(historical_precip) +
    log(historical_temp) +
    log(lat) +
    log(lon),
  data = muni_data,
  na.action = na.exclude
)

# Compute fitted values and transform back into levels
muni_data <- muni_data %>%
  mutate(co2e_ha_2017_fitted = exp(predict(reg_gamma_2017, .)))

# Match minicells with sites
site_gamma_2017 <- calib_df %>%
  st_intersection(muni_data) %>%
  select(
    id,
    muni_code,
    muni_area,
    co2e_ha_2017,
    co2e_ha_2017_fitted,
    historical_precip,
    historical_temp,
    lat,
    lon
  )

# Set area units
site_gamma_2017$muni_site_area <-
  st_area(site_gamma_2017) %>%
  set_units(ha) %>%
  unclass()

# Drop spatial feature
site_gamma_2017 <- site_gamma_2017 %>%
  st_drop_geometry()

# Average carbon density on primary forest areas by site
site_gamma_2017 <- site_gamma_2017 %>%
  group_by(id) %>%
  summarise(
    gamma = weighted.mean(
      co2e_ha_2017_fitted,
      w = muni_site_area,
      na.rm = TRUE
    )
  )

# Add site gamma to calibration variables
calib_df <- calib_df %>% left_join(site_gamma_2017)

# Clean environment
rm(site_gamma_2017)

# Identify adjacent neighbors
aux_neighbors <- st_is_within_distance(
  calib_df,
  calib_df,
  dist = 100,
  remove_self = TRUE
)

# Impute values for missing gammas using the average of adjacent neighbors
calib_df <- calib_df %>%
  mutate(gamma = if_else(is.na(gamma),
    apply(aux_neighbors, 1, function(i) {
      mean(.$gamma[i], na.rm = TRUE)
    }),
    gamma
  ))

# Estimate of alpha same as in the global model
calib_df <- calib_df %>%
  mutate(alpha = 1 - (1 - 0.99)^(1 / 100))


# Calculate average net emission factor
#   from agricultural use across years and states
avg_net_emission_factor <- state_emissions %>%
  ungroup() %>%
  summarise(
    nef = sum(net_emissions_co2e) / sum(agriculturalUse_area)
  ) %>%
  pull(nef)

# Estimate kappa
calib_df <- calib_df %>%
  mutate(kappa = avg_net_emission_factor)

# Clean environment
rm(avg_net_emission_factor)


# Zeta is calibrated such that the marginal cost of changing land use
#   (zeta*forest_to_pasture_transition_area) matches the forest to
#   pasture transition cost estimated by Araujo, Costa and Sant'Anna (2022),
#   reported on column 4 of table 4 (on the right).

# We transform to dollars using an FX rate of 4.14 (December 2019)
aux_transition_cost <- 1614.54 / 4.14

# The forest area in 2008 represented 72% of the Legal Amazon area,
#   which covers 501,506,775 ha, so the transition in hectares is
#   0.065*0.72*501,506,775 = 23,470,517,
#   resulting in an annual average of 2,347,052 ha.
aux_transition_area <- (0.065 * 0.72 * 501506775) / (2017 - 2008 + 1)
zeta <- aux_transition_cost / aux_transition_area

# Alternative value based on (https://www.otempo.com.br/brasil/investigacoes-revelam-quadrilhas-e-ganho-milionario-por-tras-do-desmate-1.2229571)
zeta_alt <- 483 / aux_transition_area

# Estimate of zeta same as in the global model
calib_df <- calib_df %>%
  mutate(zeta = zeta, zeta_alt = zeta_alt)

# x_2017 is estimated as the carbon stock stored in forest areas,
#   assuming that all forests are primary.
calib_df <- calib_df %>%
  mutate(
    x_1995 = gamma * (zbar_1995 - z_1995),
    x_2008 = gamma * (zbar_2008 - z_2008),
    x_2017 = gamma * (zbar_2017 - z_2017)
  )



# Extract average cattle price index
# BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
brl_to_usd <- 3.192

aux_price_2017 <-
  cattle_price_index %>%
  filter(year == 2017) %>%
  group_by(year) %>%
  summarise(mean_price_2017 = mean(price_real_mon_cattle) / brl_to_usd) %>%
  pull(mean_price_2017)

# Remove geometries
geo_backup <- st_geometry(muni_data)
geo_backup <- geo_backup[-c(142, 106, 112)]

# Remove rows from attribute data
muni_data <- as.data.frame(muni_data)
muni_data <- muni_data[-c(142, 106, 112), ]

# Combine back into an sf object
muni_data <- st_sf(muni_data, geometry = geo_backup)

# Convert to non-spatial dataframe for the merge
muniTheta_no_geo <- as.data.frame(muni_data)

# Merge municipal data with distance to capital data
dist_to_capital$muni_code <- as.numeric(dist_to_capital$muni_code)
merged_data <- left_join(muniTheta_no_geo, dist_to_capital, by = "muni_code")

# Reattach the geometry
muni_data <- st_sf(merged_data, geometry = geo_backup)

# Filter out observations missing distance to capital
muni_data <- muni_data %>%
  filter(!is.na(distance))

# Exclude zeros
muni_data_filtered <- muni_data %>%
  filter(slaughter_value_per_ha_2017 > 0)

# Theta regression
reg_theta <- muni_data_filtered %>%
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

# Extract fitted values
muni_data <- muni_data %>%
  mutate(slaughter_value_per_ha_fitted = exp(predict(reg_theta, .)))

# Extract minimum positive fitted value
aux_positive_theta_fitted <- muni_data %>%
  filter(slaughter_value_per_ha_fitted > 0) %>%
  pull(slaughter_value_per_ha_fitted) %>%
  min()

# Winsorize cattleSlaughter_valuePerHa_fitted
muni_data <- muni_data %>%
  mutate(d_theta_winsorized_2017 = if_else(slaughter_value_per_ha_fitted <= 0,
    1,
    0
  ))
# clean environment
rm(reg_theta, aux_positive_theta_fitted)

# Match munis with sites
site_theta_2017 <- st_intersection(
  calib_df %>% select(id),
  muni_data %>% select(
    muni_code, muni_area, slaughter_value_per_ha_fitted,
    pasture_area_2017, d_theta_winsorized_2017
  )
)


# Calculate muni areas inside each site
site_theta_2017$muni_site_area <-
  st_area(site_theta_2017) %>%
  set_units(ha) %>%
  unclass()

# Drop spatial feature
site_theta_2017 <-
  site_theta_2017 %>%
  st_drop_geometry()

# Calculate theta and pasture_area by site
# (for each muni adjust the value by the share of the muni area inside the site)
aux_theta_2017 <-
  site_theta_2017 %>%
  group_by(id) %>%
  summarise(
    theta2017_1043Sites = weighted.mean(
      slaughter_value_per_ha_fitted / aux_price_2017,
      w = muni_site_area,
      na.rm = TRUE
    ),
    pasture_area_2017 = sum(pasture_area_2017 * (muni_site_area / muni_area), na.rm = TRUE),
    d_theta_winsorized_2017 = min(d_theta_winsorized_2017, na.rm = T)
  )

# Add cattle_slaughter_value_fitted and pasture_area to spatial variables
calib_df <- left_join(calib_df, aux_theta_2017)

# Clean environment
rm(aux_theta_2017)

# Filter out missing vlaues
calib_df <- calib_df %>%
  filter(!is.na(theta2017_1043Sites))

# Calculate average theta using the values of 2006 and 2017
calib_df <-
  calib_df %>%
  group_by(id) %>%
  mutate(theta_1043Sites = rowMeans(across(starts_with("theta20")), na.rm = T)) %>%
  ungroup()

# Change from BRL to USD (commercial exchange rate - selling - average - monthly - january/2017 - ipeadata)
cattle_price_index <- cattle_price_index %>%
  mutate(price_real_mon_cattle = price_real_mon_cattle / 3.1966) %>%
  filter(year >= 1995, year <= 2017) # select time period

# Add site ID's
calib_df <- calib_df %>%
  mutate(id = seq_len(nrow(calib_df)))

# Save site boundaries
calib_df %>%
  select(id) %>%
  st_write(
    "data/calibration/grid_1043_sites.geojson",
    driver = "GeoJSON",
    delete_dsn = TRUE
  )

# Save calibration data
save(calib_df, file = "data/calibration/calibration_1043_sites.Rdata")

# Remove spatial feature
calib_df %>%
  st_drop_geometry() %>%
  write_csv(file = "data/calibration/calibration_1043_sites.csv")

# End timer
toc(log = TRUE)
