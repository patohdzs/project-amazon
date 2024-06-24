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
library(terra)
library(readxl)
library(tidyverse)
library(tictoc)

# Start timer
tic(msg = "calibration_1043_sites.R script", log = TRUE)

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
calibration_1043_sites <- as.polygons(raster_variables, dissolve = FALSE) %>%
  st_as_sf() %>%
  st_transform(5880)

# Remove sites with less than 1% overlap with the amazon biome
calibration_1043_sites <-
  calibration_1043_sites %>%
  filter(share_amazonBiome >= 0.03)

# Add id variable
calibration_1043_sites$id <- seq_len(nrow(calibration_1043_sites))

# Transform share variables in area (ha)
calibration_1043_sites <-
  calibration_1043_sites %>%
  mutate(
    amazon_biome_area_ha = share_amazonBiome * pixelArea_ha,
    forest_area_1995_ha = share_forest_1995 * pixelArea_ha,
    forest_area_2008_ha = share_forest_2008 * pixelArea_ha,
    forest_area_2017_ha = share_forest_2017 * pixelArea_ha,
    other_area_1995_ha = share_other_1995 * pixelArea_ha,
    other_area_2008_ha = share_other_2008 * pixelArea_ha,
    other_area_2017_ha = share_other_2017 * pixelArea_ha,
  ) %>%
  mutate(
    z_1995 = share_agriculturalUse_1995 * pixelArea_ha,
    z_2008 = share_agriculturalUse_2008 * pixelArea_ha,
    z_2017 = share_agriculturalUse_2017 * pixelArea_ha,
    zbar_1995 = forest_area_1995_ha + z_1995,
    zbar_2008 = forest_area_2008_ha + z_2008,
    zbar_2017 = forest_area_2017_ha + z_2017,
  ) %>%
  select(id,
    siteArea_ha_1043Sites = pixelArea_ha, amazon_biome_area_ha,
    forest_area_2017_ha, other_area_2017_ha, z_2017, zbar_2017,
    forest_area_1995_ha, other_area_1995_ha, z_1995, zbar_1995,
    forest_area_2008_ha, other_area_2008_ha, z_2008, zbar_2008
  )


# Convert biomass to CO2 equivalent
muni_data <- muni_data %>%
  mutate(co2e_ha_2017 = (agb_2017 / 2) * (44 / 12))

# Regress log-gamma on geographic covariates
reg_gamma_2017 <- lm(
  formula =
    log(co2e_ha_2017) ~
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
site_gamma_2017 <-
  calibration_1043_sites %>%
  st_intersection(muni_data) %>%
  select(id, muni_code, muni_area, co2e_ha_2017, co2e_ha_2017_fitted, historical_precip, historical_temp, lat, lon)

# Set area units
site_gamma_2017$muni_site_area <-
  st_area(site_gamma_2017) %>%
  set_units(ha) %>%
  unclass()

# Drop spatial feature
site_gamma_2017 <-
  site_gamma_2017 %>%
  st_drop_geometry()

# Average carbon density on primary forest areas by site
site_gamma_2017 <-
  site_gamma_2017 %>%
  group_by(id) %>%
  summarise(
    gamma = weighted.mean(
      co2e_ha_2017_fitted,
      w = muni_site_area,
      na.rm = TRUE
    )
  )


# Add site gamma to calibration variables
calibration_1043_sites <- calibration_1043_sites %>% left_join(site_gamma_2017)

# Clean environment
rm(site_gamma_2017)

# Identify adjacent neighbors
aux_neighbors <- st_is_within_distance(
  calibration_1043_sites,
  calibration_1043_sites,
  dist = 100,
  remove_self = TRUE
)

# Impute values for missing gammas using the average of adjacent neighbors
calibration_1043_sites <-
  calibration_1043_sites %>%
  mutate(gamma = if_else(is.na(gamma),
    apply(aux_neighbors, 1, function(i) {
      mean(.$gamma[i], na.rm = TRUE)
    }),
    gamma
  ))

# Estimate of alpha same as in the global model
calibration_1043_sites <-
  calibration_1043_sites %>%
  mutate(alpha = 1 - (1 - 0.99)^(1 / 100))


# Calculate average net emission factor from agricultural use across years and states
avg_net_emission_factor <-
  state_emissions %>%
  ungroup() %>%
  summarise(
    net_emission_factor_co2e_ha = sum(net_emissions_co2e) / sum(agricultural_use_area)
  ) %>%
  pull(net_emission_factor_co2e_ha)

# Estimate kappa
calibration_1043_sites <-
  calibration_1043_sites %>%
  mutate(kappa = avg_net_emission_factor)

# Clean environment
rm(avg_net_emission_factor)


# Zeta is calibrated such that the marginal cost of changing land use
#   (zeta*forestToPastureTransitionArea) matches the forest to
#   pasture transition cost estimated by Araujo, Costa and Sant'Anna (2022),
#   reported on column 4 of table 4 (on the right).

# We transform to dollars using an FX rate of 4.14 (December 2019) as in the paper:
#   1614.54/4.14 = 390 USD/ha. The paper also calculates that 6.5% of the forest area
#   was converted to pasture between 2008 and 2017.
aux_transition_cost <- 1614.54 / 4.14

# The forest area in 2008 represented 72% of the Legal Amazon area,
#   which covers 501,506,775 ha, so the transition in hectares is
#   0.065*0.72*501,506,775 = 23,470,517, resulting in an annual average of 2,347,052 ha.
aux_transition_area <- (0.065 * 0.72 * 501506775) / (2017 - 2008 + 1)
zeta <- aux_transition_cost / aux_transition_area

# Alternative value based on (https://www.otempo.com.br/brasil/investigacoes-revelam-quadrilhas-e-ganho-milionario-por-tras-do-desmate-1.2229571)
zeta_alt <- 483 / aux_transition_area

# Estimate of zeta same as in the global model
calibration_1043_sites <-
  calibration_1043_sites %>%
  mutate(zeta = zeta, zeta_alt = zeta_alt)

# x_2017 is estimated as the carbon stock stored in forest areas,
#   assuming that all forests are primary.
calibration_1043_sites <-
  calibration_1043_sites %>%
  mutate(
    x_1995 = gamma * (zbar_1995 - z_1995),
    x_2008 = gamma * (zbar_2008 - z_2008),
    x_2017 = gamma * (zbar_2017 - z_2017)
  )


dist_to_capital$muni_code <- as.numeric(dist_to_capital$muni_code)


# Extract average cattle price index
# BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
brl_to_usd <- 3.192

aux_price_2017 <-
  cattle_price_index %>%
  filter(year == 2017) %>%
  group_by(year) %>%
  summarise(mean_price_2017 = mean(price_real_mon_cattle) / brl_to_usd) %>%
  pull(mean_price_2017)

# Remove rows from attribute data
muni_data_df <- as.data.frame(muni_data)
muni_data_df <- muni_data_df[-c(142, 106, 112), ]

# Remove geometries
geo_backup <- st_geometry(muni_data)
geo_backup <- geo_backup[-c(142, 106, 112)]

# Combine back into an sf object
muni_data <- st_sf(muni_data_df, geometry = geo_backup)

# Convert to non-spatial dataframe for the merge
muniTheta_no_geo <- as.data.frame(muni_data)

# Perform the merge
merged_data <- left_join(muniTheta_no_geo, dist_to_capital, by = "muni_code")

# Reattach the geometry
muni_data <- st_sf(merged_data, geometry = geo_backup)

muni_data <- muni_data %>%
  filter(!is.na(distance))

# Exclude zeros
muni_data_filtered <- muni_data %>%
  filter(cattleSlaughter_valuePerHa_2017 > 0)

# Theta regression
reg_theta <-
  muni_data_filtered %>%
  lm(formula = log(cattleSlaughter_valuePerHa_2017) ~ historical_precip + historical_temp + I(historical_temp^2)
    + lat + I(lat^2) + distance + log(cattleSlaughter_farmGatePrice_2017), na.action = na.exclude, weights = pasture_area_2017)
# regression results
summary(reg_theta)

# extract fitted values
muni_data <-
  muni_data %>%
  mutate(cattleSlaughter_valuePerHa_fitted_2017 = exp(predict(reg_theta, .)))

# extract minimum positive fitted value
aux_min_positive_cattleSlaughter_value_ha_fitted_2017 <-
  muni_data %>%
  filter(cattleSlaughter_valuePerHa_fitted_2017 > 0) %>%
  pull(cattleSlaughter_valuePerHa_fitted_2017) %>%
  min()

# winsorize cattleSlaughter_valuePerHa_fitted
muni_data <-
  muni_data %>%
  mutate(d_theta_winsorized_2017 = if_else(cattleSlaughter_valuePerHa_fitted_2017 <= 0,
    1,
    0
  ))
# clean environment
rm(reg_theta, aux_min_positive_cattleSlaughter_value_ha_fitted_2017)

# match munis with sites
site_theta_2017 <- st_intersection(
  calibration_1043_sites %>% select(id),
  muni_data %>% select(
    muni_code, muni_area, cattleSlaughter_valuePerHa_fitted_2017,
    pasture_area_2017, d_theta_winsorized_2017
  )
)



# calculate muni areas inside each site
site_theta_2017$muni_site_area <-
  st_area(site_theta_2017) %>%
  set_units(ha) %>%
  unclass()

# drop spatial feature
site_theta_2017 <-
  site_theta_2017 %>%
  st_drop_geometry()

# calculate cattleSlaughter_valuePerHa_fitted and pastureArea_value by site (for each muni adjust the value by the share of the muni area inside the site)
aux_theta_2017 <-
  site_theta_2017 %>%
  group_by(id) %>%
  summarise(
    theta2017_1043Sites = weighted.mean(cattleSlaughter_valuePerHa_fitted_2017 / aux_price_2017, w = muni_site_area, na.rm = T),
    pasture_area_2017 = sum(pasture_area_2017 * (muni_site_area / muni_area), na.rm = T),
    d_theta_winsorized_2017 = min(d_theta_winsorized_2017, na.rm = T)
  )

# Add cattle_slaughter_valuePerHa_fitted and pastureArea_value to spatial variables
calibration_1043_sites <- left_join(calibration_1043_sites, aux_theta_2017)

# clean environment
rm(aux_theta_2017)


calibration_1043_sites <-
  calibration_1043_sites %>%
  filter(!is.na(theta2017_1043Sites))


# calculate average and SD theta using the values of 2006 and 2017
calibration_1043_sites <-
  calibration_1043_sites %>%
  group_by(id) %>%
  mutate(theta_1043Sites = rowMeans(across(starts_with("theta20")), na.rm = T)) %>%
  ungroup()



# PRICE INITIAL CONDITION AND TRANSITIONS ------------------------------------------------------------------------------------------------------------

cattle_price_index <-
  cattle_price_index %>%
  mutate(price_real_mon_cattle = price_real_mon_cattle / 3.1966) %>% # change from BRL to USD (commercial exchange rate - selling - average - monthly - january/2017 - ipeadata)
  filter(year >= 1995, year <= 2017) # select time period


# 2 PRICES (LOW X HIGH)
cattle_price_index <-
  cattle_price_index %>%
  mutate(
    price_low = quantile(price_real_mon_cattle, 0.25), # define low price value as the 33th percentile
    price_high = quantile(price_real_mon_cattle, 0.75), # define high price value as the 66th percentile
    price_median = quantile(price_real_mon_cattle, 0.5),
    price_mean = mean(price_real_mon_cattle)
  )


# create discretized version of the price series (high and low values only - using 25th and 75th percentiles as cut-offs)
cattle_price_index$d_high <- as.numeric(NA) # initialize dummy indicating if the price is high or low
cattle_price_index[1, "d_high"] <- 1 # initial value set to 1 because the first price is the highest of the series

for (i in 2:nrow(cattle_price_index)) {
  # a change from high to low only occurs if price reaches a value below the low
  if (cattle_price_index[i - 1, "d_high"] == 1 & cattle_price_index[i, "price_real_mon_cattle"] > cattle_price_index[i, "price_low"]) {
    cattle_price_index[i, "d_high"] <- 1
  } else if (cattle_price_index[i - 1, "d_high"] == 1 & cattle_price_index[i, "price_real_mon_cattle"] < cattle_price_index[i, "price_low"]) {
    cattle_price_index[i, "d_high"] <- 0
    # a change from low to high only occurs if price reaches a value above the high
  } else if (cattle_price_index[i - 1, "d_high"] == 0 & cattle_price_index[i, "price_real_mon_cattle"] > cattle_price_index[i, "price_high"]) {
    cattle_price_index[i, "d_high"] <- 1
  } else if (cattle_price_index[i - 1, "d_high"] == 0 & cattle_price_index[i, "price_real_mon_cattle"] < cattle_price_index[i, "price_high"]) {
    cattle_price_index[i, "d_high"] <- 0
  }
}

# construct discretized version of prices (2 states: high x low)
cattle_price_index <-
  cattle_price_index %>%
  mutate(discrete_2prices = if_else(d_high == 1, price_high, price_low)) %>%
  select(date, year, month, price_real_mon_cattle, discrete_2prices, price_high, price_median, price_mean, price_low)


# calculate probability transition matrix (give the same results as manually computing the number of consecutive prices at the same level divided by the ocurrence of that price level)
matrixTransition.2prices <- markovchainFit(cattle_price_index$discrete_2prices)$estimate@transitionMatrix


# STORE PARAMETER VALUES
calibration_1043_sites <-
  calibration_1043_sites %>%
  mutate(p_2017_1043Sites = max(cattle_price_index$price_high))





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# ORDER VARIABLES
calibration_1043_sites <-
  calibration_1043_sites %>%
  select(
    id, z_2017, zbar_2017, x_2017_1043Sites, gamma_1043Sites, theta_1043Sites,
    d_theta_winsorized_2017, pasture_area_2017, ends_with("_1043Sites")
  )



calibration_1043_sites <-
  calibration_1043_sites %>%
  mutate(id = 1:nrow(calibration_1043_sites))




id <- calibration_1043_sites %>%
  select(id)

st_write(id, "data/calibration/hmc/id_1043.geojson", driver = "GeoJSON", delete_dsn = TRUE)



# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(calibration_1043_sites,
  file = "data/calibration/hmc/hmc_1043SitesModel.Rdata"
)

# remove spatial feature
calibration_1043_sites <- calibration_1043_sites %>% st_drop_geometry()

write_csv(calibration_1043_sites,
  file = "data/calibration/hmc/hmc_1043SitesModel.csv"
)


# CLEAN TEMP DIR
tmpFiles(current = T, remove = T)
gc()



# END TIMER
toc(log = TRUE)




# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
