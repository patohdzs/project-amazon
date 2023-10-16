
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (25 Sites MODEL)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -


library(MASS)
# Install and load dplyr package
if (!"dplyr" %in% installed.packages()) {
  install.packages("dplyr")
}
library(dplyr)

# Install and load boot package
if (!"boot" %in% installed.packages()) {
  install.packages("boot")
}
library(boot)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::summarize)
# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "calibration_25SitesModel.R script", log = T)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = paste(getwd(), "data", "_temp", sep = "/"))




# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RASTER DATA (AMAZON BIOME SHARE, PIXEL AREA, AND MAPBIOMAS CATEGORIES)
raster.25Sites <- terra::rast(list.files(paste(getwd(), "data/calibration/1055SitesModel/aux_tifs", sep = "/"),
                                         pattern = "raster_",
                                         full.names = T))


# MUNI LEVEL SPATIAL SAMPLE
load(paste(getwd(), "data/calibration/prepData/sampleMuniSpatial_prepData.Rdata", sep = "/"))






# PARAMETER THETA ------------------------------------------------------------------------------------------------------------------------------------

distance_data <- read_excel("data/calibration/ipeadata[21-08-2023-01-28].xls")
distance_data$muni_code <- as.numeric(distance_data$muni_code)

# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/calibration/prepData/muniTheta_prepData.Rdata")

# load cattle price series
load("data/calibration/prepData/seriesPriceCattle_prepData.Rdata")


# DATA MANIPULATION

# EXTRACT AVERAGE 2017 PRICE (use real prices because it is normalized to 2017 )
aux.price.2017 <-
  seriesPriceCattle.prepData %>%
  dplyr::filter(year == 2017) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(mean_price_2017 = mean(price_real_mon_cattle)/3.192) %>% # BRL to USD (commercial exchange rate - selling - average - annual - 2017 - ipeadata))
  dplyr::pull(mean_price_2017)


# REGRESSION - CATTLE VALUE (2017)

# cattle value per ha
#reg.cattleValueperHa.2017 <-
 # muniTheta.prepData %>%
  #lm(formula = cattleSlaughter_valuePerHa_2017 ~ pasture_area_2017 + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
   #    lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = pasture_area_2017)

#muniTheta.prepData <- muniTheta.prepData[-c(142, 106, 112), ]




# Remove rows from attribute data
muniTheta.prepData_data <- as.data.frame(muniTheta.prepData)  # Convert to regular dataframe
muniTheta.prepData_data <- muniTheta.prepData_data[-c(142, 106, 112), ]

# Remove geometries
geo_backup <- st_geometry(muniTheta.prepData)
geo_backup <- geo_backup[-c(142, 106, 112)]


predicted_values <- read_excel("data/calibration/farm_gate_price.xlsx")

# Combine back into an sf object
muniTheta.prepData <- st_sf(muniTheta.prepData_data, geometry = geo_backup)



# Convert to non-spatial dataframe for the merge
muniTheta_no_geo <- as.data.frame(muniTheta.prepData)

# Perform the merge
merged_data <- left_join(muniTheta_no_geo, distance_data, by = "muni_code")

# Reattach the geometry
merged_data_sf <- st_sf(merged_data, geometry = geo_backup)

muniTheta.prepData<-merged_data_sf



merged_data <- muniTheta.prepData %>%
  left_join(predicted_values, by = "muni_code") %>%
  mutate(cattleSlaughter_farmGatePrice_2017 = ifelse(is.na(cattleSlaughter_farmGatePrice_2017),
                                                     average_weighted_price,
                                                     cattleSlaughter_farmGatePrice_2017))
# select(-predicted_value_column_name)  # Remove the additional column from the result


muniTheta.prepData<-merged_data


muniTheta.prepData<- muniTheta.prepData %>%
  filter(!is.na(distance))



muniTheta.prepData_filtered <- muniTheta.prepData %>%
  filter(cattleSlaughter_valuePerHa_2017 > 0) # Exclude zeros
#muniTheta.prepData_filtered <- na.omit(muniTheta.prepData_filtered)



reg.cattleValueperHa.2017 <-
  muniTheta.prepData_filtered  %>%
  lm(formula = log(cattleSlaughter_valuePerHa_2017) ~  historical_precip  + historical_temp + I(historical_temp^2)+
       lat + I(lat^2)+cattleSlaughter_farmGatePrice_2017+distance, na.action = na.exclude, weights = pasture_area_2017)

summary(reg.cattleValueperHa.2017)




prior_coe <- coef(reg.cattleValueperHa.2017)
regressand_df <- data.frame(log_cattleSlaughter_valuePerHa_2017 = log(muniTheta.prepData_filtered$cattleSlaughter_valuePerHa_2017))
regressor_df <- as.data.frame(reg.cattleValueperHa.2017$model[-1])
cols_to_scale <- setdiff(names(regressor_df), "(weights)")
scaled_data <- regressor_df
scaled_data <-scaled_data %>%
  mutate(historical_precip = (historical_precip-mean(muniTheta.prepData_filtered$historical_precip))/sd(muniTheta.prepData_filtered$historical_precip)) %>%
  mutate(historical_temp = (historical_temp-mean(muniTheta.prepData_filtered$historical_temp))/sd(muniTheta.prepData_filtered$historical_temp)) %>%
  mutate(`I(historical_temp^2)` = (`I(historical_temp^2)`-mean(muniTheta.prepData_filtered$historical_temp^2))/sd(muniTheta.prepData_filtered$historical_temp^2)) %>%
  mutate(lat = (lat-mean(muniTheta.prepData_filtered$lat))/sd(muniTheta.prepData_filtered$lat)) %>%
  mutate(`I(lat^2)` =(`I(lat^2)`-mean(muniTheta.prepData_filtered$lat^2))/sd(muniTheta.prepData_filtered$lat^2)) %>%
  mutate(cattleSlaughter_farmGatePrice_2017=(cattleSlaughter_farmGatePrice_2017-mean(muniTheta.prepData_filtered$cattleSlaughter_farmGatePrice_2017))/sd(muniTheta.prepData_filtered$cattleSlaughter_farmGatePrice_2017))%>%
  mutate(distance = (distance-mean(muniTheta.prepData_filtered$distance))/sd(muniTheta.prepData_filtered$distance))


weights <- muniTheta.prepData_filtered $pasture_area_2017

# Create the weight matrix
W <- diag(weights)
W_half <- diag(sqrt(weights))


# Combine with intercept and Y
X_df <- cbind(1, scaled_data)
X_df <- X_df %>%
  dplyr::rename(weights_var = `(weights)`)
X_df <- X_df %>%
  dplyr::select(-`weights_var`)

transformed_X_matrix <- W_half %*% as.matrix(X_df)
transformed_X_df <- as.data.frame(transformed_X_matrix)

# Set the column names again from the original X_df
colnames(transformed_X_df) <- colnames(X_df)
colnames(transformed_X_df)[1] <- "cons"

Y_matrix <- as.matrix(regressand_df)

transformed_Y_matrix <- W_half %*% Y_matrix

# Convert back to dataframe
transformed_Y_df <- data.frame(log_cattleSlaughter_valuePerHa_2017 = transformed_Y_matrix)
Y <- as.matrix(transformed_Y_df)

Var_df <- cbind(transformed_X_df, Y)

reg.prior <-
  Var_df  %>%
  lm(formula = log_cattleSlaughter_valuePerHa_2017 ~ -1+cons + historical_precip  + historical_temp + `I(historical_temp^2)`+lat
       + `I(lat^2)`+cattleSlaughter_farmGatePrice_2017+distance, na.action = na.exclude)

summary(reg.prior)



fit <- lm.fit(transformed_X_matrix, transformed_Y_matrix)
print(coef(fit))

residuals <- transformed_Y_matrix - transformed_X_matrix %*% fit$coefficients
sigma2 <- sum(residuals^2) / fit$df.residual

# Calculate variance-covariance matrix of coefficients
var_beta <- sigma2 * solve(t(transformed_X_matrix) %*% transformed_X_matrix)


#summary(reg.prior)
prior_coe <- coefficients(reg.prior)


X <- as.matrix(transformed_X_df)


# Number of observations
N <- nrow(X)
Lambda_t0 <- matrix(0, nrow=8, ncol=8)
b_t0 <- as.matrix(prior_coe)
c_t0<- 555
#zeta_ini <- rgamma(1, shape = c_t0, rate = d_t0)
#zeta=1/variances
zeta_ini <- 1/sigma2
d_t0 <-c_t0/zeta_ini

Lambda_t1 <- Lambda_t0+t(X)%*%  X
b_t1 <- tryCatch({
  # Attempt to solve the system
  solve(Lambda_t1,t(X)%*% Y)
}, error = function(e) {
  # If there's an error (because Lambda_t1 is singular), return b_t0
  cat("Matrix is singular, setting b_t1 = b_t0\n")
  b_t0
})
c_t1<- c_t0+1
d_t1<- d_t0+t(Y) %*% Y - t(b_t1) %*% Lambda_t1 %*% b_t1
d_t1_value <- d_t1[1, 1]


p <- length(b_t1)  # Assuming b_t1 is the same size as beta_sample
beta_vec <- matrix(nrow = 0, ncol = p)
zeta_vec <- NULL
index_vec<-NULL



for(j in 1: 100000 )
{

  zeta_sample <- rgamma(1, shape = c_t1/2, rate = d_t1_value/2)
  cov_matrix <- solve(zeta_sample * Lambda_t1)
  beta_sample <- mvrnorm(1, mu = as.vector(b_t1), Sigma = cov_matrix)
  d_t1new<- d_t0+t(Y) %*% Y - t(b_t1) %*% Lambda_t1 %*% b_t1
  if (d_t1new>0){
    d_t1<- d_t1new
  } else {
    print("d value is negative")
  }
  d_t1_value <- d_t1[1, 1]


  if(j>50000){
    beta_vec <- rbind(beta_vec, beta_sample)
    zeta_vec<-c(zeta_vec,zeta_sample)
    index_vec<-c(index_vec,j)
  }
}



beta_sample_df <- as.data.frame(beta_vec)
readr::write_csv(beta_sample_df, file = paste(getwd(), "data/hmc/", "theta_coe.csv", sep = "/"))




# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
