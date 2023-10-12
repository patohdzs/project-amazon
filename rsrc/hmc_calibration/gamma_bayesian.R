
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

library(MASS)
conflicts_prefer(dplyr::select)
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


# INITIAL CONDITIONS Z -------------------------------------------------------------------------------------------------------------------------------

# AGGREGATE FROM 1000 Sites TO 25 Sites
# transform shares to areas
raster.25Sites$amazonBiomeArea_ha_25Sites <- raster.25Sites$share_amazonBiome*raster.25Sites$pixelArea_ha
raster.25Sites$forestArea_1995_ha_25Sites <- raster.25Sites$share_forest_1995*raster.25Sites$pixelArea_ha
raster.25Sites$agriculturalUseArea_1995_ha_25Sites <- raster.25Sites$share_agriculturalUse_1995*raster.25Sites$pixelArea_ha
raster.25Sites$otherArea_1995_ha_25Sites <- raster.25Sites$share_other_1995*raster.25Sites$pixelArea_ha
raster.25Sites$forestArea_2017_ha_25Sites <- raster.25Sites$share_forest_2017*raster.25Sites$pixelArea_ha
raster.25Sites$agriculturalUseArea_2017_ha_25Sites <- raster.25Sites$share_agriculturalUse_2017*raster.25Sites$pixelArea_ha
raster.25Sites$otherArea_2017_ha_25Sites <- raster.25Sites$share_other_2017*raster.25Sites$pixelArea_ha

# select area variables
raster.25Sites <- terra::subset(raster.25Sites,
                                c("amazonBiomeArea_ha_25Sites", "pixelArea_ha",
                                  "forestArea_1995_ha_25Sites", "agriculturalUseArea_1995_ha_25Sites", "otherArea_1995_ha_25Sites",
                                  "forestArea_2017_ha_25Sites", "agriculturalUseArea_2017_ha_25Sites", "otherArea_2017_ha_25Sites"))

# aggregate from 1000 Sites to 25
raster.25Sites <- terra::aggregate(raster.25Sites, fact = 8, fun = sum, na.rm = T)

# extract variables as polygons, transform to sf, and project data for faster spatial manipulation
calibration.25SitesModel <- terra::as.polygons(raster.25Sites, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)

# transform share aggregate in area (ha)
calibration.25SitesModel <-
  calibration.25SitesModel %>%
  dplyr::mutate(zbar_1995_25Sites = agriculturalUseArea_1995_ha_25Sites + forestArea_1995_ha_25Sites,
                zbar_2017_25Sites = agriculturalUseArea_2017_ha_25Sites + forestArea_2017_ha_25Sites) %>%
  dplyr::select(amazonBiomeArea_ha_25Sites, siteArea_ha_25Sites = pixelArea_ha,
                forestArea_1995_ha_25Sites,
                z_1995_25Sites = agriculturalUseArea_1995_ha_25Sites, zbar_1995_25Sites,
                forestArea_2017_ha_25Sites,
                z_2017_25Sites = agriculturalUseArea_2017_ha_25Sites, zbar_2017_25Sites)

# remove Sites with less than 1% of its are intersecting with the amazon biome
calibration.25SitesModel <-
  calibration.25SitesModel %>%
  dplyr::filter(amazonBiomeArea_ha_25Sites/siteArea_ha_25Sites >= 0.01)

# add id variable
calibration.25SitesModel$id <- 1:nrow(calibration.25SitesModel)



# PARAMETER GAMMA ------------------------------------------------------------------------------------------------------------------------------------


# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/calibration/muniTheta_prepData_gamma.Rdata")




muniTheta.prepData<-muniTheta.prepData %>%
  dplyr::mutate(co2e_ha_2017 = (agb_2017/2)*(44/12))



muniTheta.prepData.filter<- muniTheta.prepData %>%
  filter(!is.na(co2e_ha_2017))

#reg.share.2017 <-
#  muniTheta.prepData  %>%
#  lm(formula = log(share)  ~ log(lat)+log(lon), na.action = na.exclude)

#summary(reg.share.2017)


#residuals_values <- residuals(reg.share.2017)
#residuals_values_df<-data.frame(residuals_values)


#muniTheta.prepData<-muniTheta.prepData%>%
#  dplyr::mutate(residuals=residuals_values_df$residuals_values)


# DATA INPUT (2017)
# load pixel sample with biomass data
# Load pixelBiomass2017_prepData.Rdata




reg.gamma.2017 <-
  muniTheta.prepData.filter  %>%
  lm(formula = log(co2e_ha_2017)  ~ log(historical_precip) + log(historical_temp) +log(lat)+log(lon), na.action = na.exclude)

summary(reg.gamma.2017)


regressand_df <- data.frame(log_co2e_ha_2017 = log(muniTheta.prepData.filter$co2e_ha_2017))
regressor_df <- as.data.frame(reg.gamma.2017$model[-1])

cols_to_scale <- setdiff(names(regressor_df), "(weights)")
scaled_data <- regressor_df
scaled_data[, cols_to_scale] <- scale(regressor_df[, cols_to_scale], center = TRUE, scale=TRUE)





X_df <- cbind(1, scaled_data)
X <- as.matrix(X_df)
Y <-as.matrix(regressand_df)
Var_df <- cbind(X_df, Y)


reg.prior <-
  Var_df  %>%
  lm(formula = log_co2e_ha_2017  ~ `log(historical_precip)` + `log(historical_temp)` +`log(lat)`+`log(lon)`, na.action = na.exclude)

summary(reg.prior)
reg_summary <- summary(reg.prior)






prior_coe <- coef(reg.prior)

#vcov_matrix <- vcov(reg.gamma.2017)
sigma2<-reg_summary$sigma^2

# Number of observations
N <- nrow(X)
Lambda_t0 <- matrix(0, nrow=5, ncol=5)
b_t0 <- as.matrix(prior_coe)
c_t0<- 533
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

  zeta_sample <- rgamma(1, shape = c_t1, rate = d_t1_value)
  cov_matrix <- solve(zeta_sample * Lambda_t1)
  beta_sample <- mvrnorm(1, mu = as.vector(b_t1), Sigma = cov_matrix)
  d_t1new<- d_t0+t(Y) %*% Y - t(beta_sample) %*% Lambda_t1 %*% beta_sample
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
readr::write_csv(beta_sample_df, file = paste(getwd(), "data/HMC_norm/", "gamma_coe.csv", sep = "/"))
