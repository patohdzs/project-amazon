
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (8 SITES MODEL)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# GROUNDHOG (REPRODUCIBILITY SOLUTION TO HANDLING DIFFERENT VERSIONS OF R AND ITS PACKAGES)

# check if groundhog is installed and load it
if ("groundhog" %in% installed.packages()) {
  library("groundhog")
} else {
  install.packages("groundhog")
  library("groundhog")
}

# define date of reference to load all packages
groundhog.date <- "2022-04-01"

# guarantee version 1.5 of groundhog is being used
groundhog::meta.groundhog(date = "2022-04-01")


# HERE
groundhog::groundhog.library("here", groundhog.date) # load package here


# TICTOC
groundhog::groundhog.library("tictoc", groundhog.date) # load package tictoc


# DECLARE LOCATION OF CURRENT SCRIPT TO SET UP PROJECT ROOT CORRECTLY
here::i_am("code/projectSpecific/8SitesModel/calibration_projectSpecific_8SitesModel.R", uuid = "0fef2915-0a76-4b7f-9a2c-c0b9830c058b")


# START TIME
tictoc::tic(msg = "calibration_projectSpecific_8SitesModel script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary
groundhog::groundhog.library("sf", groundhog.date)  # manipulate spatial data (vector format)
groundhog::groundhog.library("terra", groundhog.date)  # to calculate probability transition matrix
groundhog::groundhog.library("factoextra", groundhog.date)  # clustering
groundhog::groundhog.library("fixest", groundhog.date)  # regression


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# 1000 SITES MODEL CALIBRATION VARIABLES
load(here::here("data/projectSpecific/1000SitesModel", "calibration_1000SitesModel.Rdata"))


# CALIBRATION GLOBAL MODEL
load(here::here("data/projectSpecific/globalModel/calibration_globalModel.Rdata"))

# clean enviroment
rm(matrixTransition.2prices, matrixTransition.2prices.alt)


# MUNI LEVEL SPATIAL SAMPLE
load(here::here("data/projectSpecific/prepData/sampleMuniSpatial_prepData.Rdata"))






# CLUSTERING -----------------------------------------------------------------------------------------------------------------------------------------

# GROUP 1000 SITES INTO 8 CLUSTERS

# select variables of interest (that vary across sites)
aux.clusterData <-
  calibration.1000SitesModel %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(z_1985_1000Sites = 100*z_1985_1000Sites/zbar_1985_1000Sites,
                z_2017_1000Sites = 100*z_2017_1000Sites/zbar_2017_1000Sites,
                lon = sf::st_coordinates(sf::st_centroid(calibration.1000SitesModel))[,"X"],
                lat = sf::st_coordinates(sf::st_centroid(calibration.1000SitesModel))[,"Y"]) %>%
  dplyr::select(gamma_1000Sites, theta_1000Sites)



# STANDARDIZE VARIABLES
aux.clusterData <- scale(aux.clusterData)

# set.ssed to reproduce clustering numbers
set.seed(970)

# apply k-means clustering
aux.cluster <- kmeans(aux.clusterData, centers = 8, nstart = 25)

# inspect clusters
str(aux.cluster)
aux.cluster

# plot clusters
fviz_cluster(aux.cluster, data = aux.clusterData, pointsize = 3, stand = F, geom = "text", show.clust.cent = F) +
  ggplot2::scale_fill_brewer(name = "Cluster", palette = "Set3") +
  ggplot2::scale_color_brewer(name = "Cluster", palette = "Set3") +
  ggplot2::theme_classic(base_size = 20) +
  ggplot2::theme(legend.position = "bottom", legend.text = element_text(size = 16, face = "bold"), legend.margin=margin(t=-0.5, r=0, b=0.2, l=0, unit="cm")) +
  guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
ggplot2::ggsave(filename = here::here("data/analysis/8SitesModel/cluster_gammaTheta_1000Sites.png"), width = 10, height = 6)

# check optimal number of clusters
fviz_nbclust(aux.clusterData, kmeans, method = "wss", k.max = 15)
fviz_nbclust(aux.clusterData, kmeans, method = "silhouette", k.max = 15)

# add clusters to original 1000 sites data
calibration.1000SitesModel <-
  calibration.1000SitesModel %>%
  mutate(cluster = aux.cluster$cluster)

# add cluster numbering by gamma rank
calibration.1000SitesModel <-
  calibration.1000SitesModel %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(gamma_mean = mean(gamma_1000Sites)) %>%
  dplyr::arrange(gamma_mean) %>%
  dplyr::mutate(clusterRank = 1:nrow(.)) %>%
  dplyr::select(cluster, clusterRank) %>%
  dplyr::right_join(calibration.1000SitesModel) %>%
  sf::st_as_sf()

ggplot2::ggplot(data = calibration.1000SitesModel) +
  ggplot2::geom_sf(aes(fill = as.factor(cluster))) +
  ggplot2::scale_fill_brewer(name = "Cluster", palette = "Set3") +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(panel.grid.major = element_line(colour = "White"),
                 panel.grid.minor = element_line(colour = "white"),
                 panel.background = element_blank(),
                 strip.background = element_rect(fill = NA),
                 axis.line = element_blank(), axis.ticks = element_blank(),
                 axis.title = element_blank(), axis.text = element_blank(),
                 legend.title = element_text(hjust = 0.5, size = 22, face = "bold"),
                 legend.position = "bottom", legend.key.width = unit(1, "cm"), legend.margin=margin(t=-0.5, r=0, b=0.2, l=0, unit="cm"),
                 legend.text = element_text(size = 20, face = "bold"))

ggplot2::ggsave(filename = here::here("data/analysis/8SitesModel/map_8clusters_1000Sites.png"), width = 8, height = 6)


ggplot2::ggplot(data = calibration.1000SitesModel) +
  ggplot2::geom_sf(aes(fill = as.factor(clusterRank))) +
  ggplot2::scale_fill_viridis_d(name = "Cluster") +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(panel.grid.major = element_line(colour = "White"),
                 panel.grid.minor = element_line(colour = "white"),
                 panel.background = element_blank(),
                 strip.background = element_rect(fill = NA),
                 axis.line = element_blank(), axis.ticks = element_blank(),
                 axis.title = element_blank(), axis.text = element_blank(),
                 legend.title = element_text(hjust = 0.5, size = 22, face = "bold"),
                 legend.position = "bottom", legend.key.width = unit(1, "cm"), legend.margin=margin(t=-0.5, r=0, b=0.2, l=0, unit="cm"),
                 legend.text = element_text(size = 20, face = "bold"))

ggplot2::ggsave(filename = here::here("data/analysis/8SitesModel/map_8clustersRankGamma_1000Sites.png"), width = 8, height = 6)

# intermediate export before aggregating
save(calibration.1000SitesModel,
     file = here::here("data/projectSpecific/1000SitesModel",
                       "cluster8Sites_1000SitesModel.Rdata"))




# INITIAL CONDITIONS Z -------------------------------------------------------------------------------------------------------------------------------

# AGGREGATE FROM 1000 SITES TO 8 SITES
calibration.8SitesModel <-
  calibration.1000SitesModel %>%
  dplyr::select(cluster, clusterRank, amazonBiomeArea_ha_1000Sites,
                z_1985_1000Sites, forestArea_1985_ha_1000Sites, otherArea_1985_ha_1000Sites, otherArea_1985_ha_1000Sites, zbar_1985_1000Sites,
                z_2017_1000Sites, forestArea_2017_ha_1000Sites, otherArea_2017_ha_1000Sites, otherArea_2017_ha_1000Sites, zbar_2017_1000Sites) %>%
  dplyr::group_by(cluster, clusterRank) %>%
  dplyr::summarise(geometry = sf::st_union(geometry),
                   across(where(is.numeric), sum)) %>%
  dplyr::ungroup() %>%
  sf::st_as_sf() %>%
  dplyr::rename_with(.fn = ~str_replace(., "1000Sites", "8Sites"), .cols = ends_with("1000Sites"))

# add id variable
calibration.8SitesModel$id <- 1:nrow(calibration.8SitesModel)





# PARAMETER GAMMA ------------------------------------------------------------------------------------------------------------------------------------

# # DATA INPUT
# # load pixel sample with biomass data
# load(here::here("data/projectSpecific/prepData/pixelBiomass_prepData.Rdata"))
#
# # select minicells of primary forest with co2 information and transform to spatial points
# pixelBiomass.prepData <-
#   pixelBiomass.prepData %>%
#   dplyr::filter(mapbiomas_classAgg == "primaryForest") %>%
#   sf::st_transform(sf::st_crs(calibration.8SitesModel)) %>%
#   dplyr::mutate(co2e_ha = (agb/2)*(44/12)) %>%
#   dplyr::select(co2e_ha)
#
# # match minicells with sites
# aux.gamma <- sf::st_join(calibration.8SitesModel, pixelBiomass.prepData)
#
# # calculate average carbon density on primary forest areas by site
# aux.gamma <-
#   aux.gamma %>%
#   sf::st_drop_geometry() %>%
#   dplyr::group_by(id) %>%
#   summarise(gamma_8Sites = mean(co2e_ha, na.rm = T))
#
#
# # add gamma_8Sites to spatial variables
# calibration.8SitesModel <- left_join(calibration.8SitesModel, aux.gamma)
#
# # clean environment
# rm(pixelBiomass.prepData, aux.gamma, raster.variables)
#
# # add alternative value for gamma with a 10.1% correction factor for small trees and lianas (see Malhi et al (2006))
# calibration.8SitesModel <-
#   calibration.8SitesModel %>%
#   dplyr::mutate(gamma_alt_8Sites = gamma_8Sites*1.101)


# DATA INPUT (2017)
# load pixel sample with biomass data
load(here::here("data/projectSpecific/prepData/pixelBiomass2017_prepData.Rdata"))

# select minicells of primary forest with co2 information and transform to spatial points
pixelBiomass2017.prepData <-
  pixelBiomass2017.prepData %>%
  dplyr::filter(mapbiomas_classAgg == "primaryForest") %>%
  sf::st_transform(sf::st_crs(calibration.8SitesModel)) %>%
  dplyr::mutate(co2e_ha_2017 = (agb_2017/2)*(44/12),
                co2eSD_ha_2017 = (agbSD_2017/2)*(44/12)) %>%
  dplyr::select(lon, lat, co2e_ha_2017, co2eSD_ha_2017) %>%
  dplyr::filter(co2e_ha_2017 > 0, co2eSD_ha_2017 > 0, !is.na(co2e_ha_2017), !is.na(co2eSD_ha_2017))

# cattle value per ha
reg.carbonUncertainty <-
  pixelBiomass2017.prepData %>%
  lm(formula = log(co2eSD_ha_2017) ~ log(co2e_ha_2017) + lon*lat + I(lon^2) + I(lat^2))

# regression results
summary(reg.carbonUncertainty)

pixelBiomass2017.prepData %>%
  fixest::feols(log(co2eSD_ha_2017) ~ log(co2e_ha_2017) + lon*lat + I(lon^2) + I(lat^2)) %>%
  fixest::etable(se = "hetero", replace = T, file = here::here("data/analysis/8SitesModel/reg_gamma.tex"), style.tex = fixest::style.tex(main = "qje"))


# extract fitted values
pixelBiomass2017.prepData <-
  pixelBiomass2017.prepData %>%
  dplyr::mutate(co2eSD_ha_2017_fitted = predict(reg.carbonUncertainty, .),
                co2eSD_ha_2017_fitted = exp(co2eSD_ha_2017_fitted))

# match minicells with sites
aux.gamma2017 <- sf::st_join(calibration.8SitesModel, pixelBiomass2017.prepData)

# calculate average carbon density on primary forest areas by site
aux.gamma2017 <- aux.gamma2017 %>% sf::st_drop_geometry() %>% dplyr::group_by(id) %>% summarise(gamma_8Sites = mean(co2e_ha_2017, na.rm = T),
                                                                                                co2e_ha_2017 = mean(co2e_ha_2017, na.rm = T),
                                                                                                gammaSD_8Sites = mean(co2eSD_ha_2017_fitted, na.rm = T),
                                                                                                co2eSD_ha_2017 = mean(co2eSD_ha_2017, na.rm = T),
                                                                                                lon = mean(lon, na.rm = T),
                                                                                                lat = mean(lat, na.rm = T))

# add gamma_8Sites to spatial variables
calibration.8SitesModel <- left_join(calibration.8SitesModel, aux.gamma2017)

# clean environment
rm(pixelBiomass2017.prepData, aux.gamma2017)

calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(gammaSD_8Sites = predict(reg.carbonUncertainty, .),
                gammaSD_8Sites = exp(gammaSD_8Sites))

# clean environment
rm(reg.carbonUncertainty)





# PARAMETERS A AND B ---------------------------------------------------------------------------------------------------------------------------------

# estimate of a same as in the global model >
# estimate of b given a such that b_i/a is the average carbon in primary forest
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(a_8Sites = calibration.globalModel$a_global,
                b_8Sites = a_8Sites*gamma_8Sites*zbar_2017_8Sites)





# PARAMETER K ----------------------------------------------------------------------------------------------------------------------------------------

# estimate of k same as in the global model
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(k_8Sites = calibration.globalModel$k_global)





# PARAMETER ZETA ----------------------------------------------------------------------------------------------------------------------------------------

# estimate of zeta same as in the global model
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(zeta_8Sites = calibration.globalModel$zeta_global,
                zeta_alt_8Sites = calibration.globalModel$zeta_alt_global)




# PARAMETER RHO ----------------------------------------------------------------------------------------------------------------------------------------

# estimate of rho same as in the global model
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(rho_8Sites = calibration.globalModel$rho_global)





# PARAMETER R ----------------------------------------------------------------------------------------------------------------------------------------

# estimate of r same as in the global model
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(r_8Sites = calibration.globalModel$r_global)





# INITIAL CONDITIONS X -------------------------------------------------------------------------------------------------------------------------------

# x_2017_8Sites estimated as in the old way of global model, just considering the stock of carbon stored in forest areas assuming that all forests are primary
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(x_2017_8Sites = gamma_8Sites*forestArea_2017_ha_8Sites)





# INITIAL CONDITIONS C -------------------------------------------------------------------------------------------------------------------------------

# estimates of C same as in the global model
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(C_2017_8Sites = calibration.globalModel$C_2017_global,
                C0bern_2017_8Sites = calibration.globalModel$C0bern_2017_global,
                C1bern_2017_8Sites = calibration.globalModel$C1bern_2017_global,
                C2bern_2017_8Sites = calibration.globalModel$C2bern_2017_global,
                C3bern_2017_8Sites = calibration.globalModel$C3bern_2017_global)




# PARAMETER U ----------------------------------------------------------------------------------------------------------------------------------------

# estimates of u same as in the global model
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(u_8Sites = calibration.globalModel$u_global,
                u_bern_8Sites = calibration.globalModel$u_bern_global)





# PARAMETER THETA ------------------------------------------------------------------------------------------------------------------------------------

# DATA INPUT
# load variables at the muni level to calibrate theta
load("data/projectSpecific/prepData/muniTheta_prepData.Rdata")

# load cattle price series
load("data/projectSpecific/prepData/seriesPriceCattle_prepData.Rdata")


# DATA MANIPULATION

# EXTRACT AVERAGE 2017 PRICE (use real prices because it is normalized to 2017 )
aux.price.2017 <-
  seriesPriceCattle.prepData %>%
  dplyr::mutate(price_real_mon_cattle = price_real_mon_cattle/3.1966) %>%
  dplyr::filter(year == 2017) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(mean_price_2017 = mean(price_real_mon_cattle)) %>%
  dplyr::pull(mean_price_2017)

# EXTRACT AVERAGE 2006 PRICE (use real prices because it is normalized to 2006)
aux.price.2006 <-
  seriesPriceCattle.prepData %>%
  dplyr::mutate(price_real_mon_cattle = price_real2006_mon_cattle/2.1761) %>%
  dplyr::filter(year == 2006) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(mean_price_2006 = mean(price_real_mon_cattle)) %>%
  dplyr::pull(mean_price_2006)


# REGRESSION - CATTLE VALUE (2017)

# cattle value per ha
reg.cattleValueperHa <-
  muniTheta.prepData %>%
  lm(formula = cattleSlaughter_value_ha ~ pastureArea_value + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
       lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = pastureArea_value)

# regression results
summary(reg.cattleValueperHa)

# extract fitted values
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(cattleSlaughter_value_ha_fitted = predict(reg.cattleValueperHa, .))

# extract minimum positive fitted value
aux.min.positive.cattleSlaughter.value.ha.fitted <-
  muniTheta.prepData %>%
  dplyr::filter(cattleSlaughter_value_ha_fitted > 0) %>%
  dplyr::pull(cattleSlaughter_value_ha_fitted) %>%
  min()

# winsorize cattleSlaughter_value_ha_fitted
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(d_theta_winsorized = dplyr::if_else(cattleSlaughter_value_ha_fitted <= 0,
                                                    1,
                                                    0),
                cattleSlaughter_value_ha_fitted = dplyr::if_else(cattleSlaughter_value_ha_fitted <= 0,
                                                                 aux.min.positive.cattleSlaughter.value.ha.fitted,
                                                                 cattleSlaughter_value_ha_fitted))
# clean environment
rm(reg.cattleValueperHa, aux.min.positive.cattleSlaughter.value.ha.fitted)


# REGRESSION - CATTLE VALUE (2006)

# cattle value per ha
reg.cattleValueperHa <-
  muniTheta.prepData %>%
  lm(formula = cattleSlaughter2006_value_ha ~ pastureArea2006_value + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
       lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = pastureArea2006_value)

# regression results
summary(reg.cattleValueperHa)

# extract fitted values
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(cattleSlaughter2006_value_ha_fitted = predict(reg.cattleValueperHa, .))

# extract minimum positive fitted value
aux.min.positive.cattleSlaughter.value.ha.fitted <-
  muniTheta.prepData %>%
  dplyr::filter(cattleSlaughter2006_value_ha_fitted > 0) %>%
  dplyr::pull(cattleSlaughter2006_value_ha_fitted) %>%
  min()

# winsorize cattleSlaughter_value_ha_fitted
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(d_theta2006_winsorized = dplyr::if_else(cattleSlaughter2006_value_ha_fitted <= 0,
                                                        1,
                                                        0),
                cattleSlaughter2006_value_ha_fitted = dplyr::if_else(cattleSlaughter2006_value_ha_fitted <= 0,
                                                                     aux.min.positive.cattleSlaughter.value.ha.fitted,
                                                                     cattleSlaughter2006_value_ha_fitted))
# clean environment
rm(reg.cattleValueperHa, aux.min.positive.cattleSlaughter.value.ha.fitted)


# REGRESSION - FARM GATE PRICE (2017)

# cattle value per ha
reg.cattleFarmGatePricePerArroba <-
  muniTheta.prepData %>%
  lm(formula = cattleSlaughter_farmGatePrice ~ cattleSlaughter_head + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
       lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = cattleSlaughter_head)

# regression results
summary(reg.cattleFarmGatePricePerArroba)

# extract fitted values
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(cattleSlaughter_farmGatePrice_fitted = predict(reg.cattleFarmGatePricePerArroba, .))

# clean environment
rm(reg.cattleFarmGatePricePerArroba)



# REGRESSION - FARM GATE PRICE (2006)

# cattle value per ha
reg.cattleFarmGatePricePerArroba <-
  muniTheta.prepData %>%
  lm(formula = cattleSlaughter2006_farmGatePrice ~ cattleSlaughter2006_head + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
       lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = cattleSlaughter2006_head)

# regression results
summary(reg.cattleFarmGatePricePerArroba)

# extract fitted values
muniTheta.prepData <-
  muniTheta.prepData %>%
  dplyr::mutate(cattleSlaughter2006_farmGatePrice_fitted = predict(reg.cattleFarmGatePricePerArroba, .))

# clean environment
rm(reg.cattleFarmGatePricePerArroba)


# change crs to match with calibration.8SitesModel and select variables
muniTheta.prepData <-
  muniTheta.prepData %>%
  sf::st_transform(crs = sf::st_crs(calibration.8SitesModel)) %>%
  dplyr::select(muni_code, muni_area, cattleSlaughter_value_ha_fitted, cattleSlaughter_farmGatePrice_fitted, pastureArea_value, d_theta_winsorized,
                cattleSlaughter2006_value_ha_fitted, cattleSlaughter2006_farmGatePrice_fitted, pastureArea2006_value, d_theta2006_winsorized)

# match munis with sites
aux.theta <- sf::st_intersection(calibration.8SitesModel, muniTheta.prepData)

# calculate muni areas inside each site
aux.theta$muni_site_area <-
  sf::st_area(aux.theta) %>%
  units::set_units(ha) %>%
  unclass()

# calculate cattleSlaughter_value_ha_fitted and pastureArea_value by site (for each muni adjust the value by the share of the muni area inside the site)
aux.theta.2017 <-
  aux.theta %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(cattleSlaughter_value_ha_fitted = weighted.mean(cattleSlaughter_value_ha_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter_value_ha_scenario2 = weighted.mean((74.7/15)*cattleSlaughter_farmGatePrice_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter_value_ha_scenario3 = weighted.mean((140.2/15)*cattleSlaughter_farmGatePrice_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter_value_ha_scenario4 = weighted.mean((201.7/15)*cattleSlaughter_farmGatePrice_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter_value_ha_scenario5 = weighted.mean((221.4/15)*cattleSlaughter_farmGatePrice_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   pastureArea_value = sum(pastureArea_value*(muni_site_area/muni_area), na.rm = T),
                   d_theta_winsorized = min(d_theta_winsorized, na.rm = T))

# add cattleSlaughter_value_ha_fitted and pastureArea_value to spatial variables
calibration.8SitesModel <- left_join(calibration.8SitesModel, aux.theta.2017)

aux.theta.2006 <-
  aux.theta %>%
  sf::st_drop_geometry() %>%
  dplyr::filter(!is.na(pastureArea2006_value)) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(cattleSlaughter2006_value_ha_fitted = weighted.mean(cattleSlaughter2006_value_ha_fitted, w = pastureArea2006_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter2006_value_ha_scenario2 = weighted.mean((74.7/15)*cattleSlaughter2006_farmGatePrice_fitted, w = pastureArea2006_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter2006_value_ha_scenario3 = weighted.mean((140.2/15)*cattleSlaughter2006_farmGatePrice_fitted, w = pastureArea2006_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter2006_value_ha_scenario4 = weighted.mean((201.7/15)*cattleSlaughter2006_farmGatePrice_fitted, w = pastureArea2006_value*(muni_site_area/muni_area), na.rm = T),
                   cattleSlaughter2006_value_ha_scenario5 = weighted.mean((221.4/15)*cattleSlaughter2006_farmGatePrice_fitted, w = pastureArea2006_value*(muni_site_area/muni_area), na.rm = T),
                   pastureArea2006_value = sum(pastureArea2006_value*(muni_site_area/muni_area), na.rm = T),
                   d_theta2006_winsorized = min(d_theta2006_winsorized, na.rm = T))

# add cattleSlaughter_value_ha_fitted and pastureArea_value to spatial variables
calibration.8SitesModel <- left_join(calibration.8SitesModel, aux.theta.2006)

# clean environment
rm(muniTheta.prepData, aux.theta)

# calculate theta_8Sites
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(theta_8Sites = cattleSlaughter_value_ha_fitted/(aux.price.2017),
                theta_scenario2_8Sites = cattleSlaughter_value_ha_scenario2/(aux.price.2017),
                theta_scenario3_8Sites = cattleSlaughter_value_ha_scenario3/(aux.price.2017),
                theta_scenario4_8Sites = cattleSlaughter_value_ha_scenario4/(aux.price.2017),
                theta_scenario5_8Sites = cattleSlaughter_value_ha_scenario5/(aux.price.2017)) %>%
  dplyr::mutate(theta2006_8Sites = cattleSlaughter2006_value_ha_fitted/(aux.price.2006),
                theta2006_scenario2_8Sites = cattleSlaughter2006_value_ha_scenario2/(aux.price.2006),
                theta2006_scenario3_8Sites = cattleSlaughter2006_value_ha_scenario3/(aux.price.2006),
                theta2006_scenario4_8Sites = cattleSlaughter2006_value_ha_scenario4/(aux.price.2006),
                theta2006_scenario5_8Sites = cattleSlaughter2006_value_ha_scenario5/(aux.price.2006)) %>%
  dplyr::select(-cattleSlaughter_value_ha_fitted, -cattleSlaughter_value_ha_scenario2, -cattleSlaughter_value_ha_scenario3,
                -cattleSlaughter_value_ha_scenario4, -cattleSlaughter_value_ha_scenario5,
                -cattleSlaughter2006_value_ha_fitted, -cattleSlaughter2006_value_ha_scenario2, -cattleSlaughter2006_value_ha_scenario3,
                -cattleSlaughter2006_value_ha_scenario4, -cattleSlaughter2006_value_ha_scenario5)


# identify adjacent neighbors
aux.neighbors <- sf::st_is_within_distance(calibration.8SitesModel, calibration.8SitesModel, dist = 100, remove_self = TRUE)

# impute values for missing thetas based on the average of adjacent neighbors
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(d_imputation_theta_scenario2_8Sites = dplyr::if_else(is.na(theta_scenario2_8Sites), 1, 0),
                theta_scenario2_8Sites = dplyr::if_else(is.na(theta_scenario2_8Sites),
                                                          apply(aux.neighbors, 1, function(i){mean(.$theta_scenario2_8Sites[i], na.rm = TRUE)}),
                                                          theta_scenario2_8Sites),
                d_imputation_theta_scenario3_8Sites = dplyr::if_else(is.na(theta_scenario3_8Sites), 1, 0),
                theta_scenario3_8Sites = dplyr::if_else(is.na(theta_scenario3_8Sites),
                                                          apply(aux.neighbors, 1, function(i){mean(.$theta_scenario3_8Sites[i], na.rm = TRUE)}),
                                                          theta_scenario3_8Sites),
                d_imputation_theta_scenario4_8Sites = dplyr::if_else(is.na(theta_scenario4_8Sites), 1, 0),
                theta_scenario4_8Sites = dplyr::if_else(is.na(theta_scenario4_8Sites),
                                                          apply(aux.neighbors, 1, function(i){mean(.$theta_scenario4_8Sites[i], na.rm = TRUE)}),
                                                          theta_scenario4_8Sites),
                d_imputation_theta_scenario5_8Sites = dplyr::if_else(is.na(theta_scenario5_8Sites), 1, 0),
                theta_scenario5_8Sites = dplyr::if_else(is.na(theta_scenario5_8Sites),
                                                          apply(aux.neighbors, 1, function(i){mean(.$theta_scenario5_8Sites[i], na.rm = TRUE)}),
                                                          theta_scenario5_8Sites))

# THETA (SOYBEAN)

# DATA INPUT
# load soybean potential yield raster
raster.soybean <- terra::rast(here::here("data/raw2clean/potentialYield_faogaez/output/clean_potentialYield.tif"))



# STORE PARAMETER VALUES
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::bind_cols(dplyr::tibble(theta_soybean_8Sites = terra::extract(raster.soybean, terra::vect(calibration.8SitesModel %>% sf::st_transform(4326)), mean, na.rm = T)$soyb200b_yld))


# clean environment
rm(raster.soybean)





# EMPLOYMENT LOSS ------------------------------------------------------------------------------------------------------------------------------------

# DATA INPUT
# load variables at the muni level to calibrate theta
load(here::here("data/projectSpecific/prepData/muniAgCensusCattleRaising_prepData.Rdata"))



# DATA MANIPULATION

# REGRESSION - CATTLE VALUE

# cattle value per ha
reg.workersperHa <-
  muniAgCensusCattleRaising.prepData %>%
  lm(formula = workers_ha ~ pastureArea_value + historical_precip + I(historical_precip^2) + historical_temp + I(historical_temp^2) +
       lon*lat + I(lon^2) + I(lat^2), na.action = na.exclude, weights = pastureArea_value)

# regression results
summary(reg.workersperHa)

# extract fitted values
muniAgCensusCattleRaising.prepData <-
  muniAgCensusCattleRaising.prepData %>%
  dplyr::mutate(workers_ha_fitted = predict(reg.workersperHa, .))

# extract minimum positive fitted value
aux.min.positive.workers.ha.fitted <-
  muniAgCensusCattleRaising.prepData %>%
  dplyr::filter(workers_ha_fitted > 0) %>%
  dplyr::pull(workers_ha_fitted) %>%
  min()

# windsorize workers_ha_fitted
muniAgCensusCattleRaising.prepData <-
  muniAgCensusCattleRaising.prepData %>%
  dplyr::mutate(workers_ha_fitted = dplyr::if_else(workers_ha_fitted <= 0,
                                                   aux.min.positive.workers.ha.fitted,
                                                   workers_ha_fitted))

# clean environment
rm(reg.workersperHa, aux.min.positive.workers.ha.fitted)


# change crs to match with calibration.8SitesModel and select variables
muniAgCensusCattleRaising.prepData <-
  muniAgCensusCattleRaising.prepData %>%
  sf::st_transform(crs = sf::st_crs(calibration.8SitesModel)) %>%
  dplyr::select(muni_code, muni_area, workers_ha_fitted, pastureArea_value)

# match munis with sites
aux.worker <- sf::st_intersection(calibration.8SitesModel, muniAgCensusCattleRaising.prepData)

# calculate muni areas inside each site
aux.worker$muni_site_area <-
  sf::st_area(aux.worker) %>%
  units::set_units(ha) %>%
  unclass()

# calculate workers_ha_fitted and pastureArea_value by site (for each muni adjust the value by the share of the muni area inside the site)
aux.worker <-
  aux.worker %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(workers_ha_fitted = weighted.mean(workers_ha_fitted, w = pastureArea_value*(muni_site_area/muni_area), na.rm = T))

# add workers_ha_fitted and pastureArea_value to spatial variables
calibration.8SitesModel <- left_join(calibration.8SitesModel, aux.worker)

# clean environment
rm(muniAgCensusCattleRaising.prepData, aux.worker)

# calculate theta_8Sites
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(workerPerHa_8Sites = workers_ha_fitted) %>%
  dplyr::select(-workers_ha_fitted)




# PRICE INITIAL CONDITION ----------------------------------------------------------------------------------------------------------------------------

# estimates of p_2017 same as in the global model
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::mutate(p_2017_8Sites = calibration.globalModel$p_2017_global,
                p_2017_alt_8Sites = calibration.globalModel$p_2017_alt_global,
                pSoybean_2017_8Sites = calibration.globalModel$pSoybean_2017_global,
                pSoybean_2017_alt_8Sites = calibration.globalModel$pSoybean_2017_alt_global)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# REMOVE NAs
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::filter(!is.na(gamma_8Sites))


# ORDER VARIABLES
calibration.8SitesModel <-
  calibration.8SitesModel %>%
  dplyr::select(id, ends_with("ha_8Sites"), cluster, clusterRank,
                z_2017_8Sites, zbar_2017_8Sites, x_2017_8Sites, b_8Sites, theta_8Sites, gamma_8Sites, gammaSD_8Sites,
                d_theta_winsorized,
                pastureArea_value, ends_with("_8Sites"))


# POST-TREATMENT OVERVIEW
# summary(calibration.8SitesModel)
# View(calibration.8SitesModel)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(calibration.8SitesModel,
     file = here::here("data/projectSpecific/8SitesModel",
                       "calibration_8SitesModel.Rdata"))

# remove spatial feature
calibration.8SitesModel <- calibration.8SitesModel %>% sf::st_drop_geometry()

readr::write_csv(calibration.8SitesModel,
                 file = here::here("data/projectSpecific/8SitesModel", "calibration_8SitesModel.csv"))


# CLEAN TEMP DIR
terra::tmpFiles(current = T, remove = T)
gc()



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/8SitesModel")






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------