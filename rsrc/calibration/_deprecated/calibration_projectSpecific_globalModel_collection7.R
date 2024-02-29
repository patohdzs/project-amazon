
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (GLOBAL MODEL)
# AUTHOR: JOÃO PEDRO VIEIRA
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
here::i_am("code/calibration/globalModel/calibration_projectSpecific_globalModel.R", uuid = "59c043ea-fe79-4129-a9ad-144c141ef9a2")


# START TIME
tictoc::tic(msg = "calibration_projectSpecific_globalModel script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary
groundhog::groundhog.library("sf", groundhog.date)  # manipulate spatial data (vector format)
groundhog::groundhog.library("terra", groundhog.date)  # manipulate spatial data (raster format)
groundhog::groundhog.library("markovchain", groundhog.date)  # to calculate probability transition matrix
groundhog::groundhog.library("equatiomatic", groundhog.date)  # to convert lm models to latex




# INITIAL CONDITIONS Z -------------------------------------------------------------------------------------------------------------------------------

# Z_BAR (MAXIMUM VALUE OF Z = SUM OF FOREST + AGRICULTURAL USE AREAS IN 2017)
# load pixel sample with mapbiomas categories data
load(here::here("data/raw2clean/landUseCoverBiome_mapbiomas/output/clean_landUseCoverBiome.Rdata"))

# calculate maximum value of z (areas of forest and agricultural use in 2017 in the Amazon Biome)
zbar2017 <-
  clean.landUseCoverBiome %>%
  dplyr::filter(year == 2017) %>%
  dplyr::mutate(zbar = z + forest) %>%
  dplyr::pull(zbar)

# calculate z0 (area of agricultural use in 2017 in the Amazon Biome)
z2017 <-
  clean.landUseCoverBiome %>%
  dplyr::filter(year == 2017) %>%
  dplyr::pull(z)

calibration.globalModel <- dplyr::tibble(zbar_2017_global = zbar2017,
                                         z_2017_global = z2017)


# PARAMETER GAMMA ------------------------------------------------------------------------------------------------------------------------------------

# DATA INPUT
# load pixel sample with biomass data
load(here::here("data/calibration/prepData/pixelBiomass_prepData.Rdata"))


# DATA MANIPULATION
# drop spatial feature and change unit of measure from biomass per hectare to Carbon per hectare (50% of biomass)
pixelBiomass.prepData <-
  pixelBiomass.prepData %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(co2e_ha = (agb/2)*(44/12))

# average value of co2e_ha in primary forests
gamma_co2e_ha <- pixelBiomass.prepData %>% dplyr::filter(mapbiomas_classAgg == "primaryForest") %>% pull(co2e_ha) %>% mean(na.rm = T)


# STORE PARAMETER VALUES
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::bind_cols(dplyr::tibble(gamma_global = gamma_co2e_ha,
                                 gamma_alt_global = gamma_co2e_ha*1.101)) # alternative value for gamma with a 10.1% correction factor for small trees and lianas (see Malhi et al (2006))
                                 #gamma_up_global = gamma_global*1.25, # upper value based on 25% uncertainty from Malhi et al (2006)
                                 #gamma_low_global = gamma_global*0.75, # lower value based on 25% uncertainty from Malhi et al (2006)
                                 #gamma_alt_up_global = gamma_alt_global*1.25, # upper alternative value based on 25% uncertainty from Malhi et al (2006)
                                 #gamma_alt_low_global = gamma_alt_global*0.75)) # lower alternative value based on 25% uncertainty from Malhi et al (2006)

# clean environment
rm(gamma_co2e_ha)





# PARAMETERS A AND B ---------------------------------------------------------------------------------------------------------------------------------

# DATA MANIPULATION
# estimate of a to make convergence (0.99*b/a) happens in 100 years (time based on Heinrich et al (2021))
a <- 1 - (1-0.99)^(1/100)

# estimate of b given a such that b/a is the average carbon in primary forest
b <- calibration.globalModel$gamma_global*a*calibration.globalModel$zbar_2017_global


# STORE PARAMETER VALUES
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::bind_cols(dplyr::tibble(a_global = a,
                                 b_global = b,
                                 b_alt_global = calibration.globalModel$gamma_alt_global*a*calibration.globalModel$zbar_2017_global))
                                 #b_up_global = calibration.globalModel$gamma_up_global*a*calibration.globalModel$zbar_2017_global,
                                 #b_low_global = calibration.globalModel$gamma_low_global*a*calibration.globalModel$zbar_2017_global,
                                 #b_alt_up_global = calibration.globalModel$gamma_alt_up_global*a*calibration.globalModel$zbar_2017_global,
                                 #b_alt_low_global = calibration.globalModel$gamma_alt_low_global*a*calibration.globalModel$zbar_2017_global))

# clean environment
rm(a, b)





# PARAMETER K ----------------------------------------------------------------------------------------------------------------------------------------

# DATA INPUT
# load pixel sample with biomass data
load(here::here("data/calibration/prepData/stateEmission_prepData.Rdata"))


# DATA MANIPULATION
# calculate average of emission factor from agricultural use across years and states
avg.emissionFactor <-
  stateEmission.prepData %>%
  dplyr::ungroup() %>%
  dplyr::summarise(emissionFactor_co2e_ha = sum(emission_co2e)/sum(agriculturalUse_area)) %>%
  pull(emissionFactor_co2e_ha)

# calculate average of net emission factor from agricultural use across years and states
avg.netEmissionFactor <-
  stateEmission.prepData %>%
  dplyr::ungroup() %>%
  dplyr::summarise(netEmissionFactor_co2e_ha = sum(netEmission_co2e)/sum(agriculturalUse_area)) %>%
  pull(netEmissionFactor_co2e_ha)


# STORE PARAMETER VALUES
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::bind_cols(dplyr::tibble(k_global = avg.netEmissionFactor,
                                 k_alt_global = avg.emissionFactor))


# clean environment
rm(avg.emissionFactor, avg.netEmissionFactor)



# PARAMETER ZETA ----------------------------------------------------------------------------------------------------------------------------------------

# zeta is calibrated such that the marginal cost of changing land use (zeta*forestToPastureTransitionArea) matches the forest to pasture transition cost >
# estimated by Araujo, Costa and Sant'Anna (2022), reported on column 4 of table 4 (on the right). We transform to dollars using an FX rate of 4.14 (December 2019) >
# as in the paper: 1614.54/4.14 = 390 USD/ha. The paper also calculates that 6.5% of the forest area was converted to pasture between 2008 and 2017. >
aux.transitionCost <- 1614.54/4.14

# The forest area in 2008 represented 72% of the Legal Amazon area, which covers 501,506,775 ha, so the transition in hectares is >
#  0.065*0.72*501,506,775 = 23,470,517, resulting in an annual average of 2,347,052 ha.
aux.transitionArea <- (0.065*0.72*501506775)/(2017-2008+1)

zeta <- aux.transitionCost/aux.transitionArea
zeta_alt <- 483/aux.transitionArea # Alternative value based on a quote from (https://www.otempo.com.br/brasil/investigacoes-revelam-quadrilhas-e-ganho-milionario-por-tras-do-desmate-1.2229571)


# STORE PARAMETER VALUES
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::bind_cols(dplyr::tibble(zeta_global = zeta,
                                 zeta_alt_global = zeta_alt))


# clean environment
rm(zeta, zeta_alt, aux.transitionCost, aux.transitionArea)





# PARAMETER RHO ----------------------------------------------------------------------------------------------------------------------------------------

# rho was not calibrated, just set as 0.02

# STORE PARAMETER VALUES
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::bind_cols(dplyr::tibble(rho_global = 0.02))





# PARAMETER R ----------------------------------------------------------------------------------------------------------------------------------------

# r was defined as 0.025 per year decay rate based on http://euanmearns.com/the-half-life-of-co2-in-earths-atmosphere-part-1/

# STORE PARAMETER VALUES
calibration.globalModel <-
  calibration.globalModel %>%
  dplyr::bind_cols(dplyr::tibble(r_global = 0.025))





# INITIAL CONDITIONS X AND C -------------------------------------------------------------------------------------------------------------------------

# calculate z series
aux.zSeries <-
  clean.landUseCoverBiome %>%
  dplyr::bind_rows(tibble(year = 1972:1984, z = NA)) %>%
  dplyr::arrange(year) %>%
  dplyr::mutate(z_predicted = predict.lm(object = lm(data = ., z ~ year),
                                         newdata = data.frame(year = 1972:2021))) %>%
  dplyr::mutate(z_mixed = dplyr::case_when(year == 1972 ~ 0,
                                           year %in% 1973:1984 ~ z_predicted,
                                           TRUE ~ z)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(z_dot = z_mixed - lag(z_mixed)) # calculate forest annual changes (deforestation if positive, regeneration if negative)



# INITIAL CONDITION X

# generate variables to store values starting at b/a and 0
aux.zSeries$x_t_old <- calibration.globalModel$b_global/calibration.globalModel$a_global
aux.zSeries$x_t <- calibration.globalModel$b_global/calibration.globalModel$a_global
aux.zSeries$x_t_alt <- calibration.globalModel$b_alt_global/calibration.globalModel$a_global
aux.zSeries$e_t <- 0
aux.zSeries$e_t_alt <- 0
aux.zSeries$e_t_old <- 0

# calculate e_t and x_t
for (y in seq_along(1972:2020)) {

  aux.zSeries$x_t[y+1] <- aux.zSeries$x_t[y] - aux.zSeries$z_dot[y+1]*calibration.globalModel$gamma_global - aux.zSeries$x_t[y]*calibration.globalModel$a_global + calibration.globalModel$b_global*(1-aux.zSeries$z_mixed[y+1]/calibration.globalModel$zbar_2017_global)
  aux.zSeries$x_t_alt[y+1] <- aux.zSeries$x_t_alt[y] - aux.zSeries$z_dot[y+1]*calibration.globalModel$gamma_alt_global - aux.zSeries$x_t_alt[y]*calibration.globalModel$a_global + calibration.globalModel$b_alt_global*(1-aux.zSeries$z_mixed[y+1]/calibration.globalModel$zbar_2017_global)
  aux.zSeries$x_t_old[y+1] <- aux.zSeries$x_t_old[y] - aux.zSeries$z_dot[y+1]*calibration.globalModel$gamma_global - aux.zSeries$x_t_old[y]*calibration.globalModel$a_global + calibration.globalModel$b_global
  aux.zSeries$e_t[y+1] <- -(aux.zSeries$x_t[y+1]-aux.zSeries$x_t[y]) + aux.zSeries$z_mixed[y+1]*calibration.globalModel$k_global
  aux.zSeries$e_t_alt[y+1] <- -(aux.zSeries$x_t_alt[y+1]-aux.zSeries$x_t_alt[y]) + aux.zSeries$z_mixed[y+1]*calibration.globalModel$k_global
  aux.zSeries$e_t_old[y+1] <- -(aux.zSeries$x_t_old[y+1]-aux.zSeries$x_t_old[y]) + aux.zSeries$z_mixed[y+1]*calibration.globalModel$k_global

}



# SIMPLIFIED CARBON ACCUMULATION
aux.zSeries$C_t <- 0
aux.zSeries$C_t_alt <- 0
aux.zSeries$C_t_old <- 0

for (y in seq_along(1972:2020)) {

  aux.zSeries$C_t[y+1] <-  aux.zSeries$e_t[y+1] +aux.zSeries$C_t[y] -0.025*aux.zSeries$C_t[y]
  aux.zSeries$C_t_alt[y+1] <-  aux.zSeries$e_t_alt[y+1] +aux.zSeries$C_t_alt[y] -0.025*aux.zSeries$C_t_alt[y]
  aux.zSeries$C_t_old[y+1] <-  aux.zSeries$e_t_old[y+1] +aux.zSeries$C_t_old[y] -0.025*aux.zSeries$C_t_old[y]

}


# EXPORT TABLE WITH Z, X AND C VALUES FOR DIFFERENT YEARS (1973-2017)
variablesEvolution.globalModel <-
  aux.zSeries %>%
  dplyr::select(year, `z (ha)` = z_mixed, `x (CO2e Mg)` = x_t, `C (CO2e Mg)` = C_t) %>%
  dplyr::mutate(data_type = dplyr::if_else(year < 1985, "linear extrapolation", "observed"))

readr::write_csv(variablesEvolution.globalModel,
                 file = here::here("data/calibration/globalModel", "variablesEvolution_globalModel_collection7.csv"))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("calibration/globalModel")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
