# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: MASTERFILE SCRIPT TO SOURCE ALL CALIBRATION SCRIPTS
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(tictoc)

# Start timer
tic(msg = "_masterfile.R script", log = TRUE)

# Calibrate model with granular grid
source("rsrc/calibration/calibrate_1043_sites_model.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Calibrate model with coarse grid
source("rsrc/calibration/calibrate_78_sites_model.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Prepare gamma regression data for HMC
source("rsrc/calibration/calibrate_gamma_reg.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Prepare theta regression data for HMC
source("rsrc/calibration/calibrate_theta_reg.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# End timer
toc(log = TRUE)
