# Start timer
tictoc::tic(msg = "_masterfile.R script", log = TRUE)

source("rsrc/raw2clean/_masterfile_raw2clean.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

print("Data cleaning is done")

source("rsrc/processing/_masterfile_prep.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

print("Data processing is done")

source("rsrc/calibration/_masterfile_calibration.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

print("Calibration is done")

# End timer
tictoc::toc(log = TRUE)
