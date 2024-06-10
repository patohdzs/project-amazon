
# > THIS SCRIPT
# AIM: MASTERFILE SCRIPT TO ALL R data scripts



# START TIMER
tictoc::tic(msg = "_masterfile.R script", log = T)

source("rsrc/raw2clean/_masterfile_raw2clean.R", encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

print("raw2clean part is done")



source("rsrc/prepData/_masterfile_prep.R", encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

print("prepdata part is done")


source("rsrc/calibration/_masterfile_calibration.R", encoding = "UTF-8", echo = T)

# clear environment
rm(list = ls())

print("calibration part is done")

# END TIMER
tictoc::toc(log = T)



# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------