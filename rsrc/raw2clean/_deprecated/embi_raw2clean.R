
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW EMBI+ BRAZIL RISK (JP MORGAN)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# -





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
here::i_am("code/raw2clean/embi_raw2clean.R", uuid = "0954ec36-3be8-435c-b236-3cf550c9224b")


# START TIME
tictoc::tic(msg = "embi_raw2clean script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read xls file
raw.embi <- readr::read_csv(file      = here::here("data/raw2clean/embi_jpmorgan/input/ipeadata[05-07-2022-07-14].csv"),
                              col_names = c("date", "embi", ""),
                              col_types = "cn-",
                              skip      = 1)



# DATA EXPLORATION
#summary(raw.embi)    # object is a list of data frames
#class(raw.embi)
#View(raw.embi)      # column names indicate file of origin





# DATASET CLEANUP AND PREP --------------------------------------------------------------------------------------------------------------------------

# transform to date class
raw.embi <-
  raw.embi %>%
  dplyr::mutate(date = lubridate::dmy(date, truncated = 1))

# remove weekends (NA)
raw.embi <-
  raw.embi %>%
  dplyr::filter(!is.na(embi))



# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(raw.embi$date)            <- "calendar date (yyyy-mm-dd), daily data"
sjlabelled::set_label(raw.embi$embi)    <- "EMBI+"


# change object name for exportation
clean.embi <- raw.embi



# POST-TREATMENT OVERVIEW
# summary(clean.embi)
# View(clean.embi)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.embi,
     file = here::here("data/raw2clean/embi_jpmorgan/output",
                      "clean_embi.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT: ------------------------------------------------------------------------------------------------------------------------------------