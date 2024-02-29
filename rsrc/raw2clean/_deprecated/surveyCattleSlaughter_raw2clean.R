
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW CATTLE SLAUGHTER (HEAD COUNT AND WEIGHT) -IBGE
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
here::i_am("code/raw2clean/surveyCattleSlaughter_raw2clean.R", uuid = "eb3955a6-f041-436a-a011-19ff5d34d9aa")


# START TIME
tictoc::tic(msg = "surveyCattleSlaughter_raw2clean script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.cattle <-  readr::read_csv(file = here::here("data/raw2clean/surveyCattleSlaughter_ibge/input/surveyCattleSlaughter.csv"),
                               skip = 6,
                               na = c("...", "X"),
                               col_names = c("state_code", "state_name", "date", "cattleSlaughter_head", "cattleSlaughter_weight"),
                               col_types = "icc-in")




# DATA EXPLORATION [disabled for speed]
# summary(raw.cattle)
# View(raw.cattle)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information - coincides with NAs and 0 in state code
raw.cattle <-
  raw.cattle %>%
  dplyr::filter(!is.na(state_code) & state_code != 0)

# latin character treatment
raw.cattle <-
  raw.cattle %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.cattle$state_code)        <- "state code (2-digit, IBGE)"
sjlabelled::set_label(raw.cattle$state_name)        <- "state name"
sjlabelled::set_label(raw.cattle$date)             <- "trimester-year of reference"
sjlabelled::set_label(raw.cattle$cattleSlaughter_head) <- "number of cattle slaughtered"
sjlabelled::set_label(raw.cattle$cattleSlaughter_weight) <- "weight of cattle slaughtered"

# change object name for exportation
clean.surveyCattleSlaughter <- raw.cattle



# POST-TREATMENT OVERVIEW
# summary(clean.surveyCattleSlaughter)
# View(clean.surveyCattleSlaughter)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.surveyCattleSlaughter,
     file = here::here("data/raw2clean/surveyCattleSlaughter_ibge/output",
                      "clean_surveyCattleSlaughter.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------