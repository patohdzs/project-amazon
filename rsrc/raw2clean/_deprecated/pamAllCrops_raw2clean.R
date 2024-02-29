
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW TOTAL AllORARY CROPS DATA - PAM (IBGE)
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
here::i_am("code/raw2clean/pamAllCrops_raw2clean.R", uuid = "eb3955a6-f041-436a-a011-19ff5d34d9aa")


# START TIME
tictoc::tic(msg = "pamAllCrops_raw2clean script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary




# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.pamAllCrops <-  readr::read_csv(file = here::here("data/raw2clean/pamAllCrops_ibge/input/pam_AllCrops.csv"),
                                  skip = 4,
                                  na = c("..."),
                                  col_names = c("muni_code", "muni_name", "year",
                                                "planted_area_AllCrops", "harvested_area_AllCrops",
                                                "quantity_produced_AllCrops", "value_production_AllCrops"),
                                  col_types = list("n", "c", "i", "c", "c", "c", "c"))




# DATA EXPLORATION [disabled for speed]
# summary(raw.pamAllCrops)
# View(raw.pamAllCrops)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information - coincides with NAs in year column
raw.pamAllCrops <-
  raw.pamAllCrops %>%
  dplyr::filter(!is.na(year))

# transform "-" values to "0" as explained in the "ibge_technicalNotes_pam2019.pdf" file in documentation
raw.pamAllCrops <-
  raw.pamAllCrops %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# latin character treatment
raw.pamAllCrops <-
  raw.pamAllCrops %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))


# COLUMN CLEANUP

# remove useless column
raw.pamAllCrops <-
  raw.pamAllCrops %>%
  dplyr::select(-quantity_produced_AllCrops)

# transform column class
raw.pamAllCrops <-
  raw.pamAllCrops %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("planted"), function(x) as.numeric(x))) %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("harvested"), function(x) as.numeric(x))) %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("value"), function(x) as.numeric(x)))





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.pamAllCrops$muni_code)                   <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.pamAllCrops$muni_name)                   <- "municipality name and (state name abbreviation)"
sjlabelled::set_label(raw.pamAllCrops$year)                        <- "year of reference (calendar year)"

sjlabelled::set_label(raw.pamAllCrops$planted_area_AllCrops)           <- "planted area - total Allorary crops (ha)"
sjlabelled::set_label(raw.pamAllCrops$harvested_area_AllCrops)         <- "harvested area - total Allorary crops (ha)"
sjlabelled::set_label(raw.pamAllCrops$value_production_AllCrops)       <- "monetary value of production - total Allorary cops (BRL thousand)"

# change object name for exportation
clean.pamAllCrops <- raw.pamAllCrops



# POST-TREATMENT OVERVIEW
# summary(clean.pamAllCrops)
# View(clean.pamAllCrops)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.pamAllCrops,
     file = here::here("data/raw2clean/pamAllCrops_ibge/output",
                      "clean_pamAllCrops.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------