
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW PERMANENT CROPS SOLD - AGRICULTURAL CENSUS 2017 (IBGE)
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
here::i_am("code/raw2clean/agCensus2017PermCropsSold_raw2clean.R", uuid = "eb3955a6-f041-436a-a011-19ff5d34d9aa")


# START TIME
tictoc::tic(msg = "agCensus2017PermCropsSold_raw2clean script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.agCensus2017PermCropsSold <-  readr::read_csv(file = here::here("data/raw2clean/agCensus2017PermCropsSold_ibge/input/agCensus2017_permCropsSold.csv"),
                                  skip = 6,
                                  na = c("..."),
                                  col_names = c("muni_code", "permCropsLargeProp_value_2017", "permCropsSmallProp_value_2017"),
                                  col_types = "n--cc")




# DATA EXPLORATION [disabled for speed]
# summary(raw.agCensus2017PermCropsSold)
# View(raw.agCensus2017PermCropsSold)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information
raw.agCensus2017PermCropsSold <-  raw.agCensus2017PermCropsSold[-c(5490:5504),]

# transform "-" values to "0" as explained in the table notes
raw.agCensus2017PermCropsSold <-
  raw.agCensus2017PermCropsSold %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# transform "X" values to "NA" (here NA identify which values had to be omitted to avoid informant identification as explained in the table notes)
raw.agCensus2017PermCropsSold <-
  raw.agCensus2017PermCropsSold %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))

# latin character treatment
raw.agCensus2017PermCropsSold <-
  raw.agCensus2017PermCropsSold %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))


# transform column class
raw.agCensus2017PermCropsSold <-
  raw.agCensus2017PermCropsSold %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("permCrops"), function(x) as.numeric(x)))

# calculate total value of permanent crops sold
raw.agCensus2017PermCropsSold <-
  raw.agCensus2017PermCropsSold %>%
  dplyr::mutate(permCrops_value_2017 = rowSums(across(c("permCropsSmallProp_value_2017", "permCropsLargeProp_value_2017")), na.rm = TRUE))


# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.agCensus2017PermCropsSold$muni_code) <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.agCensus2017PermCropsSold$permCropsSmallProp_value_2017) <- "value of permanent crops sold in properties with less than 50 trees (thousand BRL, 2017 Ag Census)"
sjlabelled::set_label(raw.agCensus2017PermCropsSold$permCropsLargeProp_value_2017) <- "value of permanent crops sold in properties with more than 50 trees (thousand BRL, 2017 Ag Census)"
sjlabelled::set_label(raw.agCensus2017PermCropsSold$permCrops_value_2017) <- "value of permanent crops sold (thousand BRL, 2017 Ag Census)"

# change object name for exportation
clean.agCensus2017PermCropsSold <- raw.agCensus2017PermCropsSold



# POST-TREATMENT OVERVIEW
# summary(clean.agCensus2017PermCropsSold)
# View(clean.agCensus2017PermCropsSold)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.agCensus2017PermCropsSold,
     file = here::here("data/raw2clean/agCensus2017PermCropsSold_ibge/output",
                      "clean_agCensus2017PermCropsSold.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------