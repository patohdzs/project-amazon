
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW ANIMAL PRODUCTION - AGRICULTURAL CENSUS 2017 (IBGE)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "agCensusAnimalProduction_raw2clean script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns, prefer than Hmisc::label because has function to clear labels when necessary




# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.agCensusAnimalProduction <-  readr::read_csv(file = "data/raw2clean/agCensusAnimalProduction_ibge/input/agCensus_animalProduction.csv",
                                  skip = 8,
                                  na = c("..."),
                                  col_names = c("muni_code", "muni_name", "animalProduction_establishments", "animalProduction_value"),
                                  col_types = list("n", "c", "c", "c"))




# DATA EXPLORATION [disabled for speed]
# summary(raw.agCensusAnimalProduction)
# View(raw.agCensusAnimalProduction)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information
raw.agCensusAnimalProduction <-  raw.agCensusAnimalProduction[-c(5563:5576),]

# transform "-" values to "0" as explained in the table notes
raw.agCensusAnimalProduction <-
  raw.agCensusAnimalProduction %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# transform "X" values to "NA" (here NA identify which values had to be omitted to avoid informant identification as explained in the table notes)
raw.agCensusAnimalProduction <-
  raw.agCensusAnimalProduction %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))

# latin character treatment
raw.agCensusAnimalProduction <-
  raw.agCensusAnimalProduction %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))


# transform column class
raw.agCensusAnimalProduction <-
  raw.agCensusAnimalProduction %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("animalProduction"), function(x) as.numeric(x)))





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.agCensusAnimalProduction$muni_code)        <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.agCensusAnimalProduction$muni_name)        <- "municipality name and state name abbreviation in parenthesis"
sjlabelled::set_label(raw.agCensusAnimalProduction$animalProduction_establishments) <- "number of establishments with large animal production (count)"
sjlabelled::set_label(raw.agCensusAnimalProduction$animalProduction_value)    <- "value of large animal production (thousand BRL)"

# change object name for exportation
clean.agCensusAnimalProduction <- raw.agCensusAnimalProduction



# POST-TREATMENT OVERVIEW
# summary(clean.agCensusAnimalProduction)
# View(clean.agCensusAnimalProduction)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.agCensusAnimalProduction,
     file = file.path("data/raw2clean/agCensusAnimalProduction_ibge/output",
                      "clean_agCensusAnimalProduction.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------