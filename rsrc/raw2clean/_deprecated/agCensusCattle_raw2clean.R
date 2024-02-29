
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW LIVESTOCK CATTLE AND TOTAL (HEAD COUNT) - AGRICULTURAL CENSUS 2017 (IBGE)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "agCensusCattle_raw2clean script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns, prefer than Hmisc::label because has function to clear labels when necessary




# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.agCensusCattle <-  readr::read_csv(file = "data/raw2clean/agCensusCattle_ibge/input/agCensus_Cattle.csv",
                                  skip = 6,
                                  na = c("..."),
                                  col_names = c("muni_code", "muni_name", "useless_column", "livestock_all", "livestock_cattle"),
                                  col_types = list("n", "c", "c", "c", "c"))




# DATA EXPLORATION [disabled for speed]
# summary(raw.agCensusCattle)
# View(raw.agCensusCattle)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information - coincides with NAs in useless_column column
raw.agCensusCattle <-
  raw.agCensusCattle %>%
  dplyr::filter(!is.na(useless_column))

# transform "-" values to "0" as explained in the table notes
raw.agCensusCattle <-
  raw.agCensusCattle %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# transform "X" values to "NA" (here NA identify which values had to be omitted to avoid informant identification as explained in the table notes)
raw.agCensusCattle <-
  raw.agCensusCattle %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))

# latin character treatment
raw.agCensusCattle <-
  raw.agCensusCattle %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))


# COLUMN CLEANUP

# transform column class
raw.agCensusCattle <-
  raw.agCensusCattle %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("livestock"), function(x) as.numeric(x)))

# remove useless columns
raw.agCensusCattle <-
  raw.agCensusCattle %>%
  dplyr::select(-starts_with("useless"))




# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.agCensusCattle$muni_code)        <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.agCensusCattle$muni_name)        <- "municipality name and state name abbreviation in parenthesis"
sjlabelled::set_label(raw.agCensusCattle$livestock_cattle) <- "number of cattle (head count)"
sjlabelled::set_label(raw.agCensusCattle$livestock_all)    <- "number of all livestock species (head count)"

# change object name for exportation
clean.agCensusCattle <- raw.agCensusCattle



# POST-TREATMENT OVERVIEW
# summary(clean.agCensusCattle)
# View(clean.agCensusCattle)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.agCensusCattle,
     file = file.path("data/raw2clean/agCensusCattle_ibge/output",
                      "clean_agCensusCattle.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------