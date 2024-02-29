
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW WAGE EXPENSE - AGRICULTURAL CENSUS 2017 (IBGE)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "agCensusWage_raw2clean script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns, prefer than Hmisc::label because has function to clear labels when necessary




# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.agCensusWage <-  readr::read_csv(file = "data/raw2clean/agCensusWage_ibge/input/agCensus_wage.csv",
                                  skip = 6,
                                  na = c("..."),
                                  col_names = c("muni_code", "wage_expense"),
                                  col_types = list("n", "-", "-", "c"))




# DATA EXPLORATION [disabled for speed]
# summary(raw.agCensusWage)
# View(raw.agCensusWage)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information
raw.agCensusWage <-  raw.agCensusWage[!is.na(raw.agCensusWage$wage_expense),]

# transform "-" values to "0" as explained in the table notes
raw.agCensusWage <-
  raw.agCensusWage %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# transform "X" values to "NA" (here NA identify which values had to be omitted to avoid informant identification as explained in the table notes)
raw.agCensusWage <-
  raw.agCensusWage %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))

# transform column class
raw.agCensusWage <-
  raw.agCensusWage %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("wage_expense"), function(x) as.numeric(x)))




# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.agCensusWage$muni_code)        <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.agCensusWage$wage_expense)    <- "wage expense in 2017 (thousand BRL)"

# change object name for exportation
clean.agCensusWage <- raw.agCensusWage



# POST-TREATMENT OVERVIEW
# summary(clean.agCensusWage)
# View(clean.agCensusWage)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.agCensusWage,
     file = file.path("data/raw2clean/agCensusWage_ibge/output",
                      "clean_agCensusWage.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------