
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW TOTAL TEMPORARY CROPS DATA - PAM (IBGE)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "pamTempCrops_raw2clean script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns, prefer than Hmisc::label because has function to clear labels when necessary




# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.pamTempCrops <-  readr::read_csv(file = "data/raw2clean/pamTempCrops_ibge/input/pam_tempCrops.csv",
                                  skip = 4,
                                  na = c("..."),
                                  col_names = c("muni_code", "muni_name", "year",
                                                "planted_area_tempCrops", "harvested_area_tempCrops",
                                                "quantity_produced_tempCrops", "value_production_tempCrops"),
                                  col_types = list("n", "c", "i", "c", "c", "c", "c"))




# DATA EXPLORATION [disabled for speed]
# summary(raw.pamTempCrops)
# View(raw.pamTempCrops)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information - coincides with NAs in year column
raw.pamTempCrops <-
  raw.pamTempCrops %>%
  dplyr::filter(!is.na(year))

# transform "-" values to "0" as explained in the "ibge_technicalNotes_pam2019.pdf" file in documentation
raw.pamTempCrops <-
  raw.pamTempCrops %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# latin character treatment
raw.pamTempCrops <-
  raw.pamTempCrops %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))


# COLUMN CLEANUP

# remove useless column
raw.pamTempCrops <-
  raw.pamTempCrops %>%
  dplyr::select(-quantity_produced_tempCrops)

# transform column class
raw.pamTempCrops <-
  raw.pamTempCrops %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("planted"), function(x) as.numeric(x))) %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("harvested"), function(x) as.numeric(x))) %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("value"), function(x) as.numeric(x)))





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.pamTempCrops$muni_code)                   <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.pamTempCrops$muni_name)                   <- "municipality name and (state name abbreviation)"
sjlabelled::set_label(raw.pamTempCrops$year)                        <- "year of reference (calendar year)"

sjlabelled::set_label(raw.pamTempCrops$planted_area_tempCrops)           <- "planted area - total temporary crops (ha)"
sjlabelled::set_label(raw.pamTempCrops$harvested_area_tempCrops)         <- "harvested area - total temporary crops (ha)"
sjlabelled::set_label(raw.pamTempCrops$value_production_tempCrops)       <- "monetary value of production - total temporary cops (BRL thousand)"

# change object name for exportation
clean.pamTempCrops <- raw.pamTempCrops



# POST-TREATMENT OVERVIEW
# summary(clean.pamTempCrops)
# View(clean.pamTempCrops)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.pamTempCrops,
     file = file.path("data/raw2clean/pamTempCrops_ibge/output",
                      "clean_pamTempCrops.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------