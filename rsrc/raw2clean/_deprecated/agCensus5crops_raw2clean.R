
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW SELECTED CROPS DATA - AGRICULTURAL CENSUS 2017 (IBGE)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "agCensus5crops_raw2clean script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns, prefer than Hmisc::label because has function to clear labels when necessary




# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.agCensus5crops <-  readr::read_csv(file = "data/raw2clean/agCensus5crops_ibge/input/agCensus_5crops.csv",
                                  skip = 6,
                                  na = c("..."),
                                  col_names = c("muni_code", "muni_name", "useless_column", "useless_column2",
                                                "harvested_area_tempCrops", "harvested_area_rice", "harvested_area_sugarcane", "harvested_area_cassava",
                                                "harvested_area_corn", "harvested_area_soybean"),
                                  col_types = list("n", "c", "c", "c", "c", "c", "c", "c", "c", "c"))




# DATA EXPLORATION [disabled for speed]
# summary(raw.agCensus5crops)
# View(raw.agCensus5crops)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information - coincides with NAs in useless_column2 column
raw.agCensus5crops <-
  raw.agCensus5crops %>%
  dplyr::filter(!is.na(useless_column2))

# transform "-" values to "0" as explained in the table notes
raw.agCensus5crops <-
  raw.agCensus5crops %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# transform "X" values to "NA" (here NA identify which values had to be omitted to avoid informant identification as explained in the table notes)
raw.agCensus5crops <-
  raw.agCensus5crops %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))

# latin character treatment
raw.agCensus5crops <-
  raw.agCensus5crops %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))


# COLUMN CLEANUP

# transform column class
raw.agCensus5crops <-
  raw.agCensus5crops %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("harvested"), function(x) as.numeric(x)))

# remove useless columns
raw.agCensus5crops <-
  raw.agCensus5crops %>%
  dplyr::select(-starts_with("useless"))




# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.agCensus5crops$muni_code)                   <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.agCensus5crops$muni_name)                   <- "municipality name and state name abbreviation in parenthesis"
sjlabelled::set_label(raw.agCensus5crops$harvested_area_rice)         <- "harvested area - rice (ha)"
sjlabelled::set_label(raw.agCensus5crops$harvested_area_sugarcane)    <- "harvested area - sugarcane (ha)"
sjlabelled::set_label(raw.agCensus5crops$harvested_area_cassava)      <- "harvested area - cassava (ha)"
sjlabelled::set_label(raw.agCensus5crops$harvested_area_corn)         <- "harvested area - corn (ha)"
sjlabelled::set_label(raw.agCensus5crops$harvested_area_soybean)      <- "harvested area - soybean (ha)"
sjlabelled::set_label(raw.agCensus5crops$harvested_area_tempCrops)    <- "harvested area - total of all temporary crops (ha)"

# change object name for exportation
clean.agCensus5crops <- raw.agCensus5crops



# POST-TREATMENT OVERVIEW
# summary(clean.agCensus5crops)
# View(clean.agCensus5crops)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.agCensus5crops,
     file = file.path("data/raw2clean/agCensus5crops_ibge/output",
                      "clean_agCensus5crops.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------