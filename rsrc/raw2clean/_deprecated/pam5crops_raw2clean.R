
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW SELECTED CROPS DATA - PAM (IBGE)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "pam5crops_raw2clean script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns, prefer than Hmisc::label because has function to clear labels when necessary




# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.pam5crops <-  readr::read_csv(file = "data/raw2clean/pam5crops_ibge/input/pam_5crops.csv",
                                  skip = 4,
                                  na = c("..."),
                                  col_names = c("muni_code", "muni_name", "year",
                                                "planted_area_rice", "harvested_area_rice",
                                                "quantity_produced_rice", "value_production_rice",
                                                "planted_area_sugarcane", "harvested_area_sugarcane",
                                                "quantity_produced_sugarcane", "value_production_sugarcane",
                                                "planted_area_cassava", "harvested_area_cassava",
                                                "quantity_produced_cassava", "value_production_cassava",
                                                "planted_area_corn", "harvested_area_corn",
                                                "quantity_produced_corn", "value_production_corn",
                                                "planted_area_soybean", "harvested_area_soybean",
                                                "quantity_produced_soybean", "value_production_soybean"),
                                  col_types = list("n", "c", "i", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c",
                                                   "c", "c", "c"))




# DATA EXPLORATION [disabled for speed]
# summary(raw.pam5crops)
# View(raw.pam5crops)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information - coincides with NAs in year column
raw.pam5crops <-
  raw.pam5crops %>%
  dplyr::filter(!is.na(year))

# transform "-" values to "0" as explained in the "ibge_technicalNotes_pam2019.pdf" file in documentation
raw.pam5crops <-
  raw.pam5crops %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# latin character treatment
raw.pam5crops <-
  raw.pam5crops %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))


# COLUMN CLEANUP

# transform column class
raw.pam5crops <-
  raw.pam5crops %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("planted"), function(x) as.numeric(x))) %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("harvested"), function(x) as.numeric(x))) %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("quantity"), function(x) as.numeric(x))) %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("value"), function(x) as.numeric(x)))




# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.pam5crops$muni_code)                   <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.pam5crops$muni_name)                   <- "municipality name and (state name abbreviation)"
sjlabelled::set_label(raw.pam5crops$year)                        <- "year of reference (calendar or PRODES year)"

sjlabelled::set_label(raw.pam5crops$planted_area_rice)           <- "planted area - rice (ha)"
sjlabelled::set_label(raw.pam5crops$harvested_area_rice)         <- "harvested area - rice (ha)"
sjlabelled::set_label(raw.pam5crops$quantity_produced_rice)      <- "quantity produced - rice (t)"
sjlabelled::set_label(raw.pam5crops$value_production_rice)       <- "monetary value of production - rice (BRL thousand)"


sjlabelled::set_label(raw.pam5crops$planted_area_sugarcane)      <- "planted area - sugarcane (ha)"
sjlabelled::set_label(raw.pam5crops$harvested_area_sugarcane)    <- "harvested area - sugarcane (ha)"
sjlabelled::set_label(raw.pam5crops$quantity_produced_sugarcane) <- "quantity produced - sugarcane (t)"
sjlabelled::set_label(raw.pam5crops$value_production_sugarcane)  <- "monetary value of production - sugarcane (BRL thousand)"

sjlabelled::set_label(raw.pam5crops$planted_area_cassava)        <- "planted area - cassava (ha)"
sjlabelled::set_label(raw.pam5crops$harvested_area_cassava)      <- "harvested area - cassava (ha)"
sjlabelled::set_label(raw.pam5crops$quantity_produced_cassava)   <- "quantity produced - cassava (t)"
sjlabelled::set_label(raw.pam5crops$value_production_cassava)    <- "monetary value of production - cassava (BRL thousand)"

sjlabelled::set_label(raw.pam5crops$planted_area_corn)           <- "planted area - corn (ha)"
sjlabelled::set_label(raw.pam5crops$harvested_area_corn)         <- "harvested area - corn (ha)"
sjlabelled::set_label(raw.pam5crops$quantity_produced_corn)      <- "quantity produced - corn (t)"
sjlabelled::set_label(raw.pam5crops$value_production_corn)       <- "monetary value of production - corn (BRL thousand)"

sjlabelled::set_label(raw.pam5crops$planted_area_soybean)        <- "planted area - soybean (ha)"
sjlabelled::set_label(raw.pam5crops$harvested_area_soybean)      <- "harvested area - soybean (ha)"
sjlabelled::set_label(raw.pam5crops$quantity_produced_soybean)   <- "quantity produced - soybean (t)"
sjlabelled::set_label(raw.pam5crops$value_production_soybean)    <- "monetary value of production - soybean (BRL thousand)"

# change object name for exportation
clean.pam5crops <- raw.pam5crops



# POST-TREATMENT OVERVIEW
# summary(clean.pam5crops)
# View(clean.pam5crops)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.pam5crops,
     file = file.path("data/raw2clean/pam5crops_ibge/output",
                      "clean_pam5crops.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------