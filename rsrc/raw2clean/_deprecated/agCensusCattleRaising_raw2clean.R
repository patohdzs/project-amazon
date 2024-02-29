
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW AGRICULTURAL USE AREA + WAGE EXPENDITURE + NUMBER OF WORKERS IN CATTLE RAISING FARMS - AGRICULTURAL CENSUS 2017 (IBGE)
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
here::i_am("code/raw2clean/agCensusCattleRaising_raw2clean.R", uuid = "ea330d66-0dad-4632-b640-5c1e96a5c47b")


# START TIME
tictoc::tic(msg = "agCensusCattleRaising_raw2clean script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read agricultural use area csv file
raw.agCensusAgUse <-  readr::read_csv(file = here::here("data/raw2clean/agCensusCattleRaising_ibge/input/agCensus_agUseArea_cattleRaising.csv"),
                                      skip = 5,
                                      na = c("..."),
                                      col_names = c("muni_code", "landUse_type", "landUse_establishments", "landUse_area"),
                                      col_types = list("n", "-", "-", "-", "c", "c", "c"))

# read wage expenditure csv file
raw.agCensusWage <-  readr::read_csv(file = here::here("data/raw2clean/agCensusCattleRaising_ibge/input/agCensus_wage_cattleRaising.csv"),
                                     skip = 6,
                                     na = c("..."),
                                     col_names = c("muni_code", "wage_expenditure"),
                                     col_types = list("n", "-", "-", "-", "-", "c"))


# read number of workers csv file
raw.agCensusWorkers <-  readr::read_csv(file = here::here("data/raw2clean/agCensusCattleRaising_ibge/input/agCensus_workers_cattleRaising.csv"),
                                     skip = 6,
                                     na = c("..."),
                                     col_names = c("muni_code", "workers_number"),
                                     col_types = list("n", "-", "-", "-", "c"))




# DATA EXPLORATION [disabled for speed]
# summary(raw.agCensusAgUse)
# View(raw.agCensusAgUse)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information
raw.agCensusAgUse <-  raw.agCensusAgUse[-(4633:4646), ]
raw.agCensusWage <-  raw.agCensusWage[-(773:786), ]
raw.agCensusWorkers <-  raw.agCensusWorkers[-(773:787), ]

# transform "-" values to "0" as explained in the table notes
raw.agCensusAgUse <-
  raw.agCensusAgUse %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))
raw.agCensusWage <-
  raw.agCensusWage %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))
raw.agCensusWorkers <-
  raw.agCensusWorkers %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# transform "X" values to "NA" (here NA identify which values had to be omitted to avoid informant identification as explained in the table notes)
raw.agCensusAgUse <-
  raw.agCensusAgUse %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))
raw.agCensusWage <-
  raw.agCensusWage %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))
raw.agCensusWorkers <-
  raw.agCensusWorkers %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))

# latin character treatment
raw.agCensusAgUse <-
  raw.agCensusAgUse %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))
raw.agCensusWage <-
  raw.agCensusWage %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))
raw.agCensusWorkers <-
  raw.agCensusWorkers %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))

# transform column class
raw.agCensusAgUse <-
  raw.agCensusAgUse %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("landUse_area"), function(x) as.numeric(x)))
raw.agCensusWage <-
  raw.agCensusWage %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("wage"), function(x) as.numeric(x)))
raw.agCensusWorkers <-
  raw.agCensusWorkers %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("workers"), function(x) as.numeric(x)))

# sum pasture and agricultural use area
raw.agCensusAgUse <-
  raw.agCensusAgUse %>%
  dplyr::mutate(d_pasture = if_else(landUse_type %in% c("Pastagens - naturais",
                                                        "Pastagens - pastagens plantadas em mas condicoes",
                                                        "Pastagens - plantadas em boas condicoes"),
                                    1, 0)) %>%
  dplyr::group_by(muni_code) %>%
  dplyr::summarise(pastureArea_value = sum(landUse_area*d_pasture, na.rm = TRUE),
                   agUseArea_value = sum(landUse_area, na.rm = TRUE))


# merge all variables
raw.agCensusCattleRaising <-
  raw.agCensusAgUse %>%
  dplyr::left_join(raw.agCensusWage) %>%
  dplyr::left_join(raw.agCensusWorkers)

# clean environment
rm(raw.agCensusAgUse, raw.agCensusWage, raw.agCensusWorkers)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.agCensusCattleRaising$muni_code)        <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.agCensusCattleRaising$pastureArea_value)    <- "pasture area (ha)"
sjlabelled::set_label(raw.agCensusCattleRaising$agUseArea_value)    <- "agricultural use area (ha)"
sjlabelled::set_label(raw.agCensusCattleRaising$wage_expenditure)    <- "total wage expenditure (thousand BRL)"
sjlabelled::set_label(raw.agCensusCattleRaising$workers_number)    <- "total number of workers (count)"

# change object name for exportation
clean.agCensusCattleRaising <- raw.agCensusCattleRaising



# POST-TREATMENT OVERVIEW
# summary(clean.agCensusCattleRaising)
# View(clean.agCensusCattleRaising)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.agCensusCattleRaising,
     file = here::here("data/raw2clean/agCensusCattleRaising_ibge/output",
                      "clean_agCensusCattleRaising.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------