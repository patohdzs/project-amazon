
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW PRODUCTION EXPENDITURES - AGRICULTURAL CENSUS 2006 (IBGE)
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
here::i_am("code/raw2clean/agCensus2006ProdExpenditure_raw2clean.R", uuid = "eb3955a6-f041-436a-a011-19ff5d34d9aa")


# START TIME
tictoc::tic(msg = "agCensus2006ProdExpenditure_raw2clean script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.agCensus2006ProdExpenditure <-  readr::read_csv(file = here::here("data/raw2clean/agCensus2006ProdExpenditure_ibge/input/agCensus2006_prodExpenditure.csv"),
                                  skip = 6,
                                  na = c("..."),
                                  col_names = c("muni_code",
                                                "prodExpenditureTotal_value_2006", "prodExpenditureService_value_2006", "prodExpenditureWageFamily_value_2006", "prodExpenditureWageWorker_value_2006"),
                                  col_types = "n---cccc")




# DATA EXPLORATION [disabled for speed]
# summary(raw.agCensus2006ProdExpenditure)
# View(raw.agCensus2006ProdExpenditure)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information
raw.agCensus2006ProdExpenditure <-  raw.agCensus2006ProdExpenditure[-c(5548:5559),]

# transform "-" values to "0" as explained in the table notes
raw.agCensus2006ProdExpenditure <-
  raw.agCensus2006ProdExpenditure %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# transform "X" values to "NA" (here NA identify which values had to be omitted to avoid informant identification as explained in the table notes)
raw.agCensus2006ProdExpenditure <-
  raw.agCensus2006ProdExpenditure %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))

# latin character treatment
raw.agCensus2006ProdExpenditure <-
  raw.agCensus2006ProdExpenditure %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))


# transform column class
raw.agCensus2006ProdExpenditure <-
  raw.agCensus2006ProdExpenditure %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("prodExpenditure"), function(x) as.numeric(x)))

# calculate total expenditures without including paid wages and contracted services that represent local benefits
raw.agCensus2006ProdExpenditure <-
  raw.agCensus2006ProdExpenditure %>%
  dplyr::mutate(prodExpenditureAux_value_2006 = rowSums(across(c("prodExpenditureService_value_2006", "prodExpenditureWageFamily_value_2006", "prodExpenditureWageWorker_value_2006")), na.rm = TRUE),
                prodExpenditureAux_value_2006 = -prodExpenditureAux_value_2006, # flip sign to subtract from total in rowSums below
                prodExpenditure_value_2006 = rowSums(across(c("prodExpenditureTotal_value_2006", "prodExpenditureAux_value_2006")), na.rm = TRUE)) %>%
  dplyr::select(-prodExpenditureAux_value_2006)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.agCensus2006ProdExpenditure$muni_code)        <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.agCensus2006ProdExpenditure$prodExpenditure_value_2006) <- "value of expenditures subtracting paid wages and contracted services (thousand BRL, 2006 Ag Census)"

# change object name for exportation
clean.agCensus2006ProdExpenditure <- raw.agCensus2006ProdExpenditure



# POST-TREATMENT OVERVIEW
# summary(clean.agCensus2006ProdExpenditure)
# View(clean.agCensus2006ProdExpenditure)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.agCensus2006ProdExpenditure,
     file = here::here("data/raw2clean/agCensus2006ProdExpenditure_ibge/output",
                      "clean_agCensus2006ProdExpenditure.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------