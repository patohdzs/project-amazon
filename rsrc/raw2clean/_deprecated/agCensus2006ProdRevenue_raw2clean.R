
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW PRODUCTION REVENUES - AGRICULTURAL CENSUS 2006 (IBGE)
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
here::i_am("code/raw2clean/agCensus2006ProdRevenue_raw2clean.R", uuid = "eb3955a6-f041-436a-a011-19ff5d34d9aa")


# START TIME
tictoc::tic(msg = "agCensus2006ProdRevenue_raw2clean script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.agCensus2006ProdRevenue <-  readr::read_csv(file = here::here("data/raw2clean/agCensus2006ProdRevenue_ibge/input/agCensus2006_prodRevenue.csv"),
                                  skip = 6,
                                  na = c("..."),
                                  col_names = c("muni_code",
                                                "prodRevenueCrops_value_2006", "prodRevenueAnimal_value_2006", "prodRevenueAgIndustry_value_2006"),
                                  col_types = "n---ccc")




# DATA EXPLORATION [disabled for speed]
# summary(raw.agCensus2006ProdRevenue)
# View(raw.agCensus2006ProdRevenue)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information
raw.agCensus2006ProdRevenue <-  raw.agCensus2006ProdRevenue[-c(5548:5559),]

# transform "-" values to "0" as explained in the table notes
raw.agCensus2006ProdRevenue <-
  raw.agCensus2006ProdRevenue %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# transform "X" values to "NA" (here NA identify which values had to be omitted to avoid informant identification as explained in the table notes)
raw.agCensus2006ProdRevenue <-
  raw.agCensus2006ProdRevenue %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))

# latin character treatment
raw.agCensus2006ProdRevenue <-
  raw.agCensus2006ProdRevenue %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), iconv, from = "UTF-8", to = "ASCII//TRANSLIT"))


# transform column class
raw.agCensus2006ProdRevenue <-
  raw.agCensus2006ProdRevenue %>%
  dplyr::mutate(dplyr::across(tidyselect:::starts_with("prodRevenue"), function(x) as.numeric(x)))

# calculate total expenditures without including paid wages and contracted services that represent local benefits
raw.agCensus2006ProdRevenue <-
  raw.agCensus2006ProdRevenue %>%
  dplyr::mutate(prodRevenue_value_2006 = rowSums(across(c("prodRevenueCrops_value_2006", "prodRevenueAnimal_value_2006", "prodRevenueAgIndustry_value_2006")), na.rm = TRUE))





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.agCensus2006ProdRevenue$muni_code)        <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.agCensus2006ProdRevenue$prodRevenue_value_2006) <- "value of crop, animal, and agricultural industry revenues (thousand BRL, 2006 Ag Census)"

# change object name for exportation
clean.agCensus2006ProdRevenue <- raw.agCensus2006ProdRevenue



# POST-TREATMENT OVERVIEW
# summary(clean.agCensus2006ProdRevenue)
# View(clean.agCensus2006ProdRevenue)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.agCensus2006ProdRevenue,
     file = here::here("data/raw2clean/agCensus2006ProdRevenue_ibge/output",
                      "clean_agCensus2006ProdRevenue.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------