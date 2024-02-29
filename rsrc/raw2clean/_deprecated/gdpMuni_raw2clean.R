# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN MUNINICIPAL GDP (ADDED VALUE) - AGRICULTURAL SECTOR
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# -





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
here::i_am("code/raw2clean/gdpMuni_raw2clean.R", uuid = "eb3955a6-f041-436a-a011-19ff5d34d9aa")


# START TIME
tictoc::tic(msg = "gdpMuni_raw2clean script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary






# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# RAW DATA INPUT
raw.gdp <- readr::read_csv(file      = here::here("data/raw2clean/gdpMuni_ibge/input/gdp_muni.csv"),
                           na = c("..."),
                           skip      = 3,
                           col_names = c("muni_code", "muni_name", "year", "valueAdded_agriculture"))



# DATA EXPLORATION [disabled for speed]
# class(raw.gdp)
# summary(raw.gdp)  # missing entries; has ',' as decimal separator (despite indication of decimal delimiter in 'read.csv2' fctn)
# View(raw.gdp)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# NON-DATA TREATMENT
raw.gdp <- raw.gdp[which(!is.na(raw.gdp$year)), ] # NA entry due to incoming file table notes



# COLUMN CLEANUP

# class
lapply(raw.gdp, class)


# remove unnecessary column
raw.gdp$muni_name <- NULL





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(raw.gdp$muni_code)              <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.gdp$year)                   <- "year of reference"
sjlabelled::set_label(raw.gdp$valueAdded_agriculture) <- "value added, agriculture (current prices, thousand BRL)"



# change objects name
clean.gdpMuni <- raw.gdp



# POST-TREATMENT OVERVIEW
# summary(clean.gdpMuni)
# View(clean.gdpMuni)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.gdpMuni,
     file = here::here("data/raw2clean/gdpMuni_ibge/output", "clean_gdpMuni.Rdata"))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")




# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------