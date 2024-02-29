
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW PRODES DEFORESTATION RATE (INPE)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "prodesDeforestationRate_raw2clean script", log = T)

# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # sjlabelled::set_label columns, prefer than Hmisc::sjlabelled::set_label because has function to clear sjlabelled::set_labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read cs file
raw.prodesDef <- readr::read_csv2(file = "data/raw2clean/prodesDeforestationRate_inpe/input/terrabrasilis_legal_amazon_15_6_2021_1626370667166.csv")



# DATA EXPLORATION
#summary(raw.prodesDef)    # object is a list of data frames
#class(raw.prodesDef)
#View(raw.prodesDef)      # column names indicate file of origin





# DATASET CLEANUP AND PREP --------------------------------------------------------------------------------------------------------------------------

# transform to date class
raw.prodesDef <-
  raw.prodesDef %>%
  rename(prodes_year = year,
         prodes_deforestRate = `area km²`)




# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(raw.prodesDef$prodes_year)            <- "PRODES year (start = August of previous calendar year; end = July of current calendar year)"
sjlabelled::set_label(raw.prodesDef$prodes_deforestRate)    <- "PRODES deforestation rate in Legal Amazon (clear-cut above 6.25 ha)"


# change object name for exportation
clean.prodesDeforestationRate <- raw.prodesDef



# POST-TREATMENT OVERVIEW
# summary(clean.prodesDeforestationRate)
# View(clean.prodesDeforestationRate)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.prodesDeforestationRate,
     file = file.path("data/raw2clean/prodesDeforestationRate_inpe/output",
                      "clean_prodesDeforestationRate.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT: ------------------------------------------------------------------------------------------------------------------------------------