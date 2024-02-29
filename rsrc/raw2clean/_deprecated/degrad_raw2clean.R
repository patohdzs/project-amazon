
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CLEAN RAW DEGRADATION DATA (DEGRAD - INPE)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "degrad_raw2clean script", log = T)


# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(sf) # manipulate spatial data
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns, prefer than Hmisc::label because has function to clear labels when necessary




# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RAW DATA

# create auxiliary vector of input data years
aux.years  <- 2007:2016

# read the shapefile of each year
raw.degrad <- lapply(X   = aux.years,
                     FUN = function(x) sf::st_read(paste0("data/raw2clean/degrad_inpe/input",
                                                           "/", as.character(x))))



# EXPLORATION
# summary(raw.degrad)
# lapply(raw.degrad, summary) # columns are not consistent throughout sample
# lapply(raw.degrad, sf::st_crs)  # all missing due to reading error, 2007-2009 lack .prj file but checked with source that it has the same EPSG:4618 as 2010-2014, while 2015-2016 have EPSG:4674
# plot(raw.degrad[[1]]$geometry)





# DATASET PREP AND CLEANUP ---------------------------------------------------------------------------------------------------------------------------

# LIST ELEMENT NAMES
names(raw.degrad) <- as.list(aux.years)



# CRS ATTRIBUTION - see _metadata.txt
for (i in c(2007:2014)) {
  sf::st_crs(raw.degrad[[as.character(i)]]) <- 4618
}
for (i in c(2015:2016)) {
  sf::st_crs(raw.degrad[[as.character(i)]]) <- 4674
}



# COLUMN CLEANUP

# see all column names
lapply(X   = raw.degrad,
       FUN = names)

# standardizes lower/upper case across column names
raw.degrad$`2010` <- raw.degrad$`2010` %>% dplyr::rename_all(.funs = tolower)

# change uf to state_uf and convert it to character
for (i in c(2007:2009, 2011:2015)) {

  raw.degrad[[as.character(i)]]          <- raw.degrad[[as.character(i)]] %>% dplyr::rename(state_uf = uf)
  raw.degrad[[as.character(i)]]$state_uf <- as.character(raw.degrad[[as.character(i)]]$state_uf)


}

# translate and standardize column names
raw.degrad$`2009` <- raw.degrad$`2009` %>% dplyr::rename(year = ano)
raw.degrad$`2010` <- raw.degrad$`2010` %>% dplyr::rename(state_name = nome)
raw.degrad$`2010` <- raw.degrad$`2010` %>% dplyr::rename(state_code = codigouf)
raw.degrad$`2012` <- raw.degrad$`2012` %>% dplyr::rename(area_meters = areametros)
raw.degrad$`2013` <- raw.degrad$`2013` %>% dplyr::rename(area_meters = areameters)
raw.degrad$`2014` <- raw.degrad$`2014` %>% dplyr::rename(area_meters = areameters)
raw.degrad$`2015` <- raw.degrad$`2015` %>% dplyr::rename(area_meters = areameters)
raw.degrad$`2016` <- raw.degrad$`2016` %>% dplyr::rename(area_meters = areameters)

for (i in aux.years) {

  raw.degrad[[as.character(i)]] <- raw.degrad[[as.character(i)]] %>% dplyr::rename(degrad_type = class_name)

}


# change column class
for (i in aux.years) {

  raw.degrad[[as.character(i)]]$linkcolumn <- as.numeric(as.character(raw.degrad[[as.character(i)]]$linkcolumn))
  raw.degrad[[as.character(i)]]$pathrow <- as.numeric(as.character(raw.degrad[[as.character(i)]]$pathrow))

}

for (i in c(2010:2016)) {

  raw.degrad[[as.character(i)]]$view_date <- lubridate::ymd(as.character(raw.degrad[[as.character(i)]]$view_date)) # change to date class

}

raw.degrad$`2010`$state_name <- as.character(raw.degrad$`2010`$state_name)
raw.degrad$`2010`$state_code <- as.numeric(as.character(raw.degrad$`2010`$state_code))


# add year indicator
for (i in aux.years) {
  raw.degrad[[as.character(i)]]$degrad_year <- as.numeric(i)
}


# add missing columns to allow rbinding the data

# list all column names
columns <- unique(unlist(lapply(X = raw.degrad, FUN = names)))

# if column is not present add it with all values as NA
for (i in aux.years) {
  for (j in columns) {

  if(!any(colnames(raw.degrad[[as.character(i)]]) %in% as.character(j))) {

    raw.degrad[[as.character(i)]][as.character(j)] <- NA

  } else {
      next
    }
  }
}



# ROW CLEANUP
# standardize degrad_type elements
raw.degrad$`2007` <- raw.degrad$`2007` %>% dplyr::mutate(degrad_type = dplyr::recode(degrad_type, DEGRAD2007 = "degradation"))
raw.degrad$`2008` <- raw.degrad$`2008` %>% dplyr::mutate(degrad_type = dplyr::recode(degrad_type, DEGRAD2008 = "degradation"))
raw.degrad$`2009` <- raw.degrad$`2009` %>% dplyr::mutate(degrad_type = dplyr::recode(degrad_type, DEGRADACAO = "degradation"))
raw.degrad$`2010` <- raw.degrad$`2010` %>% dplyr::mutate(degrad_type = dplyr::recode(degrad_type, DEGRAD = "degradation"))
raw.degrad$`2011` <- raw.degrad$`2011` %>% dplyr::mutate(degrad_type = dplyr::recode(degrad_type, DEGRAD = "degradation"))
raw.degrad$`2012` <- raw.degrad$`2012` %>% dplyr::mutate(degrad_type = dplyr::recode(degrad_type, DEGRAD = "degradation"))
raw.degrad$`2013` <- raw.degrad$`2013` %>% dplyr::mutate(degrad_type = dplyr::recode(degrad_type, DEGRAD = "degradation"))
raw.degrad$`2014` <- raw.degrad$`2014` %>% dplyr::mutate(degrad_type = dplyr::recode(degrad_type, DEGRAD = "degradation"))
raw.degrad$`2015` <- raw.degrad$`2015` %>% dplyr::mutate(degrad_type = dplyr::recode(degrad_type, DEGRAD = "degradation",
                                                                                            BLOWDOWN = "blowdown",
                                                                                            CICATRIZ_QUEIMADA = "fire_scar",
                                                                                            DEGRAD_NATURAL = "degradation_natural",
                                                                                            QUEIMADA = "fire"))
raw.degrad$`2016` <- raw.degrad$`2016` %>% dplyr::mutate(degrad_type = dplyr::recode(degrad_type, DEGRAD = "degradation",
                                                                                            CICATRIZ_INCENDIO = "fire_scar"))



# PROJECTION
raw.degrad <- lapply(raw.degrad, sf::st_transform, crs = 5880) # SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)



# GEOMETRY CLEANUP
raw.degrad <- lapply(raw.degrad, sf::st_make_valid) # needed to skip this step because it was giving error afterwards when trying to convert from sf to sp



# MERGE ALL YEARS
raw.degrad <- do.call(rbind, raw.degrad)


# KEEP COLUMNS OF INTEREST
raw.degrad <-
  raw.degrad %>%
  dplyr::select(degrad_type, state_uf, area, -year, year = degrad_year, area_meters)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# change object name for exportation
clean.degrad <- raw.degrad

# clear enviroment
rm(raw.degrad)


# POST-TREATMENT OVERVIEW
# summary(clean.degrad)
# View(clean.degrad)
# plot(clean.degrad)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.degrad,
     file = file.path("data/raw2clean/degrad_inpe/output",
                      "clean_degrad.Rdata"))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
