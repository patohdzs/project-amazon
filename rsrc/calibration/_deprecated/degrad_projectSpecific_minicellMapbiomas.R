
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: ADD DEGRADATION DATA (DEGRAD/INPE - 2007-2016)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "degrad_projectSpecific_minicellMapbiomas script", log = T)



# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(sf) # manipulate spatial data
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns, prefer than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# MAPBIOMAS SAMPLE
load("data/projectSpecific/minicellMapbiomas/mapbiomasCategories_minicellMapbiomas_2000.Rdata")



# VECTOR DATA - DEGRAD
load("data/raw2clean/degrad_inpe/output/clean_degrad.Rdata")






# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# ADD DEGRAD DATA TO MAPBIOMAS SAMPLE

# transform sample to spatial data
mapbiomasCategories.minicellMapbiomas.2000 <- sf::st_as_sf(x = mapbiomasCategories.minicellMapbiomas.2000,
                                                           coords = c("lon", "lat"),
                                                           crs = sf::st_crs("+proj=longlat +datum=WGS84 +no_defs"))

# adjust crs
clean.degrad <- sf::st_transform(clean.degrad, sf::st_crs(mapbiomasCategories.minicellMapbiomas.2000))

# add vector of true or false for each year, TRUE := the point was degradated in that year FALSE := otherwise
mapbiomasCategories.minicellMapbiomas.2000$d_degrad_2007 <- lengths(sf::st_intersects(mapbiomasCategories.minicellMapbiomas.2000,
                                                                                      clean.degrad[clean.degrad$year == 2007,])) > 0
mapbiomasCategories.minicellMapbiomas.2000$d_degrad_2008 <- lengths(sf::st_intersects(mapbiomasCategories.minicellMapbiomas.2000,
                                                                                      clean.degrad[clean.degrad$year == 2008,])) > 0
mapbiomasCategories.minicellMapbiomas.2000$d_degrad_2009 <- lengths(sf::st_intersects(mapbiomasCategories.minicellMapbiomas.2000,
                                                                                      clean.degrad[clean.degrad$year == 2009,])) > 0
mapbiomasCategories.minicellMapbiomas.2000$d_degrad_2010 <- lengths(sf::st_intersects(mapbiomasCategories.minicellMapbiomas.2000,
                                                                                      clean.degrad[clean.degrad$year == 2010,])) > 0
mapbiomasCategories.minicellMapbiomas.2000$d_degrad_2011 <- lengths(sf::st_intersects(mapbiomasCategories.minicellMapbiomas.2000,
                                                                                      clean.degrad[clean.degrad$year == 2011,])) > 0
mapbiomasCategories.minicellMapbiomas.2000$d_degrad_2012 <- lengths(sf::st_intersects(mapbiomasCategories.minicellMapbiomas.2000,
                                                                                      clean.degrad[clean.degrad$year == 2012,])) > 0
mapbiomasCategories.minicellMapbiomas.2000$d_degrad_2013 <- lengths(sf::st_intersects(mapbiomasCategories.minicellMapbiomas.2000,
                                                                                      clean.degrad[clean.degrad$year == 2013,])) > 0
mapbiomasCategories.minicellMapbiomas.2000$d_degrad_2014 <- lengths(sf::st_intersects(mapbiomasCategories.minicellMapbiomas.2000,
                                                                                      clean.degrad[clean.degrad$year == 2014,])) > 0
mapbiomasCategories.minicellMapbiomas.2000$d_degrad_2015 <- lengths(sf::st_intersects(mapbiomasCategories.minicellMapbiomas.2000,
                                                                                      clean.degrad[clean.degrad$year == 2015,])) > 0
mapbiomasCategories.minicellMapbiomas.2000$d_degrad_2016 <- lengths(sf::st_intersects(mapbiomasCategories.minicellMapbiomas.2000,
                                                                                      clean.degrad[clean.degrad$year == 2016,])) > 0

# drop spatial feature
degrad.projectSpecific.minicellMapbiomas <- cbind(sf::st_drop_geometry(mapbiomasCategories.minicellMapbiomas.2000),
                                  sf::st_coordinates(mapbiomasCategories.minicellMapbiomas.2000))

# clear environment
rm(clean.degrad, mapbiomasCategories.minicellMapbiomas.2000)

# transform data into panel
degrad.projectSpecific.minicellMapbiomas <-
  degrad.projectSpecific.minicellMapbiomas %>%
  dplyr::select(lon = X, lat = Y, starts_with("d_degrad_")) %>%
  tidyr::pivot_longer(cols = -c("lon", "lat"), names_to = "year", values_to = "d_degrad") %>%
  dplyr::mutate(year = as.numeric(stringr::str_extract(year, "\\d{4}")),
                d_degrad = dplyr::if_else(d_degrad == T, 1, 0))





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(degrad.projectSpecific.minicellMapbiomas$d_degrad) <- "=1 if pixel was degradated in that year, =0 otherwise (PRODES year := from Aug/t-1 through Jul/t)"



# POST-TREATMENT OVERVIEW
# summary(degrad.projectSpecific.minicellMapbiomas)
# View(degrad.projectSpecific.minicellMapbiomas)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(degrad.projectSpecific.minicellMapbiomas,
     file = file.path("data/projectSpecific/minicellMapbiomas",
                      paste0("degrad_minicellMapbiomas", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/minicellMapbiomas")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
