
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
tictoc::tic(msg = "degrad_projectSpecific_muniLevel script", log = T)



# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(sf) # manipulate spatial data
library(tidyverse) # manipulate tables, works with sf
library(sjlabelled) # label columns, prefer than Hmisc::label because has function to clear labels when necessary





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# MAPBIOMAS SAMPLE
load("data/projectSpecific/muniLevel/sampleSpatial_muniLevel.Rdata")

load("data/projectSpecific/muniLevel/samplePanel_muniLevel.Rdata")


# VECTOR DATA - DEGRAD
load("data/raw2clean/degrad_inpe/output/clean_degrad.Rdata")






# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# ADD DEGRAD DATA TO MAPBIOMAS SAMPLE

# select columns of interest
clean.degrad <-
  clean.degrad %>%
  dplyr::select(year, degrad_area)


# check crs
if (sf::st_crs(sampleSpatial.muniLevel) == sf::st_crs(clean.degrad)) {

  # combine sample municipalities with degradation polygons
  aux.degrad <- sf::st_intersection(sampleSpatial.muniLevel, clean.degrad)

} else {

  print("Spatial objects are on different projections, need to adjust before intersection!!")

}

# calculate degradated areas inside each municipality by year
aux.degrad$degrad_area <-
  sf::st_area(aux.degrad) %>%
  units::set_units(ha) %>%
  unclass()

# aggregate total degradated area by municipality and year
aux.degrad <-
  aux.degrad %>%
  sf::st_drop_geometry() %>% # drop spatial feature
  dplyr::group_by(muni_code, year) %>%
  dplyr::summarise(degrad_area = sum(degrad_area))

# merge
spatial.degrad.muniLevel <-
  samplePanel.muniLevel %>%
  dplyr::filter(year >= 2007 & year <= 2016) %>% # select DEGRAD time period
  dplyr::left_join(aux.degrad) %>% # merge with degrad data
  dplyr::mutate(degrad_area = replace_na(degrad_area, 0)) %>% # muni-years with NA values after merge means that it does not have any degradation polygon
  dplyr::left_join(sampleSpatial.muniLevel) %>%  # add spatial feature
  sf::st_as_sf() # define object as sf (spatial)

# clear environment
rm(clean.degrad, aux.degrad, sampleSpatial.muniLevel, samplePanel.muniLevel)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(spatial.degrad.muniLevel$year)        <- "year of reference (calendar or PRODES year)"
sjlabelled::set_label(spatial.degrad.muniLevel$degrad_area) <- "degradated area (ha, calculated from shapefile under SIRGAS2000 Polyconic projection)"



# POST-TREATMENT OVERVIEW
# summary(spatial.degrad.muniLevel)
# View(spatial.degrad.muniLevel)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(spatial.degrad.muniLevel,
     file = file.path("data/projectSpecific/muniLevel",
                      paste0("spatial_degrad_muniLevel", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/muniLevel")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
