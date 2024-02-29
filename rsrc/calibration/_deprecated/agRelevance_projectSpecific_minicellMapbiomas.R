
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: CREATE AGGREGATED CATEGORIES OF INTEREST BASED ON MAPBIOMAS CATEGORIES (AGRICULTURAL RELEVANCE)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "agRelevance_projectSpecific_minicellMapbiomas script", log = T)



# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(sf) # manipulate spatial data
library(tidyverse) # manipulate tables, works with sf
library(raster) # manipulate raster data
library(sp) # manipulate spatial data
library(rgeos) # manipulate spatial data
library(rgdal) # manipulate spatial data





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# Mapbiomas sample
load("data/projectSpecific/minicellMapbiomas/sample_minicellMapbiomas.Rdata")


# agricultural relevance variables - muni level
load("data/projectSpecific/muniLevel/spatial_agRelevance_muniLevel.Rdata")



# MAPBIOMAS - RASTER
clean.mapbiomas <- raster::raster("data/raw2clean/landUseCover_mapbiomas/input/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-2016.tif")




# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# add aggregated land use/cover categories (separate farming categories from others)
agRelevance.minicellMapbiomas <-
  sample.minicellMapbiomas %>%
  dplyr::filter(year == 2017) %>%
  dplyr::mutate(mapbiomas_agRelevance = dplyr::case_when(mapbiomas_class == 15 ~ "pasture",
                                                         mapbiomas_class == 39 ~ "soybean",
                                                         mapbiomas_class == 20 ~ "sugarcane",
                                                         mapbiomas_class == 41 ~ "otherTempCrops",
                                                         mapbiomas_class %in% c(1:5, 9:13, 22:26, 29:33) ~ "nonFarming",
                                                         mapbiomas_class %in% c(0, 27) ~ "nonObserved")) %>%
  dplyr::select(-year, -mapbiomas_class)

# clear environment
rm(sample.minicellMapbiomas)


# ADJUST MAPBIOMAS CATEGORIES (MORE AGGREGATED) OF THE ORIGINAL RASTER

# if file with aggregated categories does not exist yet, adjust categories and create it. If it exists just read it.
if (!file.exists("data/projectSpecific/minicellMapbiomas/aux_tifs/raster_mapbiomas_agg_2017.tif")) {

  # aggregate categories and save intermediary file
  clean.mapbiomas <- raster::subs(clean.mapbiomas,
                                  data.frame(id = c(c(0, 27), 15, 39, 20, 41, c(1:5, 9:13, 22:26, 29:33)),
                                             v  = c(rep(NA, 2), 1, 2, 3, 4, rep(5, 20))),
                                  subsWithNA = F,
                                  progress = "text",
                                  filename = "data/projectSpecific/minicellMapbiomas/aux_tifs/raster_mapbiomas_agg_2017.tif",
                                  overwrite = T)

} else {

  # read file with adjusted categories
  clean.mapbiomas <- raster::raster("data/projectSpecific/minicellMapbiomas/aux_tifs/raster_mapbiomas_agg_2017.tif")

}



# EXTRACT NEIGHBOR VALUES RELATIVE FREQUENCY

# add ID column to match after extraction
agRelevance.minicellMapbiomas$ID <- 1:nrow(agRelevance.minicellMapbiomas)

# transform mapbiomas sample into spatial points
sp.minicellMapbiomas <- sp::SpatialPointsDataFrame(coords = agRelevance.minicellMapbiomas[, 1:2], data = data.frame(agRelevance.minicellMapbiomas),
                                                                         proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))

# generate buffers around the points to include adjacent neighbors (width size represent 3x3 cells neighborhood - "queen")
sp.minicellMapbiomas <- rgeos::gBuffer(sp.minicellMapbiomas, byid = T, width = 0.0004)

# extract neighbor values by point/buffer of interest
aux.neighbors <- list() # create empty list to be filled in the loop

# LOOP (faster than running everything at once)
for (i in seq_along(sp.minicellMapbiomas$ID)) {

  aux.neighbors[[i]] <- raster::extract(clean.mapbiomas, sp.minicellMapbiomas[i,], df = T, na.rm = F)

}

# merge all data.frames
aux.neighbors <- do.call("rbind", aux.neighbors)

# clear enviroment
rm(sp.minicellMapbiomas, clean.mapbiomas)

# adjust output table and add relative neighbor frequency (including central point) to mapbiomas sample
agRelevance.minicellMapbiomas <-
  aux.neighbors %>%
  dplyr::group_by(ID, raster_mapbiomas_agg_2017) %>%
  dplyr::tally() %>%
  dplyr::mutate(share = n/sum(n)) %>%
  dplyr::select(-n) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(raster_mapbiomas_agg_2017 = dplyr::case_when(raster_mapbiomas_agg_2017 == 1 ~ "share_pasture",
                                                             raster_mapbiomas_agg_2017 == 2 ~ "share_soybean",
                                                             raster_mapbiomas_agg_2017 == 3 ~ "share_sugarcane",
                                                             raster_mapbiomas_agg_2017 == 4 ~ "share_otherTempCrops",
                                                             raster_mapbiomas_agg_2017 == 5 ~ "share_nonFarming",
                                                             is.na(raster_mapbiomas_agg_2017) ~ "share_nonObserved")) %>%
  tidyr::pivot_wider(id_cols = ID, names_from = raster_mapbiomas_agg_2017, values_from = share, values_fill = 0) %>%
  dplyr::right_join(agRelevance.minicellMapbiomas) %>%
  dplyr::select(lon, lat, mapbiomas_agRelevance, tidyselect:::starts_with("share"))
  # adjust relative freqency for cells that have nonObserved values (new_freq = sum(freqs)-freq_nonObserved)

# clear enviroment
rm(aux.neighbors)



# SPATIAL MERGE

# transform data to spatial points (sf)
agRelevance.minicellMapbiomas <- sf::st_as_sf(sp::SpatialPointsDataFrame(coords = agRelevance.minicellMapbiomas[, 1:2], data = data.frame(agRelevance.minicellMapbiomas),
                                                                         proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs")))

# construct and select variables of interest from muni level and adjust projection to match pixel level
spatial.agRelevance.muniLevel <-
  spatial.agRelevance.muniLevel %>%
  dplyr::mutate(shareOtherTempCrops_rice_muniLevel = planted_area_rice_pam/planted_area_3crops_pam,
                shareOtherTempCrops_corn_muniLevel = planted_area_cassava_pam/planted_area_3crops_pam,
                shareOtherTempCrops_cassava_muniLevel = planted_area_corn_pam/planted_area_3crops_pam) %>%
  dplyr::select(starts_with("shareOtherTempCrops")) %>%
  sf::st_transform(sf::st_crs(agRelevance.minicellMapbiomas))

# spatial merge
agRelevance.minicellMapbiomas <-
  agRelevance.minicellMapbiomas %>%
  sf::st_join(spatial.agRelevance.muniLevel, left = T) %>% # spatial merge
  dplyr::mutate(share_rice = shareOtherTempCrops_rice_muniLevel * share_otherTempCrops,
                share_corn = shareOtherTempCrops_corn_muniLevel * share_otherTempCrops,
                share_cassava = shareOtherTempCrops_cassava_muniLevel * share_otherTempCrops) %>% # create imputed variables
  dplyr::select(lon, lat, mapbiomas_agRelevance, starts_with("share_")) # select final variables of interest

# clear environment
rm(spatial.agRelevance.muniLevel)




# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(agRelevance.minicellMapbiomas$lon)                   <- "longitude of the cell centroid (degrees)"
sjlabelled::set_label(agRelevance.minicellMapbiomas$lat)                   <- "latitude of the cell centroid (degrees)"
sjlabelled::set_label(agRelevance.minicellMapbiomas$mapbiomas_agRelevance) <- "mapbiomas land use/cover aggregated classification (separate farming categories from others)"
sjlabelled::set_label(agRelevance.minicellMapbiomas$share_pasture)         <- "pasture relative frequency across adjacent neighbors ('queen' - 3x3 cells - include central cell)"
sjlabelled::set_label(agRelevance.minicellMapbiomas$share_soybean)         <- "soybean relative frequency across adjacent neighbors ('queen' - 3x3 cells - include central cell)"
sjlabelled::set_label(agRelevance.minicellMapbiomas$share_sugarcane)       <- "sugarcane relative frequency across adjacent neighbors ('queen' - 3x3 cells - include central cell)"
sjlabelled::set_label(agRelevance.minicellMapbiomas$share_otherTempCrops)  <- "other temporary crops relative frequency across adjacent neighbors ('queen' - 3x3 cells - include central cell)"
sjlabelled::set_label(agRelevance.minicellMapbiomas$share_nonFarming)      <- "non farming relative frequency across adjacent neighbors ('queen' - 3x3 cells - include central cell)"
sjlabelled::set_label(agRelevance.minicellMapbiomas$share_rice)            <- "rice relative frequency (imputed from muni level) across adjacent neighbors ('queen' - 3x3 cells - include central cell)"
sjlabelled::set_label(agRelevance.minicellMapbiomas$share_corn)            <- "corn relative frequency (imputed from muni level) across adjacent neighbors ('queen' - 3x3 cells - include central cell)"
sjlabelled::set_label(agRelevance.minicellMapbiomas$share_cassava)         <- "cassava relative frequency (imputed from muni level) across adjacent neighbors ('queen' - 3x3 cells - include central cell)"



# POST-TREATMENT OVERVIEW
# summary(agRelevance.minicellMapbiomas)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(agRelevance.minicellMapbiomas,
     file = file.path("data/projectSpecific/minicellMapbiomas",
                      paste0("agRelevance_minicellMapbiomas", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/minicellMapbiomas")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
