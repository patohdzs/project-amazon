
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: DOWNLOAD ABOVERGROUND BIOMASS/CARBON DATA (BACCINI - GLOBAL FOREST WATCH - 2000)
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
here::i_am("code/raw2clean/abovegroundBiomass_download.R", uuid = "351b1938-6de6-4ab5-a226-02200f225232")


# START TIME
tictoc::tic(msg = "abovegroundBiomass_download script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("tidyverse", groundhog.date)  # manipulate tables, works with sf
groundhog::groundhog.library("sf", groundhog.date)  # manipulate spatial data (vector format)
groundhog::groundhog.library("raster", groundhog.date)  # manipulate spatial data (raster format)
groundhog::groundhog.library("sjlabelled", groundhog.date) # label columns, preferred than Hmisc::label because has function to clear labels when necessary





# DATA DOWNLOAD --------------------------------------------------------------------------------------------------------------------------------------

# download process is optional given that the data is provided
if (stringr::str_detect(list.files(here::here("data/raw2clean/abovegroundBiomass_gfw/input")), "_t_aboveground_biomass_ha_2000.tif") %>% sum() != 12) {

  # read aux_download data (shapefile with download links by tile)
  aux.download <- sf::st_read(dsn = here::here("data/raw2clean/abovegroundBiomass_gfw/input/aux_download"),
                              layer = "Aboveground_live_woody_biomass_density")

  # load amazon biome shapefile to have the region of interest
  load(file = here::here("data/raw2clean/amazonBiome_ibge/output/clean_amazonBiome.Rdata"))


  # project amazon biome to same crs of aux.download
  clean.amazonBiome <- sf::st_transform(clean.amazonBiome, sf::st_crs(aux.download))

  # crop aux_downnload by the region of interest (Amazon Biome)
  aux.download <- sf::st_crop(aux.download, clean.amazonBiome)

  # visual cross-check
  #plot(clean.amazonBiome$geometry)
  #plot(aux.download$geometry, add = T)

  # use links to download data
  purrr::map2(.x = aux.download$download,
              .y = here::here(paste0("data/raw2clean/abovegroundBiomass_gfw/input/", aux.download$name, ".tif")),
              .f =  download.file, mode = "wb")

}



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------