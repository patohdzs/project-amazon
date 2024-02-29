
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: RESTRICT SOYBEAN POTENTIAL YIELD DATA TO AMAZON BIOME (FAO-GAEZ)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1:




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
here::i_am("code/raw2clean/precipitation_raw2clean.R", uuid = "c7c9648a-396f-43e8-89f3-f6d86140434c")


# START TIME
tictoc::tic(msg = "precipitation_raw2clean script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# LIBRARIES
groundhog::groundhog.library("terra", groundhog.date)  # manipulate spatial data (raster format)
groundhog::groundhog.library("sf", groundhog.date)  # manipulate spatial data (vector format)


# RASTER OPTIONS
terra::terraOptions(tmpdir = here::here("data/_temp"),
                      timer  = T)





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RAW DATA
raw.raster <- terra::rast(here::here("data/raw2clean/potentialYield_faogaez/input/soyb200b_yld.tif"))

# AUX DATA (AMAZON BIOME BOUNDARY)
load(here::here("data/raw2clean/amazonBiome_ibge/output/clean_amazonBiome.Rdata"))





# DATA PROCESSING ------------------------------------------------------------------------------------------------------------------------------------

# change projection to match with raster
clean.amazonBiome <- sf::st_transform(clean.amazonBiome, 4326)

# crop to biome extent
raw.raster <- terra::crop(raw.raster, clean.amazonBiome)

# remove values outside the biome (change to NA)
raw.raster <- terra::mask(raw.raster, terra::vect(clean.amazonBiome))





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

# EXPORT
# save unified tif
terra::writeRaster(raw.raster, here::here("data/raw2clean/potentialYield_faogaez/output/clean_potentialYield.tif"), overwrite = T)




# CLEAN TEMP DIR
# showTmpFiles()
terra::tmpFiles(remove = T)
gc()



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------