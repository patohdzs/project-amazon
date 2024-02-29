
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: UNIFY TIF FILES BY YEAR - FLORESER (IMAZON - 2010)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "floreser2010_raw2clean script", log = T)



# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(raster) # for raster manipulation



# RASTER OPTIONS
raster::rasterOptions(tmpdir = file.path("data/_temp"),
                      timer  = T)





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RAW DATA
# read all parts (12) of the year and merge them together
raw.raster <- do.call(raster::merge, lapply(list.files(file.path("data/raw2clean/floreser_imazon/input/2010"),
                                                       full.names = T,
                                                       pattern = "Floreser_Biome_2010"),
                                            raster::raster))





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# RECLASSIFY
# change 0s to NAs so that any value represents only areas of secondary vegetation and its age
raw.raster <- raster::subs(raw.raster, data.frame(id = 0, v = NA), subsWithNA = F, progress = "text")





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

# EXPORT
# save unified tif
raster::writeRaster(raw.raster, "data/raw2clean/floreser_imazon/output/2010/clean_floreser.tif", overwrite = T)



# CLEAN TEMP DIR
# showTmpFiles()
raster::removeTmpFiles(h = 2)
gc()



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("raw2clean")






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------