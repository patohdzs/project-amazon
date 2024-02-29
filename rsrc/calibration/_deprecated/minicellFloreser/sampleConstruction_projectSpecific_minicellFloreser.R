
# > PROJECT INFO
# NAME: INCENTIVES AMAZON
# LEAD: JULIANO ASSUNCAO, JOSE SCHEINKMAN, AND LARS HANSEN
#
# > THIS SCRIPT
# AIM: EXTRACT RANDOM SAMPLE OF SECONDARY VEGETATION CELLS (1% - YEAR 2000)
# AUTHOR: JOAO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# START TIME
tictoc::tic(msg = "sampleConstruction_projectSpecific_minicellFloreser script", log = T)



# SOURCES
source("code/_functions/ExportTimeProcessing.R")



# LIBRARIES
library(raster)     # for raster manipulation
library(data.table) # for efficient data frame manipulation


# RASTER OPTIONS
raster::rasterOptions(tmpdir = file.path("data/_temp"),
                      timer  = T)





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RAW DATA
clean.floreser <- raster::raster("data/raw2clean/floreser_imazon/output/2000/clean_floreser.tif")





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# RANDOM SAMPLE (1% of sec veg cells)

# create empty data.table to be filled in the loop
sample.minicellFloreser <- data.table::setDT(data.frame(x = NA, y = NA, clean_floreser = NA)[-1,])

# guarantee reproducibility
set.seed(123)



# LOOP
# need to do multiple random sample because we only want sec veg cells (non-NA values) and when there is a lot of NA values sampleRandom >
# return less observations, also because we have a large raster trying to extract a lot of cells at once may consume all the RAM memory.
while (nrow(sample.minicellFloreser) <= 739200) {

  # extract random sample of sec veg cells (values > 0) and transform it to data.table
  aux.sample <- data.table::setDT(data.frame((raster::sampleRandom(clean.floreser, 10000, na.rm = T, xy = T))))

  # merge with other samples
  sample.minicellFloreser <- data.table::rbindlist(list(sample.minicellFloreser, aux.sample))

  # remove possible duplicated cells
  sample.minicellFloreser <- unique(sample.minicellFloreser)

}


# EXPORT
# save extracted random sample
save(sample.minicellFloreser, file = paste0("data/projectSpecific/minicellFloreser/sample_minicellFloreser.Rdata"))



# CLEAN TEMP DIR
# showTmpFiles()
raster::removeTmpFiles(h = 2)
gc()



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("projectSpecific/minicellFloreser")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------