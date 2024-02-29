
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: EXTRACT RANDOM SAMPLE OF MAPBIOMAS 30M-PIXELS AND RECOVER FULL PANEL (1985-2019)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -




# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("code/setup.R")


# START TIMER
tictoc::tic(msg = "sampleConstructionPixel_prepData.R script", log = T)


# SOURCE FUNCTIONS
source(here::here("code/_functions/ExportTimeProcessing.R"))


# RASTER OPTIONS
terra::terraOptions(tmpdir = here::here("data/_temp"),
                    timer  = T)





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RAW DATA
clean.mapbiomas <- terra::rast(here::here("data/raw2clean/landUseCover_mapbiomas/output/clean_landUseCover_2000.tif"))





# DATASET MANIPULATION -------------------------------------------------------------------------------------------------------------------------------

# RANDOM SAMPLE

# guarantee reproducibility
set.seed(123)

# extract random sample cells and keep only lon lat info
samplePixel.prepData <- data.frame(terra::spatSample(clean.mapbiomas, 1200000, na.rm = T, xy = T)[, 1:2])

# read all mapBiomas layers (one for each year)
aux.mapbiomas <- terra::rast(list.files(here::here("data/raw2clean/landUseCover_mapbiomas/input/"),
                                        pattern = "COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-\\d{4}",
                                        full.names = T))

# change raster stack layer names
names(aux.mapbiomas) <- c(1985:2019)

# use points sample to extract complete panel
samplePixel.prepData <- cbind(samplePixel.prepData, terra::extract(aux.mapbiomas, samplePixel.prepData, df = T)[, -1])

# reshape
samplePixel.prepData <- tidyr::pivot_longer(samplePixel.prepData,
                                            cols = tidyselect:::matches("\\d{4}"), names_to = "year", values_to = "mapbiomas_class")

# minor fix year column
samplePixel.prepData$year <- as.numeric(stringr::str_extract(samplePixel.prepData$year, "\\d{4}"))

# rename columns
samplePixel.prepData <-
  samplePixel.prepData %>%
  dplyr::rename(lon = x, lat = y)



# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(samplePixel.prepData$lon)                <- "longitude of the pixel centroid (degrees)"
sjlabelled::set_label(samplePixel.prepData$lat)                <- "latitude of the pixel centroid (degrees)"
sjlabelled::set_label(samplePixel.prepData$year)               <- "year of reference"
sjlabelled::set_label(samplePixel.prepData$mapbiomas_class)    <- "mapbiomas land use/cover classification (id)"


# POST-TREATMENT OVERVIEW
# summary(samplePixel.prepData)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(samplePixel.prepData,
     file = here::here("data/calibration/prepData",
                       paste0("samplePixel_prepData", ".Rdata")))



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("code/calibration")






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------