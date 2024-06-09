
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


# START TIMER
tictoc::tic(msg = "sampleConstructionPixel_prepData.R script", log = TRUE)


# RASTER OPTIONS
terra::terraOptions(tmpdir = "data/_temp",
                    timer  = T)


# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# RAW DATA
clean_mapbiomas <- terra::rast("data/clean/landusecover_2000.tif")


# DATASET MANIPULATION -------------------------------------------------------------------------------------------------------------------------------

# RANDOM SAMPLE

# guarantee reproducibility
set.seed(123)

# extract random sample cells and keep only lon lat info
samplePixel_prepData <- data.frame(terra::spatSample(clean_mapbiomas, 1200000, na.rm = T, xy = T)[, 1:2])

# read all mapBiomas layers (one for each year)
aux_mapbiomas <- terra::rast(list.files("data/raw/mapbiomas/land_use_cover/",
                                        pattern = "COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-\\d{4}",
                                        full.names = T))

# change raster stack layer names
names(aux_mapbiomas) <- c(1985:2019)

# use points sample to extract complete panel
samplePixel_prepData <- cbind(samplePixel_prepData, terra::extract(aux_mapbiomas, samplePixel_prepData, df = T)[, -1])

# reshape
samplePixel_prepData <- tidyr::pivot_longer(samplePixel_prepData,
                                            cols = tidyselect:::matches("\\d{4}"), names_to = "year", values_to = "mapbiomas_class")

# minor fix year column
samplePixel_prepData$year <- as.numeric(stringr::str_extract(samplePixel_prepData$year, "\\d{4}"))

# rename columns
samplePixel_prepData <-
  samplePixel_prepData %>%
  dplyr::rename(lon = x, lat = y)



# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()



# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(samplePixel_prepData$lon)                <- "longitude of the pixel centroid (degrees)"
sjlabelled::set_label(samplePixel_prepData$lat)                <- "latitude of the pixel centroid (degrees)"
sjlabelled::set_label(samplePixel_prepData$year)               <- "year of reference"
sjlabelled::set_label(samplePixel_prepData$mapbiomas_class)    <- "mapbiomas land use/cover classification (id)"



# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(samplePixel_prepData,
     file = "data/prepData/samplePixel_prepData.Rdata")

# END TIMER
tictoc::toc(log = TRUE)


# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------