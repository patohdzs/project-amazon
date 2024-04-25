
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: GENERATE AGGREGATED MAPBIOMAS VARIABLES (FOREST, AGRICULTURAL USE, OTHER) - 1055 SITES
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("code/setup.R")


# START TIMER
tictoc::tic(msg = "mapbiomas_1055SitesModel.R script", log = T)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tempdir = here::here("data", "_temp"))




# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

mapbiomas.class <- c("forest", "agriculturalUse", "other")
aux.year <- c(1995, 2008)

for (class in mapbiomas.class) {
  for (year in aux.year) {

  # DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

  # RAW DATA
  # read MapBiomas raster with land use and land cover (baseline)
  raw.raster <- terra::rast(here::here(glue::glue("data/raw2clean/landUseCover_mapbiomas/input/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-{year}.tif")))



  # CHANGE RASTER VALUES - see "documentation/mapbiomasClass_id_legend.pdf"
  if(class == "forest") {
    # create dummy for forest class
    raw.raster[raw.raster != 3] <- 0
    raw.raster[raw.raster == 3] <- 1
  } else if (class == "agriculturalUse") {
    # create dummy for agricultural Use classes
    raw.raster[!(raw.raster %in%  c(15, 20, 39, 41))] <- 0
    raw.raster[raw.raster %in% c(15, 20, 39, 41)] <- 1
  } else if (class == "other") {
    # create dummy for other classes (non-forest neither agricultural use)
    raw.raster[raw.raster %in%  c(3, 15, 20, 39, 41)] <- 0
    raw.raster[!(raw.raster %in%  c(0, 3, 15, 20, 39, 41))] <- 1
  }



  # aggregate raster calculating the share of minicells with the specified category
  #raw.raster <- terra::aggregate(raw.raster, fact = 2250, fun = sum, na.rm = T)/(2250^2) # (2250^2) is the total number of minicells
  raw.raster <- terra::aggregate(raw.raster, fact = 100, fun = sum, na.rm = T)/(100^2) 
  

  # add category name
  names(raw.raster) <- glue::glue("share_{class}_{year}")



  # EXPORT
  # save unified tif
  terra::writeRaster(raw.raster, here::here(glue::glue("data/calibration/1055SitesModel/aux_tifs2/raster_mapbiomas_{class}_{year}_1055SitesModel.tif")),
                    overwrite = T)

  #terra::writeRaster(raw.raster, paste0("C:/Users/pengyu/Desktop/code_data_20230628/data/calibration/1055SitesModel/aux_tifs2/raster_mapbiomas_", class, "_", year, "_1055SitesModel.tif"), overwrite = T)
  
  
  # clean environment
  rm(raw.raster)



  # CLEAN TEMP DIR
  terra::tmpFiles(current = T, remove = T)
  gc()

  }
}

stop()

# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("code/calibration")






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------