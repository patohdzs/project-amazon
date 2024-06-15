
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




# START TIMER
tictoc::tic(msg = "mapbiomas_1043SitesModel.R script", log = TRUE)


# TERRA OPTIONS (specify temporary file location)
terra::terraOptions(tmpdir = "data/_temp",
                      timer  = T)




# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

mapbiomas_class <- c("forest", "agriculturalUse", "other")
aux_year <- c(1995, 2008,2017)

for (class in mapbiomas_class) {
  for (year in aux_year) {

  # DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

  # RAW DATA
  # read MapBiomas raster with land use and land cover (baseline)
  raw_raster <- terra::rast(glue::glue("data/raw/mapbiomas/land_use_cover/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-{year}.tif"))



  # CHANGE RASTER VALUES - see "documentation/mapbiomasClass_id_legend.pdf"
  if(class == "forest") {
    # create dummy for forest class
    raw_raster[raw_raster != 3] <- 0
    raw_raster[raw_raster == 3] <- 1
  } else if (class == "agriculturalUse") {
    # create dummy for agricultural Use classes
    raw_raster[!(raw_raster %in%  c(15, 20, 39, 41))] <- 0
    raw_raster[raw_raster %in% c(15, 20, 39, 41)] <- 1
  } else if (class == "other") {
    # create dummy for other classes (non-forest neither agricultural use)
    raw_raster[raw_raster %in%  c(3, 15, 20, 39, 41)] <- 0
    raw_raster[!(raw_raster %in%  c(0, 3, 15, 20, 39, 41))] <- 1
  }



  # aggregate raster calculating the share of minicells with the specified category
  raw_raster <- terra::aggregate(raw_raster, fact = 2250, fun = sum, na.rm = T)/(2250^2) # (2250^2) is the total number of minicells

  # add category name
  names(raw_raster) <- glue::glue("share_{class}_{year}")



  # EXPORT
  # save unified tif

  #terra::writeRaster(raw_raster, paste0("C:/Users/pengyu/Desktop/code_data_20230628/data/calibration/1055SitesModel/aux_tifs/raster_mapbiomas_", class, "_", year, "_1055SitesModel.tif"), overwrite = T)
  terra::writeRaster(raw_raster, glue::glue("data/calibration/1043SitesModel/aux_tifs/raster_mapbiomas_{class}_{year}_1043SitesModel.tif"),
                     overwrite = T)
  
  # clean environment
  rm(raw_raster)



  # CLEAN TEMP DIR
  terra::tmpFiles(current = T, remove = T)
  gc()

  }
}



# END TIMER
tictoc::toc(log = TRUE)






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------