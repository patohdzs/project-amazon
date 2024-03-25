
load("data/calibration/raw_muni.Rdata")

raster.variables <- terra::rast(list.files(here::here("data/calibration/1055SitesModel/aux_tifs2"),
                                           pattern = "raster_mapbiomas",
                                           full.names = T))
z_muni <- terra::as.polygons(raster.variables, dissolve = F) %>% sf::st_as_sf() %>% sf::st_transform(5880)
z_muni_2017 <- st_join(z_muni, raw.muni, join = st_within)



file_path <- paste(getwd(), "data/calibration", paste0("z_muni_2017", ".Rdata"), sep = "/")

save(z_muni_2017, file = file_path)