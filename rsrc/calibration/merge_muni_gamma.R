

load("data/calibration/raw_muni.Rdata")
# 2010
load("data/calibration/prepData/pixelBiomass2017_prepData.Rdata")

sfdata_2017 <- sf::st_transform(x = pixelBiomass2017.prepData, crs = 5880) # SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)
gamma_muni_2017 <- st_join(sfdata_2017, raw.muni, join = st_within)

file_path <- paste(getwd(), "data/calibration", paste0("muni_Biomass2017_prepData", ".Rdata"), sep = "/")


save(gamma_muni_2017, file = file_path)




