

source("rsrc/setup.R")

# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# Municipalitalities coordinates
load(here::here("data/raw2clean/muniDivision2015_ibge/output/raw_muni.Rdata"))

# Pixel gammas
load(here::here("data/calibration/prepData/pixelBiomass2017_prepData.Rdata"))


# transform to Brazil Polyconic
sfdata_2017 <- sf::st_transform(x = pixelBiomass2017.prepData, crs = 5880) # SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)
gamma_muni_2017 <- st_join(sfdata_2017, raw.muni, join = st_within)


# file_path <- paste(getwd(), "data/calibration", paste0("muni_Biomass2017_prepData", ".Rdata"), sep = "/")


# save(gamma_muni_2017, file = file_path)

save(gamma_muni_2017,
     file = here::here("data/calibration/prepData",
                       paste0("muni_Biomass2017_prepData", ".Rdata")))


