

# START TIMER
tictoc::tic(msg = "merge_muni_gamma.R script", log = TRUE)

# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# Municipalitalities coordinates
load("data/clean/raw_muni.Rdata")

# Pixel gammas
load("data/prepData/pixelBiomass2017_prepData.Rdata")


# transform to Brazil Polyconic
sfdata_2017 <- sf::st_transform(x = pixelBiomass2017_prepData, crs = 5880) # SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)
gamma_muni_2017 <- st_join(sfdata_2017, raw_muni, join = st_within)


save(gamma_muni_2017,
     file = "data/prepData/muni_Biomass2017_prepData.Rdata")



# END TIMER
tictoc::toc(log = TRUE)



# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------