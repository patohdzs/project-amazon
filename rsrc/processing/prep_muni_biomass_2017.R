library(sf)
library(tictoc)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter())

# Start timer
tic(msg = "prep_muni_biomass_2017.R script", log = TRUE)

# Load municipalitality coordinates
load("data/clean/raw_muni.Rdata")

# Load pixel-level biomass data for 2017
load("data/processed/pixel_biomass_2017.Rdata")

# Transform to SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880),
# Join to municipalities data
muni_biomass_2017 <- pixel_biomass_2017 %>%
     st_transform(, crs = 5880) %>%
     st_join(raw_muni, join = st_within)

save(muni_biomass_2017, file = "data/processed/muni_biomass_2017.Rdata")

# End timer
toc(log = TRUE)
