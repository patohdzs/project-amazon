# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: TREAT RAW DATA - AMAZON BIOME BONUDARY (IBGE - 2019)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

# START TIMER
tictoc::tic(msg = "amazonBiome_raw2clean.R script", log = T)

# DATA INPUT

# read shapefile
raw.biome <- sf::st_read(
  dsn = here::here("data/raw2clean/amazonBiome_ibge/input"),
  layer = "lm_bioma_250"
)

# DATA EXPLORATION [disabled for speed]
# summary(raw.biome)
# View(raw.biome@data)
# plot(raw.biome)

# DATASET CLEANUP AND PREP

# COLUMN CLEANUP
# names
colnames(raw.biome)

# translate column names
raw.biome <-
  raw.biome %>%
  dplyr::rename(
    biome_code = CD_Bioma,
    biome_name = Bioma
  )

# class - no change needed
lapply(raw.biome, class)

# TRANSLATION
# 'grepl' used to avoid encoding trouble with latin characters
raw.biome$biome_name[which(grepl(pattern = "Amazônia", x = raw.biome$biome_name))] <- "Amazon"
raw.biome$biome_name[which(grepl(pattern = "Mata Atlântica", x = raw.biome$biome_name))] <- "Atlantic Forest"

# LETTERS CAPITALIZATION
raw.biome <-
  raw.biome %>%
  dplyr::mutate(biome_name = toupper(biome_name))

# FILTER BIOME OF INTEREST (AMAZON)
raw.biome <-
  raw.biome %>%
  dplyr::filter(biome_name == "AMAZON")


output_path <- here::here("data/hmc/", "map.geojson")
raw.biome2 <- sf::st_transform(x = raw.biome, crs = 4326) # SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)
st_write(raw.biome2, output_path, driver = "GeoJSON", delete_layer = TRUE)

# PROJECTION
raw.biome <- sf::st_transform(x = raw.biome, crs = 5880) # SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)

# GEOMETRY CLEANUP
raw.biome <- sf::st_make_valid(raw.biome)

# EXPORT PREP

# LABELS
sjlabelled::set_label(raw.biome$biome_code) <- "biome code"
sjlabelled::set_label(raw.biome$biome_name) <- "biome name"

# change object name for exportation
clean.amazonBiome <- raw.biome

# POST-TREATMENT OVERVIEW
# summary(clean.amazonBiome)
# View(clean.amazonBiome@data)
# plot(clean.amazonBiome$geometry)

# EXPORT

save(clean.amazonBiome,
  file = here::here(
    "data/raw2clean/amazonBiome_ibge/output",
    "clean_amazonBiome.Rdata"
  )
)

# # END TIMER
# tictoc::toc(log = T)

# # export time to csv table
# ExportTimeProcessing("code/raw2clean")
