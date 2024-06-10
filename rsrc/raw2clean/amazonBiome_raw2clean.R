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
tictoc::tic(msg = "amazonBiome_raw2clean.R script", log = TRUE)

# DATA INPUT

# read shapefile
raw_biome <- sf::st_read(
  dsn = "data/raw/ibge/amazon_biome/",
  layer = "lm_bioma_250"
)


# Column names
colnames(raw_biome)

# Translate column names
raw_biome <-
  raw_biome %>%
  dplyr::rename(
    biome_code = CD_Bioma,
    biome_name = Bioma
  )

# Class - no change needed
lapply(raw_biome, class)

# TRANSLATION
# 'grepl' used to avoid encoding trouble with latin characters
raw_biome$biome_name[which(grepl(pattern = "Amazônia", x = raw_biome$biome_name))] <- "Amazon"
raw_biome$biome_name[which(grepl(pattern = "Mata Atlântica", x = raw_biome$biome_name))] <- "Atlantic Forest"

# LETTERS CAPITALIZATION
raw_biome <-
  raw_biome %>%
  dplyr::mutate(biome_name = toupper(biome_name))

# FILTER BIOME OF INTEREST (AMAZON)
raw_biome <-
  raw_biome %>%
  dplyr::filter(biome_name == "AMAZON")


output_path <- "data/calibration/hmc/map.geojson"
biome_output <- sf::st_transform(x = raw_biome, crs = 4326) # SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)
st_write(biome_output, output_path, driver = "GeoJSON", delete_layer = TRUE)

# PROJECTION
# SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)
raw_biome <- sf::st_transform(x = raw_biome, crs = 5880)

# GEOMETRY CLEANUP
raw_biome <- sf::st_make_valid(raw_biome)

# LABELS
sjlabelled::set_label(raw_biome$biome_code) <- "biome code"
sjlabelled::set_label(raw_biome$biome_name) <- "biome name"

# change object name for exportation
clean_amazon_biome <- raw_biome

# Save data set
save(clean_amazon_biome,
  file = "data/clean/amazon_biome.Rdata"
)

# END TIMER
tictoc::toc(log = TRUE)
