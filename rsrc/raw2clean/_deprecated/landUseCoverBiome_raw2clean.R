
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CLEAN RAW BRAZILIAN AMAZON LAND USE DATA TO PROVIDE BIOME-BY-YEAR PANEL (MAPBIOMAS COLLECTION 7)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("code/setup.R")


# START TIMER
tictoc::tic(msg = "landUseCoverBiome_raw2clean.R script", log = T)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.mapbiomas <- readxl::read_xlsx(path = here::here("data/raw2clean/landUseCoverBiome_mapbiomas/input/1_-_TABELA_GERAL_COL7_MAPBIOMAS_BIOMAS_UF_FINAL.xlsx"),
                                   sheet = 2)



# DATA EXPLORATION
#summary(raw.mapbiomas)
#class(raw.mapbiomas)
#View(raw.mapbiomas)      # column names indicate file of origin





# DATASET CLEANUP AND PREP --------------------------------------------------------------------------------------------------------------------------

# FILTER AMAZON BIOME AND REMOVE UNNECESSARY COLUMNS
raw.mapbiomas <-
  raw.mapbiomas %>%
  dplyr::filter(biome == "AMAZÔNIA") %>%
  dplyr::select(class_id, as.character(1985:2021))

# adjust NAs to 0s
raw.mapbiomas[is.na(raw.mapbiomas)] <- 0



# AGGREGATE TO AMAZON-CLASS-YEAR
raw.mapbiomas <-
  raw.mapbiomas %>%
  dplyr::mutate(mapbiomas_classAgg = dplyr::case_when(class_id == 3 ~ "forest",
                                                      class_id %in% c(15, 20, 21, 39, 41, 48, 62) ~ "z",
                                                      class_id %in% c(0, 4, 5, 9, 11:13, 22:27, 29:33) ~ "otherCategories")) %>%
  dplyr::select(-class_id) %>%
  dplyr::group_by(mapbiomas_classAgg) %>%
  dplyr::summarise_all(sum)



# RESHAPE
raw.mapbiomas <-
  raw.mapbiomas %>%
  tidyr::pivot_longer(cols = "1985":"2021", names_to = "year", values_to = "mapbiomas_area") %>% # transform year columns to long format
  tidyr::pivot_wider(id_cols = c("year"), names_from = "mapbiomas_classAgg",
              values_from = "mapbiomas_area", values_fill = list(mapbiomas_area = 0)) %>% # transform mapbiomas_id rows to wide format filling NA values with 0
  dplyr::mutate(year = as.numeric(year)) # change year class from character to numeric





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
# change object name for exportation
clean.landUseCoverBiome <- raw.mapbiomas



# POST-TREATMENT OVERVIEW
# summary(clean.landUseCoverMuni)
# View(clean.landUseCoverMuni)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.landUseCoverBiome,
     file = here::here("data/raw2clean/landUseCoverBiome_mapbiomas/output",
                      "clean_landUseCoverBiome.Rdata"))


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("code/raw2clean")





# END OF SCRIPT: ------------------------------------------------------------------------------------------------------------------------------------