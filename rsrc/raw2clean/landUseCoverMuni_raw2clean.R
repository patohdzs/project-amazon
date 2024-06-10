# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CLEAN RAW BRAZILIAN AMAZON LAND USE DATA TO PROVIDE MUNICIPALITY-BY-YEAR PANEL (MAPBIOMAS COLLECTION 5)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# -
library(tidyverse)
library(tictoc)
library(sjlabelled)
library(conflicted)
library(readxl)

# Resolve conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

# START TIMER
tic(msg = "landUseCoverMuni_raw2clean.R script", log = TRUE)

# Read in data
raw_mapbiomas <- read_xlsx(
  path = "data/raw/mapbiomas/landuse_cover_muni/Dados_Cobertura_MapBiomas_5.0_UF-MUN_SITE_v2.xlsx",
  sheet = 3
)

# Clean column names
raw_mapbiomas <-
  raw_mapbiomas %>%
  rename(
    state_uf = state,
    muni_code = territory_id,
    muni_name = municipality
  )

# Reshape
raw_mapbiomas <-
  raw_mapbiomas %>%
  select(-muni_name, -state_uf) %>% # remove unnecessary columns
  unite("mapbiomas_class", level_0, level_1, level_2, level_3, level_4, sep = "_") %>%
  mutate(mapbiomas_id = case_when(
    mapbiomas_class == "Anthropic_3 - Farming_Agriculture_Temporary Crops_Mosaic of Crops" ~ "mapbiomasLandCoverId_41",
    mapbiomas_class == "Anthropic_3 - Farming_Agriculture_Temporary Crops_Soy Beans" ~ "mapbiomasLandCoverId_39",
    mapbiomas_class == "Anthropic_3 - Farming_Pasture_Pasture_Pasture" ~ "mapbiomasLandCoverId_15",
    mapbiomas_class == "Anthropic_4 - Non Vegetated Area_Urban Infrastructure_Urban Infrastructure_Urban Infrastructure" ~ "mapbiomasLandCoverId_24",
    mapbiomas_class == "Natural_1 - Forest_Natural Forest_Forest Formation_Forest Formation" ~ "mapbiomasLandCoverId_3",
    mapbiomas_class == "Natural_2 - Non Forest Natural Formation_Grassland_Grassland_Grassland" ~ "mapbiomasLandCoverId_12",
    mapbiomas_class == "Natural_5 - Water_River, Lake and Ocean_River, Lake and Ocean_River, Lake and Ocean" ~ "mapbiomasLandCoverId_33",
    mapbiomas_class == "Not applied_Non Observed_Non Observed_Non Observed_Non Observed" ~ "mapbiomasLandCoverId_27",
    mapbiomas_class == "Natural_1 - Forest_Natural Forest_Savanna Formation_Savanna Formation" ~ "mapbiomasLandCoverId_4",
    mapbiomas_class == "Anthropic_4 - Non Vegetated Area_Other Non Vegetated Area_Other Non Vegetated Area_Other Non Vegetated Area" ~ "mapbiomasLandCoverId_25",
    mapbiomas_class == "Anthropic_1 - Forest_Forest Plantation_Forest Plantation_Forest Plantation" ~ "mapbiomasLandCoverId_9",
    mapbiomas_class == "Anthropic_4 - Non Vegetated Area_Mining_Mining_Mining" ~ "mapbiomasLandCoverId_30",
    mapbiomas_class == "Anthropic_3 - Farming_Agriculture_Temporary Crops_Sugar Cane" ~ "mapbiomasLandCoverId_20",
    mapbiomas_class == "Natural_1 - Forest_Natural Forest_Magrove_Magrove" ~ "mapbiomasLandCoverId_5",
    mapbiomas_class == "Natural_2 - Non Forest Natural Formation_Salt flat_Salt flat_Salt flat" ~ "mapbiomasLandCoverId_32",
    mapbiomas_class == "Natural_4 - Non Vegetated Area_Beach and Dune_Beach and Dune_Beach and Dune" ~ "mapbiomasLandCoverId_23",
    mapbiomas_class == "Natural_2 - Non Forest Natural Formation_Wetland_Wetland_Wetland" ~ "mapbiomasLandCoverId_11",
    mapbiomas_class == "Anthropic_3 - Farming_Agriculture_Perennial Crops_Perennial Crops" ~ "mapbiomasLandCoverId_36",
    mapbiomas_class == "Anthropic_3 - Farming_Mosaic of Agriculture and Pasture_Mosaic of Agriculture and Pasture_Mosaic of Agriculture and Pasture" ~ "mapbiomasLandCoverId_21",
    mapbiomas_class == "Natural_2 - Non Forest Natural Formation_Rocky outcrop_Rocky outcrop_Rocky outcrop" ~ "mapbiomasLandCoverId_29",
    mapbiomas_class == "Anthropic_5 - Water_Aquaculture_Aquaculture_Aquaculture" ~ "mapbiomasLandCoverId_31",
    mapbiomas_class == "Natural_2 - Non Forest Natural Formation_Other Non Forest Natural Formation_Other Non Forest Natural Formation_Other Non Forest Natural Formation" ~ "mapbiomasLandCoverId_13"
  )) %>% # add mapbiomas id to simplify future column names - see "_PT-BR__Códigos_da_legenda_Coleção_5" on documentation folder
  select(-mapbiomas_class) %>% # remove mapbiomas_class column
  pivot_longer(
    cols = "1985":"2019",
    names_to = "year",
    values_to = "mapbiomas_id_area"
  ) %>% # transform year columns to long format
  pivot_wider(
    id_cols = c("muni_code", "year"),
    names_from = "mapbiomas_id",
    values_from = "mapbiomas_id_area",
    values_fill = list(mapbiomas_id_area = 0)
  ) %>% # transform mapbiomas_id rows to wide format filling NA values with 0
  mutate(year = as.numeric(year)) # change year class from character to numeric


# EXPORT PREP
# set_labelS
set_label(raw_mapbiomas$muni_code) <- "municipality code (7-digit, IBGE)"
set_label(raw_mapbiomas$year) <- "year of reference (calendar or PRODES year)"
set_label(raw_mapbiomas$mapbiomasLandCoverId_41) <- "(calendar year) area (ha) mapbiomas category  = anthropic_farming_agriculture_temporaryCrops_otherTemporaryCrops"
set_label(raw_mapbiomas$mapbiomasLandCoverId_39) <- "(calendar year) area (ha) mapbiomas category  = anthropic_farming_agriculture_temporaryCrops_soybeans"
set_label(raw_mapbiomas$mapbiomasLandCoverId_15) <- "(calendar year) area (ha) mapbiomas category  = anthropic_farming_pasture"
set_label(raw_mapbiomas$mapbiomasLandCoverId_24) <- "(calendar year) area (ha) mapbiomas category  = anthropic_nonvevegetatedArea_UrbanInfrastructure"
set_label(raw_mapbiomas$mapbiomasLandCoverId_3) <- "(calendar year) area (ha) mapbiomas category  = natural_forest_naturalForest_forestFormation"
set_label(raw_mapbiomas$mapbiomasLandCoverId_12) <- "(calendar year) area (ha) mapbiomas category  = natural_nonForestNaturalFormation_grassland"
set_label(raw_mapbiomas$mapbiomasLandCoverId_33) <- "(calendar year) area (ha) mapbiomas category  = natural_water_riverLakeOcean"
set_label(raw_mapbiomas$mapbiomasLandCoverId_27) <- "(calendar year) area (ha) mapbiomas category  = notApplied_nonObserved"
set_label(raw_mapbiomas$mapbiomasLandCoverId_4) <- "(calendar year) area (ha) mapbiomas category  = natural_forest_naturalForest_savannaFormatio"
set_label(raw_mapbiomas$mapbiomasLandCoverId_25) <- "(calendar year) area (ha) mapbiomas category  = anthropic_nonVegetatedArea_otherNonVegetatedAreas"
set_label(raw_mapbiomas$mapbiomasLandCoverId_9) <- "(calendar year) area (ha) mapbiomas category  = anthropic_forest_forestPlantation"
set_label(raw_mapbiomas$mapbiomasLandCoverId_30) <- "(calendar year) area (ha) mapbiomas category  = anthropic_nonVegetatedArea_mining"
set_label(raw_mapbiomas$mapbiomasLandCoverId_20) <- "(calendar year) area (ha) mapbiomas category  = anthropic_farming_agriculture_temporaryCrops_sugarcane"
set_label(raw_mapbiomas$mapbiomasLandCoverId_5) <- "(calendar year) area (ha) mapbiomas category  = natural_forest_naturalForest_mangrove"
set_label(raw_mapbiomas$mapbiomasLandCoverId_32) <- "(calendar year) area (ha) mapbiomas category  = natural_nonForestNaturalFormation_saltFlat"
set_label(raw_mapbiomas$mapbiomasLandCoverId_23) <- "(calendar year) area (ha) mapbiomas category  = natural_nonVegetatedArea_beachAndDune"
set_label(raw_mapbiomas$mapbiomasLandCoverId_11) <- "(calendar year) area (ha) mapbiomas category  = natural_nonForestNaturalFormation_wetland"
set_label(raw_mapbiomas$mapbiomasLandCoverId_36) <- "(calendar year) area (ha) mapbiomas category  = anthropic_farming_agriculture_perennialCrops"
set_label(raw_mapbiomas$mapbiomasLandCoverId_21) <- "(calendar year) area (ha) mapbiomas category  = anthropic_farming_mosaicOfAgriculturaAndPasture"
set_label(raw_mapbiomas$mapbiomasLandCoverId_29) <- "(calendar year) area (ha) mapbiomas category  = natural_nonForestNaturalFormation_rockyOutcrop"
set_label(raw_mapbiomas$mapbiomasLandCoverId_31) <- "(calendar year) area (ha) mapbiomas category  = anthropic_water_aquaculture"
set_label(raw_mapbiomas$mapbiomasLandCoverId_13) <- "(calendar year) area (ha) mapbiomas category  = natural_nonForestNaturalFormation_otherNonForestFormations"

# Change object name before saving
land_use_cover_muni <- raw_mapbiomas

# Save data set
out_path <- "data/clean/land_use_cover_muni.Rdata"
save(land_use_cover_muni, file = out_path)

# END TIMER
toc(log = TRUE)
