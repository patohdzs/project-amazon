# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CLEAN RAW AGRICULTURAL USE AREA - AGRICULTURAL CENSUS 2017 (IBGE)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

# SETUP

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")

# START TIMER
tictoc::tic(msg = "agCensus2017AgUseArea_raw2clean.R script", log = T)

# DATA INPUT

# read csv file
raw.agCensus2017AgUseArea <- readr::read_csv(
  file = here::here("data/raw2clean/agCensus2017AgUseArea_ibge/input/agCensus2017_agUseArea.csv"),
  skip = 7,
  na = c("..."),
  col_names = c(
    "muni_code",
    "cropPerm_area_2017", "cropTemp_area_2017", "cropFlower_area_2017",
    "pastureNatural_area_2017", "pasturePlantedGood_area_2017", "pasturePlantedBad_area_2017"
  ),
  col_types = "n--cccccc"
)

# DATA EXPLORATION [disabled for speed]
# summary(raw.agCensus2017AgUseArea)
# View(raw.agCensus2017AgUseArea)

# DATASET CLEANUP AND PREP

# ROW CLEANUP
# remove last rows of the table with notes information
raw.agCensus2017AgUseArea <- raw.agCensus2017AgUseArea[-(5564:5577), ]

# transform "-" values to "0" as explained in the table notes
raw.agCensus2017AgUseArea <-
  raw.agCensus2017AgUseArea %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# transform "X" values to "NA" (here NA identify which values had to be omitted to avoid informant identification as explained in the table notes)
raw.agCensus2017AgUseArea <-
  raw.agCensus2017AgUseArea %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))

# latin character treatment
raw.agCensus2017AgUseArea <-
  raw.agCensus2017AgUseArea %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), \(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")))

# transform column class
raw.agCensus2017AgUseArea <-
  raw.agCensus2017AgUseArea %>%
  dplyr::mutate(dplyr::across(tidyselect:::ends_with("_area_2017"), function(x) as.numeric(x)))

# sum pasture, crop, and agricultural use area
raw.agCensus2017AgUseArea <-
  raw.agCensus2017AgUseArea %>%
  dplyr::mutate(
    pasture_area_2017 = rowSums(across(c("pastureNatural_area_2017", "pasturePlantedGood_area_2017", "pasturePlantedBad_area_2017")), na.rm = TRUE),
    crop_area_2017 = rowSums(across(c("cropPerm_area_2017", "cropTemp_area_2017", "cropFlower_area_2017")), na.rm = TRUE),
    agUse_area_2017 = rowSums(across(c("pasture_area_2017", "crop_area_2017")), na.rm = TRUE)
  )

# EXPORT PREP

# sjlabelled::set_labelS
sjlabelled::set_label(raw.agCensus2017AgUseArea$muni_code) <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.agCensus2017AgUseArea$pasture_area_2017) <- "pasture area (ha, 2017 Ag Census)"
sjlabelled::set_label(raw.agCensus2017AgUseArea$crop_area_2017) <- "crop area (ha, 2017 Ag Census)"
sjlabelled::set_label(raw.agCensus2017AgUseArea$agUse_area_2017) <- "agricultural use area (ha, 2017 Ag Census)"

# change object name for exportation
clean.agCensus2017AgUseArea <- raw.agCensus2017AgUseArea

# POST-TREATMENT OVERVIEW
# summary(clean.agCensus2017AgUseArea)
# View(clean.agCensus2017AgUseArea)

# EXPORT

save(clean.agCensus2017AgUseArea,
  file = here::here(
    "data/raw2clean/agCensus2017AgUseArea_ibge/output",
    "clean_agCensus2017AgUseArea.Rdata"
  )
)

# END TIMER
tictoc::toc(log = T)

# export time to csv table
# ExportTimeProcessing("code/raw2clean")
