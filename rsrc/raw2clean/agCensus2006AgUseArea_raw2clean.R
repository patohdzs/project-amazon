# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CLEAN RAW AGRICULTURAL USE AREA - AGRICULTURAL CENSUS 2006 (IBGE)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -


# START TIMER
tictoc::tic(msg = "agCensus2006AgUseArea_raw2clean.R script", log = TRUE)

# DATA INPUT
# read csv file
raw_agCensus2006AgUseArea <- readr::read_csv(
  file = "data/raw/ibge/ag_census_2006_ag_use_area/agCensus2006_agUseArea.csv",
  skip = 6,
  na = c("..."),
  col_names = c(
    "muni_code",
    "cropPerm_area_2006", "cropTemp_area_2006", "cropCut_area_2006", "cropFlower_area_2006",
    "pastureNatural_area_2006", "pasturePlantedGood_area_2006", "pasturePlantedBad_area_2006"
  ),
  col_types = "n---cccccc"
)

# DATASET CLEANUP AND PREP
# ROW CLEANUP
# remove last rows of the table with notes information
raw_agCensus2006AgUseArea <- raw_agCensus2006AgUseArea[-(5549:5561), ]

# transform "-" values to "0" as explained in the table notes
raw_agCensus2006AgUseArea <-
  raw_agCensus2006AgUseArea %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# transform "X" values to "NA" (here NA identify which values had to be omitted to avoid informant identification as explained in the table notes)
raw_agCensus2006AgUseArea <-
  raw_agCensus2006AgUseArea %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))

# latin character treatment
raw_agCensus2006AgUseArea <-
  raw_agCensus2006AgUseArea %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), \(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")))

# transform column class
raw_agCensus2006AgUseArea <-
  raw_agCensus2006AgUseArea %>%
  dplyr::mutate(dplyr::across(tidyselect:::ends_with("_area_2006"), function(x) as.numeric(x)))

# sum pasture, crop, and agricultural use area
raw_agCensus2006AgUseArea <-
  raw_agCensus2006AgUseArea %>%
  dplyr::mutate(
    pasture_area_2006 = rowSums(across(c("pastureNatural_area_2006", "pasturePlantedGood_area_2006", "pasturePlantedBad_area_2006")), na.rm = TRUE),
    crop_area_2006 = rowSums(across(c("cropPerm_area_2006", "cropTemp_area_2006", "cropCut_area_2006", "cropFlower_area_2006")), na.rm = TRUE),
    agUse_area_2006 = rowSums(across(c("pasture_area_2006", "crop_area_2006")), na.rm = TRUE)
  )

# EXPORT PREP

# sjlabelled::set_labelS
sjlabelled::set_label(raw_agCensus2006AgUseArea$muni_code) <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw_agCensus2006AgUseArea$pasture_area_2006) <- "pasture area (ha, 2006 Ag Census)"
sjlabelled::set_label(raw_agCensus2006AgUseArea$crop_area_2006) <- "crop area (ha, 2006 Ag Census)"
sjlabelled::set_label(raw_agCensus2006AgUseArea$agUse_area_2006) <- "agricultural use area (ha, 2006 Ag Census)"

# change object name for exportation
clean_agCensus2006AgUseArea <- raw_agCensus2006AgUseArea


# EXPORT
save(clean_agCensus2006AgUseArea,
  file = 
    "data/clean/agcensus2006_ag_usearea.Rdata"
)

# END TIMER
tictoc::toc(log = TRUE)


