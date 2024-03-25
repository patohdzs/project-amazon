
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





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "agCensus2006AgUseArea_raw2clean.R script", log = T)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read csv file
raw.agCensus2006AgUseArea <-  readr::read_csv(file = here::here("data/raw2clean/agCensus2006AgUseArea_ibge/input/agCensus2006_agUseArea.csv"),
                                  skip = 6,
                                  na = c("..."),
                                  col_names = c("muni_code",
                                                "cropPerm_area_2006", "cropTemp_area_2006", "cropCut_area_2006", "cropFlower_area_2006",
                                                "pastureNatural_area_2006", "pasturePlantedGood_area_2006", "pasturePlantedBad_area_2006"),
                                  col_types ="n---cccccc")




# DATA EXPLORATION [disabled for speed]
# summary(raw.agCensus2006AgUseArea)
# View(raw.agCensus2006AgUseArea)





# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# ROW CLEANUP
# remove last rows of the table with notes information
raw.agCensus2006AgUseArea <-  raw.agCensus2006AgUseArea[-(5549:5561), ]

# transform "-" values to "0" as explained in the table notes
raw.agCensus2006AgUseArea <-
  raw.agCensus2006AgUseArea %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# transform "X" values to "NA" (here NA identify which values had to be omitted to avoid informant identification as explained in the table notes)
raw.agCensus2006AgUseArea <-
  raw.agCensus2006AgUseArea %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))

# latin character treatment
raw.agCensus2006AgUseArea <-
  raw.agCensus2006AgUseArea %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), \(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")))

# transform column class
raw.agCensus2006AgUseArea <-
  raw.agCensus2006AgUseArea %>%
  dplyr::mutate(dplyr::across(tidyselect:::ends_with("_area_2006"), function(x) as.numeric(x)))

# sum pasture, crop, and agricultural use area
raw.agCensus2006AgUseArea <-
  raw.agCensus2006AgUseArea %>%
  dplyr::mutate(pasture_area_2006 = rowSums(across(c("pastureNatural_area_2006", "pasturePlantedGood_area_2006", "pasturePlantedBad_area_2006")), na.rm = TRUE),
                crop_area_2006 = rowSums(across(c("cropPerm_area_2006", "cropTemp_area_2006", "cropCut_area_2006", "cropFlower_area_2006")), na.rm = TRUE),
                agUse_area_2006 =  rowSums(across(c("pasture_area_2006", "crop_area_2006")), na.rm = TRUE))



# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# sjlabelled::set_labelS
sjlabelled::set_label(raw.agCensus2006AgUseArea$muni_code)        <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.agCensus2006AgUseArea$pasture_area_2006)    <- "pasture area (ha, 2006 Ag Census)"
sjlabelled::set_label(raw.agCensus2006AgUseArea$crop_area_2006)    <- "crop area (ha, 2006 Ag Census)"
sjlabelled::set_label(raw.agCensus2006AgUseArea$agUse_area_2006)    <- "agricultural use area (ha, 2006 Ag Census)"

# change object name for exportation
clean.agCensus2006AgUseArea <- raw.agCensus2006AgUseArea



# POST-TREATMENT OVERVIEW
# summary(clean.agCensus2006AgUseArea)
# View(clean.agCensus2006AgUseArea)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.agCensus2006AgUseArea,
     file = here::here("data/raw2clean/agCensus2006AgUseArea_ibge/output",
                       "clean_agCensus2006AgUseArea.Rdata"))


# END TIMER
tictoc::toc(log = T)

# # export time to csv table
# ExportTimeProcessing("code/raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------