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


# START TIMER
tictoc::tic(msg = "agCensus2017AgUseArea_raw2clean.R script", log = TRUE)

# Read csv file
use_area_2017 <- readr::read_csv(
  file = "data/raw/ibge/ag_census_2017_ag_use_area/agCensus2017_agUseArea.csv",
  skip = 7,
  na = c("..."),
  col_names = c(
    "muni_code",
    "cropPerm_area_2017",
    "cropTemp_area_2017",
    "cropFlower_area_2017",
    "pastureNatural_area_2017",
    "pasturePlantedGood_area_2017",
    "pasturePlantedBad_area_2017"
  ),
  col_types = "n--cccccc"
)


# Remove last rows of the table with notes information
use_area_2017 <- use_area_2017[-(5564:5577), ]

# Transform "-" values to "0" as explained in the table notes
use_area_2017 <-
  use_area_2017 %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect:::where(is.character),
      function(x) dplyr::if_else(x == "-", "0", x)
    )
  )

# Transform "X" values to "NA"
# (NA's identify values that had to be omitted to avoid informant identification)
use_area_2017 <-
  use_area_2017 %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect:::where(is.character),
      function(x) dplyr::if_else(x == "X", NA_character_, x)
    )
  )

# Latin character treatment
use_area_2017 <-
  use_area_2017 %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect:::where(is.character),
      \(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
    )
  )

# Transform column class
use_area_2017 <-
  use_area_2017 %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect:::ends_with("_area_2017"),
      function(x) as.numeric(x)
    )
  )

# Sum pasture, crop, and agricultural use area
use_area_2017 <-
  use_area_2017 %>%
  dplyr::mutate(
    pasture_area_2017 = rowSums(across(c("pastureNatural_area_2017", "pasturePlantedGood_area_2017", "pasturePlantedBad_area_2017")), na.rm = TRUE),
    crop_area_2017 = rowSums(across(c("cropPerm_area_2017", "cropTemp_area_2017", "cropFlower_area_2017")), na.rm = TRUE),
    agUse_area_2017 = rowSums(across(c("pasture_area_2017", "crop_area_2017")), na.rm = TRUE)
  )

# Set labels
sjlabelled::set_label(use_area_2017$muni_code) <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(use_area_2017$pasture_area_2017) <- "pasture area (ha, 2017 Ag Census)"
sjlabelled::set_label(use_area_2017$crop_area_2017) <- "crop area (ha, 2017 Ag Census)"
sjlabelled::set_label(use_area_2017$agUse_area_2017) <- "agricultural use area (ha, 2017 Ag Census)"

# Save data set
save(use_area_2017, file = "data/clean/use_area_2017.Rdata")

# END TIMER
tictoc::toc(log = TRUE)
