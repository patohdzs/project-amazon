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

library(sf)
library(tidyverse)
library(tictoc)
library(sjlabelled)
library(conflicted)

# Resolve conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

# START TIMER
tic(msg = "agCensus2006AgUseArea_raw2clean.R script", log = TRUE)

# Read csv file
use_area_2006 <- read_csv(
  file = "data/raw/ibge/ag_census_2006_ag_use_area/agCensus2006_agUseArea.csv",
  skip = 6,
  na = c("..."),
  col_names = c(
    "muni_code",
    "cropPerm_area_2006",
    "cropTemp_area_2006",
    "cropCut_area_2006",
    "cropFlower_area_2006",
    "pastureNatural_area_2006",
    "pasturePlantedGood_area_2006",
    "pasturePlantedBad_area_2006"
  ),
  col_types = "n---cccccc"
)

# Remove last rows of the table with notes information
use_area_2006 <- use_area_2006[-(5549:5561), ]

# Transform "-" values to "0" as explained in the table notes
use_area_2006 <-
  use_area_2006 %>%
  mutate(
    across(
      where(is.character),
      function(x) if_else(x == "-", "0", x)
    )
  )

# Transform "X" values to "NA"
# (NA's identify values that were omitted to avoid informant identification)
use_area_2006 <-
  use_area_2006 %>%
  mutate(
    across(
      where(is.character),
      function(x) if_else(x == "X", NA_character_, x)
    )
  )

# Latin character treatment
use_area_2006 <-
  use_area_2006 %>%
  mutate(
    across(
      where(is.character),
      \(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
    )
  )

# Transform column class
use_area_2006 <-
  use_area_2006 %>%
  mutate(
    across(
      ends_with("_area_2006"),
      function(x) as.numeric(x)
    )
  )

# Sum pasture, crop, and agricultural use area
use_area_2006 <-
  use_area_2006 %>%
  mutate(
    pasture_area_2006 = rowSums(across(c("pastureNatural_area_2006", "pasturePlantedGood_area_2006", "pasturePlantedBad_area_2006")), na.rm = TRUE),
    crop_area_2006 = rowSums(across(c("cropPerm_area_2006", "cropTemp_area_2006", "cropCut_area_2006", "cropFlower_area_2006")), na.rm = TRUE),
    agUse_area_2006 = rowSums(across(c("pasture_area_2006", "crop_area_2006")), na.rm = TRUE)
  )

# Set labels
set_label(use_area_2006$muni_code) <- "municipality code (7-digit, IBGE)"
set_label(use_area_2006$pasture_area_2006) <- "pasture area (ha, 2006 Ag Census)"
set_label(use_area_2006$crop_area_2006) <- "crop area (ha, 2006 Ag Census)"
set_label(use_area_2006$agUse_area_2006) <- "agricultural use area (ha, 2006 Ag Census)"

# Save data set
save(use_area_2006, file = "data/clean/use_area_2006.Rdata")

# END TIMER
toc(log = TRUE)
