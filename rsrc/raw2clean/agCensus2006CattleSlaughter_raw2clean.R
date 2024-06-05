# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CLEAN RAW CATTLE FOR SLAUGHTER - AGRICULTURAL CENSUS 2006 (IBGE)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

# START TIMER
tictoc::tic(msg = "agCensus2006CattleSlaughter_raw2clean.R script", log = T)

# DATA INPUT

# read csv file
raw.agCensus2006CattleSlaughter <- readr::read_csv(
  file = "data/raw2clean/agCensus2006CattleSlaughter_ibge/input/agCensus2006_cattleSlaughter.csv",
  skip = 9,
  na = c("..."),
  col_names = c("muni_code", "cattleSlaughter_type", "cattleSlaughter_value", "cattleSlaughter_head"),
  col_types = "ncccc"
)

# DATA EXPLORATION [disabled for speed]
# summary(raw.agCensus2006CattleSlaughter)
# View(raw.agCensus2006CattleSlaughter)

# DATASET CLEANUP AND PREP

# ROW CLEANUP
# remove last rows of the table with notes information
raw.agCensus2006CattleSlaughter <- raw.agCensus2006CattleSlaughter[-c(16030:16041), ]

# transform "-" values to "0" as explained in the table notes
raw.agCensus2006CattleSlaughter <-
  raw.agCensus2006CattleSlaughter %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "-", "0", x)))

# transform "X" values to "NA" (here NA identify which values had to be omitted to avoid informant identification as explained in the table notes)
raw.agCensus2006CattleSlaughter <-
  raw.agCensus2006CattleSlaughter %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), function(x) dplyr::if_else(x == "X", NA_character_, x)))

# transform column class
raw.agCensus2006CattleSlaughter <-
  raw.agCensus2006CattleSlaughter %>%
  dplyr::mutate(dplyr::across(c("cattleSlaughter_value", "cattleSlaughter_head"), function(x) as.numeric(x)))

# sum all cattle sold for slaughter valus
raw.agCensus2006CattleSlaughter <-
  raw.agCensus2006CattleSlaughter %>%
  dplyr::group_by(muni_code) %>%
  dplyr::summarise(
    cattleSlaughter2006_value = sum(cattleSlaughter_value, na.rm = T),
    cattleSlaughter2006_head = sum(cattleSlaughter_head, na.rm = T)
  ) %>%
  dplyr::select(muni_code, cattleSlaughter2006_value, cattleSlaughter2006_head)

# EXPORT PREP

# sjlabelled::set_labelS
sjlabelled::set_label(raw.agCensus2006CattleSlaughter$muni_code) <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.agCensus2006CattleSlaughter$cattleSlaughter2006_value) <- "value of cattle sold for slaughter in 2006 (thousand BRL)"

# change object name for exportation
clean.agCensus2006CattleSlaughter <- raw.agCensus2006CattleSlaughter

# POST-TREATMENT OVERVIEW
# summary(clean.agCensus2006CattleSlaughter)
# View(clean.agCensus2006CattleSlaughter)

# EXPORT

save(clean.agCensus2006CattleSlaughter,
  file = file.path(
    "data/raw2clean/agCensus2006CattleSlaughter_ibge/output",
    "clean_agCensus2006CattleSlaughter.Rdata"
  )
)

# END TIMER
tictoc::toc(log = T)

# # export time to csv table
# ExportTimeProcessing("code/raw2clean")
