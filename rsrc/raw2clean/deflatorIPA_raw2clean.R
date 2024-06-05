# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CLEAN RAW DEFLATOR IPA-DI - FGV
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

# START TIMER
tictoc::tic(msg = "deflatorIPA_raw2clean.R script", log = T)

# DATA INPUT

# read csv file
raw.deflator <- readr::read_csv2(
  file = here::here("data/raw2clean/deflatorIPA_fgv/input/ipeadata[23-02-2021-10-02].csv"),
  col_names = c("date", "deflator_ipa"),
  col_types = c("c", "n"),
  skip = 1
)

# DATA EXPLORATION [disabled for speed]
# summary(raw.deflator)
# View(raw.deflator)

# DATASET CLEANUP AND PREP

# transform to date class
raw.deflator <-
  raw.deflator %>%
  dplyr::mutate(date = lubridate::ymd(date, truncated = 1))

# EXPORT PREP

# sjlabelled::set_labelS
sjlabelled::set_label(raw.deflator$date) <- "calendar date (yyyy-mm-dd), monthly data, all 'dd' set to 01"
sjlabelled::set_label(raw.deflator$deflator_ipa) <- "deflator (IPA-DI; 1994-08-01 = 100)"

# change object name for exportation
clean.deflatorIPA <- raw.deflator

# POST-TREATMENT OVERVIEW
# summary(clean.deflatorIPA)
# View(clean.deflatorIPA)

# EXPORT

save(clean.deflatorIPA,
  file = here::here(
    "data/raw2clean/deflatorIPA_fgv/output",
    "clean_deflatorIPA.Rdata"
  )
)

# END TIMER
tictoc::toc(log = T)

# export time to csv table
# ExportTimeProcessing("code/raw2clean")
