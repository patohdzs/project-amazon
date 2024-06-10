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
tictoc::tic(msg = "deflatorIPA_raw2clean.R script", log = TRUE)

# Read CSV file
raw_deflator <- readr::read_csv2(
  file = "data/raw/fgv/deflator_ipa/ipeadata[23-02-2021-10-02].csv",
  col_names = c("date", "deflator_ipa"),
  col_types = c("c", "n"),
  skip = 1
)

# Transform to date class
raw_deflator <-
  raw_deflator %>%
  dplyr::mutate(date = lubridate::ymd(date, truncated = 1))


# Set labels
sjlabelled::set_label(raw_deflator$date) <- "calendar date (yyyy-mm-dd), monthly data, all 'dd' set to 01"
sjlabelled::set_label(raw_deflator$deflator_ipa) <- "deflator (IPA-DI; 1994-08-01 = 100)"

# Change object name before saving
clean_deflator <- raw_deflator

# Save data set
save(clean_deflator,
  file =
    "data/clean/deflator.Rdata"
)

# END TIMER
tictoc::toc(log = TRUE)
