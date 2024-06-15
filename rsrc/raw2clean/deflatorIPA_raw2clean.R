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
tic(msg = "deflatorIPA_raw2clean.R script", log = TRUE)

# Read CSV file
deflator <- read_csv2(
  file = "data/raw/fgv/deflator_ipa/ipeadata[23-02-2021-10-02].csv",
  col_names = c("date", "deflator_ipa"),
  col_types = c("c", "n"),
  skip = 1
)

# Transform to date class
deflator <-
  deflator %>%
  mutate(date = ymd(date, truncated = 1))


# Set labels
set_label(deflator$date) <- "calendar date (yyyy-mm-dd), monthly data, all 'dd' set to 01"
set_label(deflator$deflator_ipa) <- "deflator (IPA-DI; 1994-08-01 = 100)"


# Save data set
save(deflator, file = "data/clean/deflator.Rdata")

# END TIMER
toc(log = TRUE)
