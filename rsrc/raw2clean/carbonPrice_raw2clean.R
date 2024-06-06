# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: TREAT RAW DATA CARBON PRICE (WORLD BANK)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -


# START TIMER
tictoc::tic(msg = "carbonPrice_raw2clean.R script", log = TRUE)

# DATA INPUT

# read input file
raw_carbonPrice <- readr::read_csv(file = "data/raw/worldbank/carbon_price/carbonPrice_price.csv", skip = 2, na = "N/A")
raw_carbonPriceRevenue <- readr::read_csv(file = "data/raw/worldbank/carbon_price/carbonPrice_revenue.csv", skip = 2, na = "N/A")


# DATASET CLEANUP AND PREP

raw_carbonPrice <-
  raw_carbonPrice %>%
  dplyr::left_join(raw_carbonPriceRevenue) %>%
  dplyr::select(
    initiative_name = `Name of the initiative`, instrument_type = `Instrument Type`,
    revenue_2020 = `...34`, price1_2020 = Price_rate_1_2020, price2_2020 = Price_rate_2_2020,
    price1_2021 = Price_rate_1_2021, price2_2021 = Price_rate_2_2021
  )

# EXPORT PREP

# LABELS
sjlabelled::set_label(raw_carbonPrice$initiative_name) <- "name of the carbon initiative"
sjlabelled::set_label(raw_carbonPrice$instrument_type) <- "type of instrunment (ETS, Carbon Tax or Undecided)"
sjlabelled::set_label(raw_carbonPrice$revenue_2020) <- "revenue raised by governments from carbon pricing initiatives in 2020 (USD million)"
sjlabelled::set_label(raw_carbonPrice$price1_2020) <- "Nominal carbon prices on February, 01 2020 (USD/tCO2e)"
sjlabelled::set_label(raw_carbonPrice$price2_2020) <- "Nominal carbon prices rate 2 on February, 01 2020 (USD/tCO2e)"
sjlabelled::set_label(raw_carbonPrice$price1_2021) <- "Nominal carbon prices on April, 01 2021 (USD/tCO2e)"
sjlabelled::set_label(raw_carbonPrice$price2_2021) <- "Nominal carbon prices rate 2 on April, 01 2021 (USD/tCO2e)"

# change object name for exportation
clean_carbonPrice <- raw_carbonPrice

# EXPORT

save(clean_carbonPrice,
  file = 
    "data/clean/carbon_price.Rdata"
)

# END TIMER
tictoc::toc(log = TRUE)


