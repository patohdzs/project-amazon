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

# Read CSV files
raw_carbon_price <- readr::read_csv(file = "data/raw/worldbank/carbon_price/carbonPrice_price.CSV", skip = 2, na = "N/A")
raw_carbon_price_revenue <- readr::read_csv(file = "data/raw/worldbank/carbon_price/carbonPrice_revenue.CSV", skip = 2, na = "N/A")


# DATASET CLEANUP AND PREP
raw_carbon_price <-
  raw_carbon_price %>%
  dplyr::left_join(raw_carbon_price_revenue) %>%
  dplyr::select(
    initiative_name = `Name of the initiative`,
    instrument_type = `Instrument Type`,
    revenue_2020 = `...34`,
    price1_2020 = Price_rate_1_2020,
    price2_2020 = Price_rate_2_2020,
    price1_2021 = Price_rate_1_2021,
    price2_2021 = Price_rate_2_2021
  )


# LABELS
sjlabelled::set_label(raw_carbon_price$initiative_name) <- "name of the carbon initiative"
sjlabelled::set_label(raw_carbon_price$instrument_type) <- "type of instrunment (ETS, Carbon Tax or Undecided)"
sjlabelled::set_label(raw_carbon_price$revenue_2020) <- "revenue raised by governments from carbon pricing initiatives in 2020 (USD million)"
sjlabelled::set_label(raw_carbon_price$price1_2020) <- "Nominal carbon prices on February, 01 2020 (USD/tCO2e)"
sjlabelled::set_label(raw_carbon_price$price2_2020) <- "Nominal carbon prices rate 2 on February, 01 2020 (USD/tCO2e)"
sjlabelled::set_label(raw_carbon_price$price1_2021) <- "Nominal carbon prices on April, 01 2021 (USD/tCO2e)"
sjlabelled::set_label(raw_carbon_price$price2_2021) <- "Nominal carbon prices rate 2 on April, 01 2021 (USD/tCO2e)"

# Change object name before saving
carbon_price <- raw_carbon_price

# Save data set
save(carbon_price, file = "data/clean/carbon_price.Rdata")

# END TIMER
tictoc::toc(log = TRUE)
