# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CLEAN RAW AGRICULTURAL COMMODITY PRICES FROM THE PARANA SECRETARIAT OF SUPPLY AND AGRICULTURE (SEAB-PR)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# -
library(sf)
library(tidyverse)
library(tictoc)
library(sjlabelled)
library(conflicted)
library(readxl)

# Resolve conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

# START TIMER
tic(msg = "commodityPrices_raw2clean.R script", log = TRUE)

# Read XLS file
raw_commodity_prices <- read_xls(
  path = "data/raw/seabpr/commodity_prices/ipeadata[22-02-2021-10-04].xls",
  sheet = 1,
  col_names = c(
    "date",
    "price_rice",
    "price_cattle",
    "price_sugarcane",
    "price_cassava",
    "price_corn",
    "price_soybean"
  ),
  col_types = c("text", rep("numeric", 6)),
  skip = 1
)

# Transform to date class
raw_commodity_prices <-
  raw_commodity_prices %>%
  mutate(date = ymd(date, truncated = 1))

# Set labels
set_label(raw_commodity_prices$date) <- "calendar date (yyyy-mm-dd), monthly data, all 'dd' set to 01"
set_label(raw_commodity_prices$price_cattle) <- "price (nominal), average received by producer - cattle (1@; PR)"
set_label(raw_commodity_prices$price_cassava) <- "price (nominal), average received by producer - cassava (1t; PR)"
set_label(raw_commodity_prices$price_corn) <- "price (nominal), average received by producer - corn (60kg; PR)"
set_label(raw_commodity_prices$price_rice) <- "price (nominal), average received by producer - rice (50kg; PR)"
set_label(raw_commodity_prices$price_soybean) <- "price (nominal), average received by producer - soybean (60kg; PR)"
set_label(raw_commodity_prices$price_sugarcane) <- "price (nominal), average received by producer - sugarcane (1t; PR)"

# Change object name before saving
commodity_prices <- raw_commodity_prices

# Save data set
save(commodity_prices, file = "data/clean/commodity_prices.Rdata")

# END TIMER
toc(log = TRUE)
