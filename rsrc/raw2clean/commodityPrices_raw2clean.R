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


# START TIMER
tictoc::tic(msg = "commodityPrices_raw2clean.R script", log = TRUE)

# DATA INPUT

# read xls file
raw_prices <- readxl::read_xls(
  path = "data/raw/seabpr/commodity_prices/ipeadata[22-02-2021-10-04].xls",
  sheet = 1,
  col_names = c(
    "date", "price_rice", "price_cattle", "price_sugarcane",
    "price_cassava", "price_corn", "price_soybean"
  ),
  col_types = c("text", rep("numeric", 6)),
  skip = 1
)


# DATASET CLEANUP AND PREP

# transform to date class
raw_prices <-
  raw_prices %>%
  mutate(date = lubridate::ymd(date, truncated = 1))

# EXPORT PREP

# LABELS
sjlabelled::set_label(raw_prices$date) <- "calendar date (yyyy-mm-dd), monthly data, all 'dd' set to 01"
sjlabelled::set_label(raw_prices$price_cattle) <- "price (nominal), average received by producer - cattle (1@; PR)"
sjlabelled::set_label(raw_prices$price_cassava) <- "price (nominal), average received by producer - cassava (1t; PR)"
sjlabelled::set_label(raw_prices$price_corn) <- "price (nominal), average received by producer - corn (60kg; PR)"
sjlabelled::set_label(raw_prices$price_rice) <- "price (nominal), average received by producer - rice (50kg; PR)"
sjlabelled::set_label(raw_prices$price_soybean) <- "price (nominal), average received by producer - soybean (60kg; PR)"
sjlabelled::set_label(raw_prices$price_sugarcane) <- "price (nominal), average received by producer - sugarcane (1t; PR)"

# change object name for exportation
clean_commodityPrices <- raw_prices

# EXPORT

save(clean_commodityPrices,
  file = 
    "data/clean/commodity_prices.Rdata"
)

# END TIMER
tictoc::toc(log = TRUE)


