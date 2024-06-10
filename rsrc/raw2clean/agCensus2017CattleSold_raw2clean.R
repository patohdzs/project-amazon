# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CLEAN RAW CATTLE SOLD - AGRICULTURAL CENSUS 2017 (IBGE)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -


# START TIMER
tictoc::tic(msg = "agCensus2017CattleSold_raw2clean.R script", log = TRUE)

# Read CSV file
cattle_sold_2017 <- readr::read_csv(
  file = "data/raw/ibge/ag_census_2017_cattle_sold/agCensus2017_cattleSold.csv",
  skip = 6,
  na = c("..."),
  col_names = c(
    "muni_code",
    "cattleSoldSmallProp_head_2017",
    "cattleSoldSmallProp_value_2017",
    "cattleSoldSlaughterLargeProp_head_2017",
    "cattleSoldSlaughterLargeProp_value_2017"
  ),
  col_types = "n--cccc"
)

# Remove last rows of the table with notes information
cattle_sold_2017 <- cattle_sold_2017[-c(5550:5563), ]

# Transform "-" values to "0" as explained in the table notes
cattle_sold_2017 <-
  cattle_sold_2017 %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect:::where(is.character),
      function(x) dplyr::if_else(x == "-", "0", x)
    )
  )

# Transform "X" values to "NA"
# (NA identifies values that were omitted to avoid informant identification)
cattle_sold_2017 <-
  cattle_sold_2017 %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect:::where(is.character),
      function(x) dplyr::if_else(x == "X", NA_character_, x)
    )
  )

# Latin character treatment
cattle_sold_2017 <-
  cattle_sold_2017 %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect:::where(is.character),
      \(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
    )
  )

# Transform column class
cattle_sold_2017 <-
  cattle_sold_2017 %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect:::starts_with("cattleSold"),
      function(x) as.numeric(x)
    )
  )

# Calculate total value of cattle sold
cattle_sold_2017 <-
  cattle_sold_2017 %>%
  dplyr::mutate(cattleSold_value_2017 = rowSums(across(c("cattleSoldSmallProp_value_2017", "cattleSoldSlaughterLargeProp_value_2017")), na.rm = TRUE))

# Set labels
sjlabelled::set_label(cattle_sold_2017$muni_code) <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(cattle_sold_2017$cattleSoldSmallProp_head_2017) <- "number of cattle head sold in properties with less than 50 cattle heads (count, 2017 Ag Census)"
sjlabelled::set_label(cattle_sold_2017$cattleSoldSmallProp_value_2017) <- "value of cattle sold in properties with less than 50 cattle heads (thousand BR, 2017 Ag CensusL)"
sjlabelled::set_label(cattle_sold_2017$cattleSoldSlaughterLargeProp_head_2017) <- "number of cattle head sold for slaughter in properties with more than 50 cattle heads (count, 2017 Ag Census)"
sjlabelled::set_label(cattle_sold_2017$cattleSoldSlaughterLargeProp_value_2017) <- "value of cattle sold for slaughter in properties with more than 50 cattle heads (thousand BRL, 2017 Ag Census)"
sjlabelled::set_label(cattle_sold_2017$cattleSold_value_2017) <- "value of cattle sold (thousand BRL, 2017 Ag Census)"

# Save data set
save(cattle_sold_2017, file = "data/clean/cattle_sold_2017.Rdata")

# END TIMER
tictoc::toc(log = TRUE)
