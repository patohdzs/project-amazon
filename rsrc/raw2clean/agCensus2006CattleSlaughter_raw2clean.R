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
tictoc::tic(msg = "agCensus2006CattleSlaughter_raw2clean.R script", log = TRUE)

# Read csv file
in_path <- "data/raw/ibge/ag_census_2006_cattle_slaughter/agCensus2006_cattleSlaughter.csv"

cattle_slaughter_2006 <- readr::read_csv(
  file = in_path,
  skip = 9,
  na = c("..."),
  col_names = c(
    "muni_code",
    "cattleSlaughter_type",
    "cattleSlaughter_value",
    "cattleSlaughter_head"
  ),
  col_types = "ncccc"
)


# Remove last rows of the table with notes information
cattle_slaughter_2006 <- cattle_slaughter_2006[-c(16030:16041), ]

# Transform "-" values to "0" as explained in the table notes
cattle_slaughter_2006 <-
  cattle_slaughter_2006 %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect:::where(is.character),
      function(x) dplyr::if_else(x == "-", "0", x)
    )
  )

# Transform "X" values to "NA"
# (NA's identify values that were omitted to avoid informant identification)
cattle_slaughter_2006 <-
  cattle_slaughter_2006 %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect:::where(is.character),
      function(x) dplyr::if_else(x == "X", NA_character_, x)
    )
  )

# Transform column class
cattle_slaughter_2006 <-
  cattle_slaughter_2006 %>%
  dplyr::mutate(
    dplyr::across(
      c("cattleSlaughter_value", "cattleSlaughter_head"),
      function(x) as.numeric(x)
    )
  )

# Sum all cattle sold for slaughter valus
cattle_slaughter_2006 <-
  cattle_slaughter_2006 %>%
  dplyr::group_by(muni_code) %>%
  dplyr::summarise(
    cattleSlaughter2006_value = sum(cattleSlaughter_value, na.rm = TRUE),
    cattleSlaughter2006_head = sum(cattleSlaughter_head, na.rm = TRUE)
  ) %>%
  dplyr::select(
    muni_code,
    cattleSlaughter2006_value,
    cattleSlaughter2006_head
  )

# Set labels
sjlabelled::set_label(cattle_slaughter_2006$muni_code) <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(cattle_slaughter_2006$cattleSlaughter2006_value) <- "value of cattle sold for slaughter in 2006 (thousand BRL)"

# Save data set
save(cattle_slaughter_2006, file = "data/clean/cattle_slaughter_2006.Rdata")

# END TIMER
tictoc::toc(log = TRUE)
