# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CLEAN RAW SEEG AGRICULTURAL EMISSION DATA - LEGAL AMAZON STATES 1990-2019
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# -

# START TIMER
tictoc::tic(msg = "emission_raw2clean.R script", log = T)

# DATA INPUT

# read csv file
raw.emission <- readxl::read_xlsx(
  path = here::here("data/raw2clean/emission_seeg/input/emission_agriculture_states.xlsx"),
  sheet = 1, col_names = c("state_uf", 1990:2019), skip = 1
)

raw.removal <- readxl::read_xlsx(
  path = here::here("data/raw2clean/emission_seeg/input/removalNCI_agriculture_states.xlsx"),
  sheet = 1, col_names = c("state_uf", 1990:2019), skip = 1
)

# DATA EXPLORATION
# summary(raw.emission)    # object is a list of data frames
# class(raw.emission)
# View(raw.emission)      # column names indicate file of origin

# DATASET CLEANUP AND PREP

# RESHAPE
raw.emission <-
  raw.emission %>%
  tidyr::pivot_longer(-state_uf, names_to = "year", values_to = "emission_co2e")

raw.removal <-
  raw.removal %>%
  tidyr::pivot_longer(-state_uf, names_to = "year", values_to = "removal_co2e")

# MERGE
raw.emission <-
  raw.emission %>%
  dplyr::left_join(raw.removal)

# ADD NET EMISSION VARIABLE
raw.emission <-
  raw.emission %>%
  dplyr::mutate(netEmission_co2e = emission_co2e + removal_co2e) %>%
  dplyr::mutate(year = as.numeric(year))

# clean environmnet
rm(raw.removal)

# EXPORT PREP

# sjlabelled::set_labelS
sjlabelled::set_label(raw.emission$state_uf) <- "state name abbreviation"
sjlabelled::set_label(raw.emission$year) <- "year of reference (calendar or PRODES year)"
sjlabelled::set_label(raw.emission$emission_co2e) <- "total emissions from agricultural land (CO2e-GWP-AR5)"
sjlabelled::set_label(raw.emission$removal_co2e) <- "total removals from agricultural land (CO2e-GWP-AR5)"
sjlabelled::set_label(raw.emission$netEmission_co2e) <- "total net emissions from agricultural land (CO2e-GWP-AR5)"

# change object name for exportation
clean.emission <- raw.emission

# POST-TREATMENT OVERVIEW
# summary(clean.emission)
# View(clean.emission)

# EXPORT

save(clean.emission,
  file = here::here(
    "data/raw2clean/emission_seeg/output",
    "clean_emission.Rdata"
  )
)

# END TIMER
tictoc::toc(log = T)

# export time to csv table
# ExportTimeProcessing("code/raw2clean")
