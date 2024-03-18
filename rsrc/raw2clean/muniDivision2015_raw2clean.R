
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CLEAN RAW DATA SHAPEFILE OF MUNI DIVISION (2015)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "muniDivision2015_raw2clean.R script", log = T)





# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# read shapefile
raw.muni <- sf::st_read(dsn   = here::here("data/raw2clean/muniDivision2015_ibge/input"),
                        layer = "BRMUE250GC_SIR")



# DATA EXPLORATION [disabled for speed]
# summary(raw.muni)
# View(raw.muni)
# plot(raw.muni$geometry)




# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# COLUMN CLEANUP
# names
colnames(raw.muni)

# translate column names
raw.muni <-
  raw.muni %>%
  dplyr::rename(muni_code = CD_GEOCMU,
                muni_name = NM_MUNICIP)

# add state_uf column
raw.muni <-
  raw.muni %>%
  dplyr::mutate(state_uf = case_when(str_sub(muni_code, 1, 2) == 11 ~  "RO",
                                     str_sub(muni_code, 1, 2) == 12 ~  "AC",
                                     str_sub(muni_code, 1, 2) == 13 ~  "AM",
                                     str_sub(muni_code, 1, 2) == 14 ~  "RR",
                                     str_sub(muni_code, 1, 2) == 15 ~  "PA",
                                     str_sub(muni_code, 1, 2) == 16 ~  "AP",
                                     str_sub(muni_code, 1, 2) == 17 ~  "TO",
                                     str_sub(muni_code, 1, 2) == 21 ~  "MA",
                                     str_sub(muni_code, 1, 2) == 22 ~  "PI",
                                     str_sub(muni_code, 1, 2) == 23 ~  "CE",
                                     str_sub(muni_code, 1, 2) == 24 ~  "RN",
                                     str_sub(muni_code, 1, 2) == 25 ~  "PB",
                                     str_sub(muni_code, 1, 2) == 26 ~  "PE",
                                     str_sub(muni_code, 1, 2) == 27 ~  "AL",
                                     str_sub(muni_code, 1, 2) == 28 ~  "SE",
                                     str_sub(muni_code, 1, 2) == 29 ~  "BA",
                                     str_sub(muni_code, 1, 2) == 31 ~  "MG",
                                     str_sub(muni_code, 1, 2) == 32 ~  "ES",
                                     str_sub(muni_code, 1, 2) == 33 ~  "RJ",
                                     str_sub(muni_code, 1, 2) == 35 ~  "SP",
                                     str_sub(muni_code, 1, 2) == 41 ~  "PR",
                                     str_sub(muni_code, 1, 2) == 42 ~  "SC",
                                     str_sub(muni_code, 1, 2) == 43 ~  "RS",
                                     str_sub(muni_code, 1, 2) == 50 ~  "MS",
                                     str_sub(muni_code, 1, 2) == 51 ~  "MT",
                                     str_sub(muni_code, 1, 2) == 52 ~  "GO",
                                     str_sub(muni_code, 1, 2) == 53 ~  "DF"))

# class - muni_code should be numeric
lapply(raw.muni, class)

raw.muni <- raw.muni %>% dplyr::mutate(muni_code = as.numeric(muni_code))



# LATIN CHARACTER TREATMENT
raw.muni <-
  raw.muni %>%
  dplyr::mutate(dplyr::across(tidyselect:::where(is.character), \(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")))



# LETTERS CAPITALIZATION
raw.muni <-
  raw.muni %>%
  dplyr::mutate(muni_name = toupper(muni_name))



# PROJECTION
raw.muni <- sf::st_transform(x = raw.muni, crs = 5880) # SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)



file_path <- paste(getwd(), "data/raw2clean/muniDivision2015_ibge/output", paste0("raw_muni", ".Rdata"), sep = "/")

save(raw.muni, file = file_path)

# save(raw.muni,
#      file = here::here("data/calibration/",
#                       "raw_muni.Rdata"))





# REMOVE POLYGONS IDENTIFIED AS BODY OF WATERS AND NOT MUNICIPALITIES - see muniDivision2007_raw2clean
raw.muni <- raw.muni %>% dplyr::filter(!muni_code %in% c(4300001, 4300002))




# GEOMETRY CLEANUP
raw.muni <- sf::st_make_valid(raw.muni)





# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(raw.muni$muni_code)            <- "municipality code (7-digit, IBGE)"
sjlabelled::set_label(raw.muni$muni_name)            <- "municipality name"
sjlabelled::set_label(raw.muni$state_uf)             <- "state name (abbreviation)"




# change object name for exportation
clean.muniDivision2015 <- raw.muni



# POST-TREATMENT OVERVIEW
# summary(clean.muniDivision2015)
# View(clean.muniDivision2015)
# plot(clean.muniDivision2015$geoemtry)





# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(clean.muniDivision2015,
     file = here::here("data/raw2clean/muniDivision2015_ibge/output",
                      "clean_muniDivision2015.Rdata"))


# END TIMER
tictoc::toc(log = T)

# # export time to csv table
# ExportTimeProcessing("code/raw2clean")





# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------