# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: DOWNLOAD ABOVERGROUND BIOMASS/CARBON DATA (ESA BIOMASS - 2010, 2017, 2018)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: DOWNLOAD PROCESS IS OPTIONAL GIVEN THAT THE DATA IS PROVIDED

# START TIMER
tictoc::tic(msg = "abovegroundBiomassESA_download.R script", log = T)

options(timeout = 300)

# DATA DOWNLOAD

# define download parameters
aux.params <- tidyr::expand_grid(
  year = c(2010, 2017, 2018),
  tile = c(
    "N00W050", "N00W060", "N00W070", "N00W080",
    "N10W050", "N10W060", "N10W070", "N10W080",
    "S10W050", "S10W060", "S10W070", "S10W080"
  )
)

# download AGB layer
purrr::map2(
  .x = aux.params$tile,
  .y = aux.params$year,
  .f = ~ download.file(
    url = glue::glue("https://dap.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/v3.0/geotiff/{.y}/{.x}_ESACCI-BIOMASS-L4-AGB-MERGED-100m-{.y}-fv3.0.tif?download=1"),
    destfile = here::here(glue::glue("data/raw2clean/abovegroundBiomass_esa/input/{.x}_ESACCI-BIOMASS-L4-AGB-MERGED-100m-{.y}-fv3.0.tif")),
    mode = "wb"
  )
)

# download uncertainty layer
purrr::map2(
  .x = aux.params$tile,
  .y = aux.params$year,
  .f = ~ download.file(
    url = glue::glue("https://dap.ceda.ac.uk/neodc/esacci/biomass/data/agb/maps/v3.0/geotiff/{.y}/{.x}_ESACCI-BIOMASS-L4-AGB_SD-MERGED-100m-{.y}-fv3.0.tif?download=1"),
    destfile = here::here(glue::glue("data/raw2clean/abovegroundBiomass_esa/input/{.x}_ESACCI-BIOMASS-L4-AGB_SD-MERGED-100m-{.y}-fv3.0.tif")),
    mode = "wb"
  )
)

# END TIMER
tictoc::toc(log = T)

# # export time to csv table
# ExportTimeProcessing("code/raw2clean")
