
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CREATE AGGREGATED CATEGORIES OF INTEREST BASED ON MAPBIOMAS 30M-PIXELS VALUES
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -


# START TIMER
tictoc::tic(msg = "pixelCategories_prepData.R script", log = TRUE)


# DATA INPUT -----------------------------------------------------------------------------------------------------------------------------------------

# Mapbiomas sample
load("data/prepData/pixelArea_prepData.Rdata")


# DATA MANIPULATION ----------------------------------------------------------------------------------------------------------------------------------

# add aggregated land use/cover categories
pixelCategories_prepData <-
  pixelArea_prepData %>%
  dplyr::mutate(mapbiomas_classAgg = dplyr::case_when(mapbiomas_class == 3 ~ "forest",
                                                      mapbiomas_class == 15 ~ "pasture",
                                                      mapbiomas_class == 39 ~ "soybean",
                                                      mapbiomas_class %in% c(20, 41, 21) ~ "otherCrops", # perennial crops not present in amazon, mosaic only small area in 1985
                                                      mapbiomas_class %in% c(4, 5, 9, 11:13, 22:27, 29:33) ~ "otherCategories"))

# clear environment
rm(pixelArea_prepData)

# EXPORT ONLY PIXELS WITH PRIMARY FOREST BY YEAR
purrr::map(.x = c(2018),
           .f = function(.x) {
             pixelPrimaryForest_prepData <-
               pixelCategories_prepData %>%
               dplyr::filter(year <= .x) %>%
               dplyr::mutate(d_forest = dplyr::if_else(mapbiomas_classAgg == "forest", 1, 0)) %>%
               dplyr::group_by(lon, lat) %>%
               dplyr::mutate(d_primaryForest = dplyr::if_else(sum(d_forest) == (.x-1985+1), 1, 0)) %>%
               dplyr::ungroup() %>%
               dplyr::filter(year == .x, d_primaryForest == 1)

             save(pixelPrimaryForest_prepData,
                  file = file.path("data/prepData",
                                   glue::glue("pixelPrimaryForest{.x}_prepData.Rdata")))
           }

           )


# CLEAN TEMP DIR
terra::tmpFiles(current = TRUE, remove = TRUE)
gc()



# EXPORT PREP ----------------------------------------------------------------------------------------------------------------------------------------

# LABELS
sjlabelled::set_label(pixelCategories_prepData$mapbiomas_classAgg) <- "mapbiomas land use/cover aggregated classification"



# EXPORT ---------------------------------------------------------------------------------------------------------------------------------------------

save(pixelCategories_prepData,
     file = "data/prepData/pixelCategories_prepData.Rdata")


# END TIMER
tictoc::toc(log = TRUE)



# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------
