
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: MAPS OF 78 SITES MODEL PREDICTED VALUES
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -





# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "mapsPrediction_78SitesModel.R script", log = T)


# OPTIONS
options(scipen = 999)





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# # GLOBAL MODEL CALIBRATION VARIABLES
# load(here::here("data/calibration/globalModel", "calibration_globalModel.Rdata"))

# extract high and low price
p_high <- 42.03


# # clear unnecessary objects
# rm(matrixTransition.2prices, calibration.globalModel)


# 78 SITES MODEL CALIBRATION VARIABLES
load(here::here("data/hmc", "hmc_78SitesModel.Rdata"))




# 78 SITES AGGREGATE PREDICTION
aux.prices <- c(6.5, 16.5, 21.5, 26.5)



# read the site-level predictions with high cattle price to calculate aggregated agricultural output
prediction.78SitesModel <- purrr::map_df(.x = aux.prices,
                                           .f = ~readr::read_delim(here::here("data/prediction/78-hmc",
                                                                              glue::glue("p_e_{.x}/amazon_data_z.dat")),
                                                                   col_types = cols(.default = "n")) %>%
                                             dplyr::mutate(p_e = .x, p_a = p_high)) 





# # read the site-level predictions with low cattle price to calculate aggregated agricultural output
# prediction.78SitesModel <- purrr::map_df(.x = aux.prices,
#                                            .f = ~readr::read_delim(here::here("data/prediction/78SitesModel",
#                                                                               glue::glue("pl_p_e_{.x}/amazon_data_z.dat")),
#                                                                    col_types = cols(.default = "n")) %>%
#                                              dplyr::mutate(p_e = .x, p_a = p_low)) %>%
#   bind_rows(prediction.78SitesModel)


# AMAZON BIOME VECTOR DATA
load(here::here("data/raw2clean/amazonBiome_ibge/output/clean_amazonBiome.Rdata"))




# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# select relevant calibrated variables
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::mutate(rank_theta_78Sites = dense_rank(desc(theta_78Sites)), # rank values
                rank_gamma_78Sites = dense_rank(desc(gamma_78Sites))) %>%  # rank values
  dplyr::select(id, z_2017_78Sites, theta_78Sites,  rank_theta_78Sites, rank_gamma_78Sites, x_2017_78Sites, zbar_2017_78Sites)

# adjust predicted data to site-year panel + add calibrated variables
prediction.78SitesModel <-
  prediction.78SitesModel %>%
  tidyr::pivot_longer(cols = -c("T/R ", "p_a", "p_e"), names_to = "id", values_to = "z_t") %>%
  dplyr::mutate(id = as.numeric(stringr::str_trim(id))) %>%
  dplyr::rename(time = "T/R ") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
  dplyr::right_join(calibration.78SitesModel, by = c("id" = "id")) %>% # match by id, guarantee that prediction data uses the same id than the calibrated data
  sf::st_as_sf() %>%
  dplyr::mutate(z_t = 100*1e11*z_t/zbar_2017_78Sites,
                z_2017_78Sites = 100*z_2017_78Sites/zbar_2017_78Sites) # transform share to %


# clean environment
rm(calibration.78SitesModel)

# adjust projection
prediction.78SitesModel <- sf::st_as_sf(prediction.78SitesModel)
prediction.78SitesModel <- sf::st_transform(prediction.78SitesModel, sf::st_crs(4326))
clean.amazonBiome <- sf::st_transform(clean.amazonBiome, sf::st_crs(prediction.78SitesModel))





# GENERATE MAPS --------------------------------------------------------------------------------------------------------------------------------------

# PLOT z_2017, GAMMA, NAD THETA CALIBRATED VALUES

# z_2017_78Sites
z_2017_78Sites <-
    ggplot2::ggplot(data = prediction.78SitesModel %>%
                      dplyr::filter(time == 0, p_e == 21.5) %>%
                      dplyr::mutate(z_t = cut(round(z_t,5),
                                              breaks = c(0,0.5,20,40,60,80,105),
                                              include.lowest = T,
                                              dig.lab = 3,
                                              labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")))) +
    ggplot2::geom_sf(aes(fill = z_t)) +
    ggplot2::scale_fill_manual(name = expr(paste("Z"[!!2017]^"i", ~"(%)")), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
    ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
    ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
    ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                   panel.grid.minor = element_line(colour = "white"),
                   panel.background = element_blank(),
                   strip.background = element_rect(fill = NA),
                   axis.line = element_blank(), axis.ticks = element_blank(),
                   axis.title = element_blank(), axis.text = element_blank(),
                   legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
                   legend.position = "bottom", legend.margin=margin(t=-1, r=-0, b=0.3, l=-0, unit="cm"),
                   legend.text = element_text(size = 28, face = "bold"),
                   plot.margin=unit(c(t=-0.5,r=-1,b=-0,l=-1), "cm"))


# gamma_78Sites
gamma_78Sites <-
  ggplot2::ggplot(data = prediction.78SitesModel %>% dplyr::filter(time == 0, p_e == 6.5) %>%
                    dplyr::mutate(rank_gamma_78Sites = cut(round(rank_gamma_78Sites), breaks = c(1, 17, 33, 49, 65, 78), include.lowest = T, dig.lab = 4))) +
  ggplot2::geom_sf(aes(fill = rank_gamma_78Sites)) +
  ggplot2::scale_fill_brewer(name = expression(paste(gamma^"i", ~"(rank)")), palette = "YlOrRd", direction = -1) +
  ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                 panel.grid.minor = element_line(colour = "white"),
                 panel.background = element_blank(),
                 strip.background = element_rect(fill = NA),
                 axis.line = element_blank(), axis.ticks = element_blank(),
                 axis.title = element_blank(), axis.text = element_blank(),
                 legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
                 legend.position = "bottom", legend.margin=margin(t=-1, r=-0, b=0.3, l=-0, unit="cm"),
                 legend.text = element_text(size = 28, face = "bold"),
                 plot.margin=unit(c(t=-0.5,r=-1,b=-0,l=-1), "cm"))

# theta_78Sites
theta_78Sites <-
  ggplot2::ggplot(data = prediction.78SitesModel %>% dplyr::filter(time == 0, p_e == 6.5) %>%
                    dplyr::mutate(rank_theta_78Sites = cut(round(rank_theta_78Sites), breaks = c(1, 17, 33, 49, 65, 78), include.lowest = T, dig.lab = 4))) +
  ggplot2::geom_sf(aes(fill = rank_theta_78Sites)) +
  ggplot2::scale_fill_brewer(name = expression(paste(theta^"i", ~"(rank)")), palette = "YlOrRd", direction = -1) +
  ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                 panel.grid.minor = element_line(colour = "white"),
                 panel.background = element_blank(),
                 strip.background = element_rect(fill = NA),
                 axis.line = element_blank(), axis.ticks = element_blank(),
                 axis.title = element_blank(), axis.text = element_blank(),
                 legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
                 legend.position = "bottom", legend.margin=margin(t=-1, r=-0, b=0.3, l=-0, unit="cm"),
                 legend.text = element_text(size = 28, face = "bold"),
                 plot.margin=unit(c(t=-0.5,r=-1,b=-0,l=-1), "cm"))


# LOOP ACROSS PRICES TO GENERATE Z50 MAP AND SAVE Z0, Z50, GAMMA, AND THETA MAPS BY PRICE
mapList <- list()
mapIndex <- 1
for (price in aux.prices) {
  transfer <- price-6.5
  # z50 map (vary by model)
  mapList[[mapIndex]] <-
    ggplot2::ggplot(data = prediction.78SitesModel %>%
                            dplyr::filter(time == 50, p_e == price) %>%
                            dplyr::mutate(z_t = cut(round(z_t,5),
                                                    breaks = c(0,0.5,20,40,60,80,105),
                                                    include.lowest = T,
                                                    dig.lab = 3,
                                                    labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")))) +
    ggplot2::geom_sf(aes(fill = z_t)) +
    ggplot2::scale_fill_manual(name = expr(paste("Z"[2067]^"i", ~"(%), ", "b", "=", !!transfer)), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
    ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
    ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
    ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                   panel.grid.minor = element_line(colour = "white"),
                   panel.background = element_blank(),
                   strip.background = element_rect(fill = NA),
                   axis.line = element_blank(), axis.ticks = element_blank(),
                   axis.title = element_blank(), axis.text = element_blank(),
                   legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
                   legend.position = "bottom", legend.margin=margin(t=-1, r=-0, b=0.3, l=-0, unit="cm"),
                   legend.text = element_text(size = 28, face = "bold"),
                   plot.margin=unit(c(t=-0.5,r=-1,b=-0,l=-1), "cm"))

    # SAVE MAP Z0, Z50, GAMMA, THETA
    ggpubr::ggarrange(mapList[[mapIndex]]) %>%
      ggpubr::ggexport(filename = here::here(glue::glue("results/prediction/78-hmc/map_z50_78Sites_pe{price}.png")),
                       width = 2400, height = 1500)

    mapIndex <- mapIndex + 1

}

# FIGURE WTIH ALL MAPS
# SAVE MAP Z0, GAMMA, THETA, with varying prices
ggpubr::ggarrange(gamma_78Sites, theta_78Sites, z_2017_78Sites, "",
                  mapList[[1]], mapList[[2]], mapList[[3]], mapList[[4]],
                  ncol = 4, nrow = 2) %>%
  ggpubr::ggexport(filename = here::here(glue::glue("results/prediction/78-hmc/map_z0z50GammaTheta_78Sites_allPrices_det.png")),
                   width = 2700, height = 1500)



# LOOP ACROSS PRICE AND YEARS TO GENERATE Z DYNAMICS FOR EACH PRICE IN A SINGLE FIGURE
# years
aux.years <- c(2017, 2022, 2027, 2032, 2037, 2067)

# generate empty map list
for (p in seq_along(aux.prices)) {

  aux.mapList <- list()
  transfer <- aux.prices[p]-6.5

  for (y in seq_along(aux.years)) {

    # z50 map (vary by model)
    aux.mapList[[y]] <-
      ggplot2::ggplot(data = prediction.78SitesModel %>%
                        dplyr::filter(time == aux.years[y]-2017, p_e == aux.prices[p]) %>%
                        dplyr::mutate(z_t = cut(round(z_t,5),
                                                breaks = c(0,0.5,20,40,60,80,105),
                                                include.lowest = T,
                                                dig.lab = 3,
                                                labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")))) +
      ggplot2::geom_sf(aes(fill = z_t)) +
      ggplot2::scale_fill_manual(name = expr(paste("Z"[!!aux.years[y]]^"i", ~"(%), ", "b", "=", !!transfer)), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
      ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
      ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
      ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                     panel.grid.minor = element_line(colour = "white"),
                     panel.background = element_blank(),
                     strip.background = element_rect(fill = NA),
                     axis.line = element_blank(), axis.ticks = element_blank(),
                     axis.title = element_blank(), axis.text = element_blank(),
                     legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
                     legend.position = "bottom", legend.margin=margin(t=-1, r=-0, b=0.3, l=-0, unit="cm"),
                     legend.text = element_text(size = 28, face = "bold"),
                     plot.margin=unit(c(t=-0.5,r=-1,b=-0,l=-1), "cm"))

  }


  # SAVE MAP Z0, Z50, GAMMA, THETA
  ggpubr::ggarrange(z_2017_78Sites, aux.mapList[[2]], aux.mapList[[3]],
                    aux.mapList[[4]], aux.mapList[[5]], aux.mapList[[6]],
                    ncol = 3, nrow = 2) %>%
    ggpubr::ggexport(filename = here::here(glue::glue("results/prediction/78-hmc/map_zDecades_78Sites_pe{aux.prices[p]}_det.png")),
                     width = 2700, height = 1500)

}



# LOOP ACROSS YEARS TO GENERATE Z DYNAMICS IN MULTIPLE FIGURES WITH PRICE = 20.76
# years by decade
aux.years <- c(2017, 2022, 2027, 2032, 2037, 2042, 2047, 2052, 2057, 2062, 2067)

for (y in seq_along(aux.years)) {

  # z50 map (vary by model)
    ggplot2::ggplot(data = prediction.78SitesModel %>%
                      dplyr::filter(time == aux.years[y]-2017, p_e == 21.5) %>%
                      dplyr::mutate(z_t = cut(round(z_t,5),
                                              breaks = c(0,0.5,20,40,60,80,105),
                                              include.lowest = T,
                                              dig.lab = 3,
                                              labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")))) +
    ggplot2::geom_sf(aes(fill = z_t)) +
    ggplot2::scale_fill_manual(name = expr(paste("Z"[!!aux.years[y]]^"i", ~"(%), ", "P"^"e", "=", 22.5)), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
    ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
    ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
    ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                   panel.grid.minor = element_line(colour = "white"),
                   panel.background = element_blank(),
                   strip.background = element_rect(fill = NA),
                   axis.line = element_blank(), axis.ticks = element_blank(),
                   axis.title = element_blank(), axis.text = element_blank(),
                   legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
                   legend.position = "bottom", legend.margin=margin(t=-1, r=-0, b=0.3, l=-0, unit="cm"),
                   legend.text = element_text(size = 28, face = "bold"),
                   plot.margin=unit(c(t=-0.5,r=-1,b=-0,l=-1), "cm"))
  ggplot2::ggsave(here::here(glue::glue("results/prediction/78-hmc/map_z{aux.years[y]}_78Sites_pe20.76.png")), width = 12, height = 9, dpi = 300)

}



# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("code/analysis")






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------