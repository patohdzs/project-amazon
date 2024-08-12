# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: MAPS OF 1043 SITES MODEL PREDICTED VALUES
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

# SETUP

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")

# START TIMER
tictoc::tic(msg = "mapsPrediction_1043SitesModel.R script", log = T)

# OPTIONS
options(scipen = 999)

# DATA INPUT

# # GLOBAL MODEL CALIBRATION VARIABLES
# load(here::here("data/calibration/globalModel", "calibration_globalModel.Rdata"))

# extract high and low price
p_high <- 44.76
p_low <- 38.30

# # clear unnecessary objects
# rm(matrixTransition.2prices, calibration.globalModel)

# 1043 SITES MODEL CALIBRATION VARIABLES
load(here::here("data/hmc", "hmc_1043SitesModel.Rdata"))

# 1043 SITES AGGREGATE PREDICTION
aux.prices <- c(7.5, 17.5, 22.5, 27.5)

# read the site-level predictions with high cattle price to calculate aggregated agricultural output
prediction.1043SitesModel <- purrr::map_df(
  .x = aux.prices,
  .f = ~ readr::read_delim(
    here::here(
      "data/prediction/1043SitesModel",
      glue::glue("ph_p_e_{.x}/amazon_data_z.dat")
    ),
    col_types = cols(.default = "n")
  ) %>%
    dplyr::mutate(p_e = .x, p_a = p_high)
)

# # read the site-level predictions with low cattle price to calculate aggregated agricultural output
# prediction.1043SitesModel <- purrr::map_df(.x = aux.prices,
#                                            .f = ~readr::read_delim(here::here("data/prediction/1043SitesModel",
#                                                                               glue::glue("pl_p_e_{.x}/amazon_data_z.dat")),
#                                                                    col_types = cols(.default = "n")) %>%
#                                              dplyr::mutate(p_e = .x, p_a = p_low)) %>%
#   bind_rows(prediction.1043SitesModel)

# AMAZON BIOME VECTOR DATA
load(here::here("data/raw2clean/amazonBiome_ibge/output/clean_amazonBiome.Rdata"))

# DATASET CLEANUP AND PREP

# select relevant calibrated variables
calibration.1043SitesModel <-
  calibration.1043SitesModel %>%
  dplyr::mutate(
    rank_theta_1043Sites = dense_rank(desc(theta_1043Sites)), # rank values
    rank_gamma_1043Sites = dense_rank(desc(gamma_1043Sites))
  ) %>% # rank values
  dplyr::select(id, z_2017_1043Sites, theta_1043Sites, rank_theta_1043Sites, rank_gamma_1043Sites, x_2017_1043Sites, zbar_2017_1043Sites)

# adjust predicted data to site-year panel + add calibrated variables
prediction.1043SitesModel <-
  prediction.1043SitesModel %>%
  tidyr::pivot_longer(cols = -c("T/R ", "p_a", "p_e"), names_to = "id", values_to = "z_t") %>%
  dplyr::mutate(id = as.numeric(stringr::str_trim(id))) %>%
  dplyr::rename(time = "T/R ") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
  dplyr::right_join(calibration.1043SitesModel, by = c("id" = "id")) %>% # match by id, guarantee that prediction data uses the same id than the calibrated data
  sf::st_as_sf() %>%
  dplyr::mutate(
    z_t = 100 * 1e11 * z_t / zbar_2017_1043Sites,
    z_2017_1043Sites = 100 * z_2017_1043Sites / zbar_2017_1043Sites
  ) # transform share to %

# clean environment
rm(calibration.1043SitesModel)

# adjust projection
prediction.1043SitesModel <- sf::st_as_sf(prediction.1043SitesModel)
prediction.1043SitesModel <- sf::st_transform(prediction.1043SitesModel, sf::st_crs(4326))
clean.amazonBiome <- sf::st_transform(clean.amazonBiome, sf::st_crs(prediction.1043SitesModel))

# GENERATE MAPS

# PLOT z_2017, GAMMA, NAD THETA CALIBRATED VALUES

# z_2017_1043Sites
z_2017_1043Sites <-
  ggplot2::ggplot(data = prediction.1043SitesModel %>% dplyr::filter(time == 0, p_e == 7.5) %>%
    dplyr::mutate(z_2017_1043Sites = cut(round(z_2017_1043Sites, 1),
      breaks = c(0, 0.5, 20, 40, 60, 80, 105), include.lowest = T, dig.lab = 3,
      labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")
    ))) +
  ggplot2::geom_sf(aes(fill = z_2017_1043Sites)) +
  ggplot2::scale_fill_manual(name = expression(paste("Z"[2017]^"i", ~"(%)")), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
  ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
    legend.position = "bottom", legend.margin = margin(t = -1, r = -0, b = 0.3, l = -0, unit = "cm"),
    legend.text = element_text(size = 28, face = "bold"),
    plot.margin = unit(c(t = -0.5, r = -1, b = -0, l = -1), "cm")
  )

# gamma_1043Sites
gamma_1043Sites <-
  ggplot2::ggplot(data = prediction.1043SitesModel %>% dplyr::filter(time == 0, p_e == 7.5) %>%
    dplyr::mutate(rank_gamma_1043Sites = cut(round(rank_gamma_1043Sites), breaks = c(1, 212, 423, 634, 845, 1043), include.lowest = T, dig.lab = 4))) +
  ggplot2::geom_sf(aes(fill = rank_gamma_1043Sites)) +
  ggplot2::scale_fill_brewer(name = expression(paste(gamma^"i", ~"(rank)")), palette = "YlOrRd", direction = -1) +
  ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
    legend.position = "bottom", legend.margin = margin(t = -1, r = -0, b = 0.3, l = -0, unit = "cm"),
    legend.text = element_text(size = 28, face = "bold"),
    plot.margin = unit(c(t = -0.5, r = -1, b = -0, l = -1), "cm")
  )

# theta_1043Sites
theta_1043Sites <-
  ggplot2::ggplot(data = prediction.1043SitesModel %>% dplyr::filter(time == 0, p_e == 7.5) %>%
    dplyr::mutate(rank_theta_1043Sites = cut(round(rank_theta_1043Sites), breaks = c(1, 212, 423, 634, 845, 1043), include.lowest = T, dig.lab = 4))) +
  ggplot2::geom_sf(aes(fill = rank_theta_1043Sites)) +
  ggplot2::scale_fill_brewer(name = expression(paste(theta^"i", ~"(rank)")), palette = "YlOrRd", direction = -1) +
  ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
    legend.position = "bottom", legend.margin = margin(t = -1, r = -0, b = 0.3, l = -0, unit = "cm"),
    legend.text = element_text(size = 28, face = "bold"),
    plot.margin = unit(c(t = -0.5, r = -1, b = -0, l = -1), "cm")
  )

# LOOP ACROSS PRICES TO GENERATE Z50 MAP AND SAVE Z0, Z50, GAMMA, AND THETA MAPS BY PRICE
mapList <- list()
mapIndex <- 1
for (price in aux.prices) {
  # z50 map (vary by model)
  mapList[[mapIndex]] <-
    ggplot2::ggplot(data = prediction.1043SitesModel %>%
      dplyr::filter(time == 50, p_e == price) %>%
      dplyr::mutate(z_t = cut(round(z_t, 2),
        breaks = c(0, 0.01, 20, 40, 60, 80, 105),
        include.lowest = T,
        dig.lab = 3,
        labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")
      ))) +
    ggplot2::geom_sf(aes(fill = z_t)) +
    ggplot2::scale_fill_manual(name = expr(paste("Z"[2067]^"i", ~"(%), ", "P"^"e", "=", !!price)), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
    ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
    ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
    ggplot2::theme(
      panel.grid.major = element_line(colour = "white"),
      panel.grid.minor = element_line(colour = "white"),
      panel.background = element_blank(),
      strip.background = element_rect(fill = NA),
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(), axis.text = element_blank(),
      legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
      legend.position = "bottom", legend.margin = margin(t = -1, r = -0, b = 0.3, l = -0, unit = "cm"),
      legend.text = element_text(size = 28, face = "bold"),
      plot.margin = unit(c(t = -0.5, r = -1, b = -0, l = -1), "cm")
    )

  # SAVE MAP Z0, Z50, GAMMA, THETA
  ggpubr::ggarrange(mapList[[mapIndex]]) %>%
    ggpubr::ggexport(
      filename = here::here(glue::glue("plots/prediction/1043SitesModel/map_z50_1043Sites_pe{price}.png")),
      width = 2400, height = 1500
    )

  mapIndex <- mapIndex + 1
}

# FIGURE WTIH ALL MAPS
# SAVE MAP Z0, GAMMA, THETA, with varying prices
ggpubr::ggarrange(gamma_1043Sites, theta_1043Sites, z_2017_1043Sites, "",
  mapList[[1]], mapList[[2]], mapList[[3]], mapList[[4]],
  ncol = 4, nrow = 2
) %>%
  ggpubr::ggexport(
    filename = here::here(glue::glue("plots/prediction/1043SitesModel/map_z0z50GammaTheta_1043Sites_allPrices.png")),
    width = 2700, height = 1500
  )

# LOOP ACROSS PRICE AND YEARS TO GENERATE Z DYNAMICS FOR EACH PRICE IN A SINGLE FIGURE
# years
aux.years <- c(2017, 2022, 2027, 2032, 2037, 2067)

# generate empty map list
for (p in seq_along(aux.prices)) {
  aux.mapList <- list()

  for (y in seq_along(aux.years)) {
    # z50 map (vary by model)
    aux.mapList[[y]] <-
      ggplot2::ggplot(data = prediction.1043SitesModel %>%
        dplyr::filter(time == aux.years[y] - 2017, p_e == aux.prices[p]) %>%
        dplyr::mutate(z_t = cut(round(z_t, 5),
          breaks = c(0, 0.01, 20, 40, 60, 80, 105),
          include.lowest = T,
          dig.lab = 3,
          labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")
        ))) +
      ggplot2::geom_sf(aes(fill = z_t)) +
      ggplot2::scale_fill_manual(name = expr(paste("Z"[!!aux.years[y]]^"i", ~"(%), ", "P"^"e", "=", !!aux.prices[p])), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
      ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
      ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
      ggplot2::theme(
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA),
        axis.line = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), axis.text = element_blank(),
        legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
        legend.position = "bottom", legend.margin = margin(t = -1, r = -0, b = 0.3, l = -0, unit = "cm"),
        legend.text = element_text(size = 28, face = "bold"),
        plot.margin = unit(c(t = -0.5, r = -1, b = -0, l = -1), "cm")
      )
  }

  # SAVE MAP Z0, Z50, GAMMA, THETA
  ggpubr::ggarrange(z_2017_1043Sites, aux.mapList[[2]], aux.mapList[[3]],
    aux.mapList[[4]], aux.mapList[[5]], aux.mapList[[6]],
    ncol = 3, nrow = 2
  ) %>%
    ggpubr::ggexport(
      filename = here::here(glue::glue("plots/prediction/1043SitesModel/map_zDecades_1043Sites_pe{aux.prices[p]}.png")),
      width = 2700, height = 1500
    )
}

# LOOP ACROSS YEARS TO GENERATE Z DYNAMICS IN MULTIPLE FIGURES WITH PRICE = 20.76
# years by decade
aux.years <- c(2017, 2022, 2027, 2032, 2037, 2042, 2047, 2052, 2057, 2062, 2067)

for (y in seq_along(aux.years)) {
  # z50 map (vary by model)
  plot <- ggplot2::ggplot(data = prediction.1043SitesModel %>%
    dplyr::filter(time == aux.years[y] - 2017, p_e == 22.5) %>%
    dplyr::mutate(z_t = cut(round(z_t, 5),
      breaks = c(0, 0.01, 20, 40, 60, 80, 105),
      include.lowest = T,
      dig.lab = 3,
      labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")
    ))) +
    ggplot2::geom_sf(aes(fill = z_t)) +
    ggplot2::scale_fill_manual(name = expr(paste("Z"[!!aux.years[y]]^"i", ~"(%), ", "P"^"e", "=", 22.5)), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
    ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
    ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
    ggplot2::theme(
      panel.grid.major = element_line(colour = "white"),
      panel.grid.minor = element_line(colour = "white"),
      panel.background = element_blank(),
      strip.background = element_rect(fill = NA),
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(), axis.text = element_blank(),
      legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
      legend.position = "bottom", legend.margin = margin(t = -1, r = -0, b = 0.3, l = -0, unit = "cm"),
      legend.text = element_text(size = 28, face = "bold"),
      plot.margin = unit(c(t = -0.5, r = -1, b = -0, l = -1), "cm")
    )
  ggplot2::ggsave(here::here(glue::glue("plots/prediction/1043SitesModel/map_z{aux.years[y]}_1043Sites_pe20.76.png")), plot = plot, width = 12, height = 9, dpi = 300)
}

# END TIMER
tictoc::toc(log = T)
