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



# START TIMER
tictoc::tic(msg = "mapsPrediction_1043SitesModel.R script", log = T)

# OPTIONS
options(scipen = 999)

# DATA INPUT

# # GLOBAL MODEL CALIBRATION VARIABLES
# load(here::here("data/calibration/globalModel", "calibration_globalModel.Rdata"))

# extract high and low price
p_high <- 41.11

# # clear unnecessary objects
# rm(matrixTransition.2prices, calibration.globalModel)

# 1043 SITES MODEL CALIBRATION VARIABLES
load(here::here("data/calibration/", "calibration_1043_sites.Rdata"))


# 1043 SITES AGGREGATE PREDICTION
aux.prices <- c(6.6, 16.6, 21.6, 26.6, 31.6)




prediction.1043SitesModel <- purrr::map_df(
  .x = aux.prices,
  .f = ~ {
    file_path <- here::here("output/optimization/det/gams/1043sites/pa_41.11",
                            paste0("pe_", .x, "/Z.txt"))
    
    # Read the .txt file directly
    data <- readr::read_delim(file_path, 
                              delim = ",", 
                              col_names = FALSE,  # No column names
                              col_types = cols(.default = "n"))  # All numeric
    
    # Add a time column (T/R) ranging from 0 to 201
    time_period <- 0:(nrow(data) - 1)  # Assuming 202 time steps
    
    # Rename the columns to numeric ids (1, 2, ..., 1043)
    colnames(data) <- as.character(1:ncol(data))
    
    # Add the additional variables (T/R, p_e and p_a)
    data %>%
      dplyr::mutate(`T/R` = time_period, 
                    p_e = .x, 
                    p_a = p_high) %>%
      dplyr::relocate(`T/R`, .before = 1)  # Place T/R as the first column
  }
)




# AMAZON BIOME VECTOR DATA
load(here::here("data/clean/amazon_biome.Rdata"))

map_basin <- st_read(
  dsn = "data/raw/mapbiomas/basin",
  layer = "BASIN_LEVEL_2_PNRH"
)
calib_df <- calib_df %>%
  st_transform(st_crs(map_basin))

amazon_biome <- amazon_biome %>%
  st_transform(st_crs(map_basin))




gamma_fit <- read.csv(here::here("data/calibration/hmc", "gamma_fit_1043.csv"))%>%
  mutate(id = row_number())
theta_fit <- read.csv(here::here("data/calibration/hmc", "theta_fit_1043.csv"))%>%
  mutate(id = row_number())
calib_df <- calib_df %>%
  left_join(gamma_fit, by = "id")%>%
  left_join(theta_fit, by="id")
calib_df <- calib_df %>%
  mutate(
    x_1995 = calib_df$gamma_fit * (zbar_1995 - z_1995),
    x_2008 = calib_df$gamma_fit * (zbar_2008 - z_2008),
    x_2017 = calib_df$gamma_fit * (zbar_2017 - z_2017)
  )



# DATASET CLEANUP AND PREP

# select relevant calibrated variables
calib_df <-
  calib_df %>%
  dplyr::mutate(
    rank_theta_1043Sites = dense_rank(desc(theta_fit)), # rank values
    rank_gamma_1043Sites = dense_rank(desc(gamma_fit))
  ) %>% # rank values
  dplyr::select(id, z_2017, theta, rank_theta_1043Sites, rank_gamma_1043Sites, x_2017, zbar_2017)


# adjust predicted data to site-year panel + add calibrated variables
prediction.1043SitesModel <-
  prediction.1043SitesModel %>%
  tidyr::pivot_longer(cols = -c("T/R", "p_a", "p_e"), names_to = "id", values_to = "z_t") %>%
  dplyr::mutate(id = as.numeric(stringr::str_trim(id))) %>%
  dplyr::rename(time = "T/R") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
  dplyr::right_join(calib_df, by = c("id" = "id")) %>% # match by id, guarantee that prediction data uses the same id than the calibrated data
  sf::st_as_sf() %>%
  dplyr::mutate(
    z_t =  1e11 * z_t / zbar_2017,
    z_2017_1043Sites =  100* z_2017 / zbar_2017
  ) # transform share to %

# clean environment
rm(calib_df)

# adjust projection
prediction.1043SitesModel <- sf::st_as_sf(prediction.1043SitesModel)
prediction.1043SitesModel <- sf::st_transform(prediction.1043SitesModel, sf::st_crs(4326))
amazon_biome <- sf::st_transform(amazon_biome, sf::st_crs(prediction.1043SitesModel))


# GENERATE MAPS

# PLOT z_2017, GAMMA, NAD THETA CALIBRATED VALUES

# z_2017_1043Sites
z_2017_1043Sites <-
  ggplot2::ggplot(data = prediction.1043SitesModel %>%
    dplyr::filter(time == 0, p_e == 21.6) %>%
    dplyr::mutate(z_t = cut(round(z_t, 1),
      breaks = c(0, 0.5, 20, 40, 60, 80, 105),
      include.lowest = T,
      dig.lab = 3,
      labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")
    ))) +
  ggplot2::geom_sf(aes(fill = z_t)) +
  ggplot2::scale_fill_manual(name = expr(paste("Z"[!!2017]^"i", ~"(%)")), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
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
  ggplot2::ggplot(data = prediction.1043SitesModel %>% dplyr::filter(time == 0, p_e == 6.6) %>%
    dplyr::mutate(rank_gamma_1043Sites = cut(round(rank_gamma_1043Sites), breaks = c(1, 212, 423, 634, 845, 1043), include.lowest = T, dig.lab = 4))) +
  ggplot2::geom_sf(aes(fill = rank_gamma_1043Sites)) +
  ggplot2::scale_fill_brewer(name = expression(paste(gamma^"i", ~"(rank)")), palette = "YlOrRd", direction = -1) +
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
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
  ggplot2::ggplot(data = prediction.1043SitesModel %>% dplyr::filter(time == 0, p_e == 6.6) %>%
    dplyr::mutate(rank_theta_1043Sites = cut(round(rank_theta_1043Sites), breaks = c(1, 212, 423, 634, 845, 1043), include.lowest = T, dig.lab = 4))) +
  ggplot2::geom_sf(aes(fill = rank_theta_1043Sites)) +
  ggplot2::scale_fill_brewer(name = expression(paste(theta^"i", ~"(rank)")), palette = "YlOrRd", direction = -1) +
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
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
  transfer <- price - 6.6
  # z50 map (vary by model)
  mapList[[mapIndex]] <-
    ggplot2::ggplot(data = prediction.1043SitesModel %>%
      dplyr::filter(time == 30, p_e == price) %>%
      dplyr::mutate(z_t = cut(round(z_t, 5),
        breaks = c(0, 0.5, 20, 40, 60, 80, 105),
        include.lowest = T,
        dig.lab = 3,
        labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")
      ))) +
    ggplot2::geom_sf(aes(fill = z_t)) +
    ggplot2::scale_fill_manual(name = expr(paste("Z"[2047]^"i", ~"(%), ", "b", "=", !!transfer)), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
    ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
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
      filename = here::here(glue::glue("plots/1043-det/map_z50_1043Sites_pe{price}.png")),
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
    filename = here::here(glue::glue("plots/1043-det/map_z0z50GammaTheta_1043Sites_allPrices_det.png")),
    width = 2700, height = 1500
  )


# FIGURE WTIH ALL MAPS
# SAVE MAP Z0, GAMMA, THETA, with varying prices
ggpubr::ggarrange(z_2017_1043Sites,
  mapList[[1]], mapList[[2]], mapList[[5]],
  ncol = 4, nrow = 1
) %>%
  ggpubr::ggexport(
    filename = here::here(glue::glue("plots/1043-det/map_z0z30GammaTheta_1043Sites_allPrices_det.png")),
    width = 2700, height = 800
  )

# LOOP ACROSS PRICE AND YEARS TO GENERATE Z DYNAMICS FOR EACH PRICE IN A SINGLE FIGURE
# years
aux.years <- c(2017, 2022, 2027, 2032, 2037, 2047)

# generate empty map list
for (p in seq_along(aux.prices)) {
  aux.mapList <- list()
  transfer <- aux.prices[p] - 6.6

  for (y in seq_along(aux.years)) {
    # z50 map (vary by model)
    aux.mapList[[y]] <-
      ggplot2::ggplot(data = prediction.1043SitesModel %>%
        dplyr::filter(time == aux.years[y] - 2017, p_e == aux.prices[p]) %>%
        dplyr::mutate(z_t = cut(round(z_t, 5),
          breaks = c(0, 0.5, 20, 40, 60, 80, 105),
          include.lowest = T,
          dig.lab = 3,
          labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")
        ))) +
      ggplot2::geom_sf(aes(fill = z_t)) +
      ggplot2::scale_fill_manual(name = expr(paste("Z"[!!aux.years[y]]^"i", ~"(%), ", "b", "=", !!transfer)), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
      ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
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
      filename = here::here(glue::glue("plots/1043-det/map_zDecades_1043Sites_pe{aux.prices[p]}_det.png")),
      width = 2700, height = 1500
    )
}



# END TIMER
tictoc::toc(log = T)
