
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
# GLOBAL MODEL CALIBRATION VARIABLES
# load(here::here("data/calibration/globalModel", "calibration_globalModel.Rdata"))

# extract high and low price
p_high <- 42.03

# # clear unnecessary objects
# rm(matrixTransition.2prices, calibration.globalModel)


# 78 SITES MODEL CALIBRATION VARIABLES
load(here::here("data/hmc", "hmc_78SitesModel.Rdata"))


# 78 SITES AGGREGATE PREDICTION
aux.prices <- c(7, 17, 22, 27)

prediction.78SitesModel <- purrr::map_df(.x = aux.prices,
                                         .f = ~readr::read_delim(here::here("data/prediction/78-mpc",
                                                                            glue::glue("p_e_{.x}/amazon_data_z.dat")), skip = 3,
                                                                 col_types = cols(.default = "n"), col_names = c("Y",  "z_t")) %>%
                                           dplyr::mutate(p_e = .x)) %>%
                            dplyr::mutate(id = c(1:78, rep(0:78, times = 199),1:78, rep(0:78, times = 199),1:78, rep(0:78, times = 199),1:78, rep(0:78, times = 199)))


# AMAZON BIOME VECTOR DATA
load(here::here("data/raw2clean/amazonBiome_ibge/output/clean_amazonBiome.Rdata"))



# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# select relevant calibrated variables
calibration.78SitesModel <-
  calibration.78SitesModel %>%
  dplyr::mutate(rank_theta_78Sites = dense_rank(desc(theta_78Sites)), # rank values
                rank_gamma_78Sites = dense_rank(desc(gamma_78Sites))) %>%  # rank values
  dplyr::select(id, z_2017_78Sites, rank_theta_78Sites, rank_gamma_78Sites, theta_78Sites, x_2017_78Sites, zbar_2017_78Sites)

# adjust predicted data to site-year panel + add calibrated variables
prediction.78SitesModel <-
  prediction.78SitesModel %>%
  dplyr::filter(is.na(Y))  %>% # remove rows indicating change in Y (time-period)
  dplyr::mutate(time = rep(rep(1:200, each = 78), length(aux.prices))) %>% # manually add time variable
  dplyr::mutate(time = time-1) %>% # adjust time variable to match with 1000 sites results
  dplyr::select(id, time, z_t, p_e) %>% # select columns of interest
  dplyr::arrange(p_e, time, id) %>%
  dplyr::right_join(calibration.78SitesModel, by = c("id" = "id")) %>% # match by id, guarantee that prediction data uses the same id than the calibrated data
  sf::st_as_sf() %>%
  dplyr::mutate(z_t = 100*z_t/zbar_2017_78Sites,
                z_2017_78Sites = 100*z_2017_78Sites/zbar_2017_78Sites) # transform share to %

# clean environment
rm(calibration.78SitesModel)

# adjust projection
prediction.78SitesModel <- sf::st_as_sf(prediction.78SitesModel)
prediction.78SitesModel <- sf::st_transform(prediction.78SitesModel, sf::st_crs(4326))
clean.amazonBiome <- sf::st_transform(clean.amazonBiome, sf::st_crs(prediction.78SitesModel))

# add centroids
prediction.78SitesModel <- cbind(prediction.78SitesModel, st_coordinates(st_centroid(prediction.78SitesModel)))





# GENERATE MAPS --------------------------------------------------------------------------------------------------------------------------------------

# PLOT z_2017, GAMMA, THETA CALIBRATED VALUES AND z_50 PREDICTED VALUE

# z_2017_78Sites
z_2017_78Sites <-
  ggplot2::ggplot(data = prediction.78SitesModel %>% dplyr::filter(time == 1, p_e == 7) %>%
                    dplyr::mutate(z_2017_78Sites = cut(round(z_2017_78Sites, 2), breaks = c(0,0.5,20,40,60,80,105), include.lowest = T, dig.lab = 3,
                                                        labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")))) +
  ggplot2::geom_sf(aes(fill = z_2017_78Sites)) +
  ggplot2::scale_fill_manual(name = expression(paste("Z"[2017]^"i", ~"(%)")), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
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

zNumber_2017_78Sites <-
  ggplot2::ggplot(data = prediction.78SitesModel %>% dplyr::filter(time == 1, p_e == 7) %>%
                    dplyr::mutate(z_2017_78Sites = round(z_2017_78Sites))) +
  ggplot2::geom_sf(fill = NA) +
  ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1) +
  ggplot2::geom_text(aes(X, Y, label = z_2017_78Sites), size = 10) +
  ggplot2::geom_text(data = prediction.78SitesModel %>%
                       dplyr::filter(time == 50, p_e == 7, rank_theta_78Sites == 8, rank_gamma_78Sites == 78) %>%
                       dplyr::mutate(z_2017_78Sites = round(z_2017_78Sites)),
                     aes(X, Y, label = z_2017_78Sites), size = 10, fontface = "bold") +
  ggplot2::ggtitle(expression(paste("Z"[2017]^"i", ~"(%)"))) +
  ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                 panel.grid.minor = element_line(colour = "white"),
                 panel.background = element_blank(),
                 strip.background = element_rect(fill = NA),
                 axis.line = element_blank(), axis.ticks = element_blank(),
                 axis.title = element_blank(), axis.text = element_blank(),
                 plot.title = element_text(hjust = 0.5, size = 36, face = "bold"),
                 plot.margin=unit(c(t=0.5,r=-1,b=-0,l=-1), "cm"))

# gamma_78Sites
gamma_78Sites <-
  ggplot2::ggplot(data = prediction.78SitesModel %>% dplyr::filter(time == 1, p_e == 7) %>%
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

gammaNumber_78Sites <-
  ggplot2::ggplot(data = prediction.78SitesModel %>% dplyr::filter(time == 1, p_e == 7) %>%
                    dplyr::mutate(rank_gamma_78Sites = round(rank_gamma_78Sites))) +
  ggplot2::geom_sf(fill = NA) +
  ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1) +
  ggplot2::geom_text(aes(X, Y, label = rank_gamma_78Sites), size = 10) +
  ggplot2::geom_text(data = prediction.78SitesModel %>%
                       dplyr::filter(time == 50, p_e ==7, rank_theta_78Sites == 8, rank_gamma_78Sites == 78) %>%
                       dplyr::mutate(rank_gamma_78Sites = round(rank_gamma_78Sites)),
                     aes(X, Y, label = rank_gamma_78Sites), size = 10, fontface = "bold") +
  ggplot2::ggtitle(expression(paste(gamma^"i", ~"(rank)"))) +
  ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                 panel.grid.minor = element_line(colour = "white"),
                 panel.background = element_blank(),
                 strip.background = element_rect(fill = NA),
                 axis.line = element_blank(), axis.ticks = element_blank(),
                 axis.title = element_blank(), axis.text = element_blank(),
                 plot.title = element_text(hjust = 0.5, size = 36, face = "bold"),
                 plot.margin=unit(c(t=0.5,r=-1,b=-0,l=-1), "cm"))

# theta_78Sites
theta_78Sites <-
  ggplot2::ggplot(data = prediction.78SitesModel %>% dplyr::filter(time == 1, p_e == 7) %>%
                    dplyr::mutate(rank_theta_78Sites = ggplot2::cut_interval(round(rank_theta_78Sites),  n = 5, dig.lab = 3))) +
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

thetaNumber_78Sites <-
  ggplot2::ggplot(data = prediction.78SitesModel %>% dplyr::filter(time == 1, p_e == 7) %>%
                    dplyr::mutate(rank_theta_78Sites = round(rank_theta_78Sites))) +
  ggplot2::geom_sf(fill = NA) +
  ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1) +
  ggplot2::geom_text(aes(X, Y, label = rank_theta_78Sites), size = 10) +
  ggplot2::geom_text(data = prediction.78SitesModel %>%
                       dplyr::filter(time == 50, p_e == 7, rank_theta_78Sites == 8, rank_gamma_78Sites == 78) %>%
                       dplyr::mutate(rank_theta_78Sites = round(rank_theta_78Sites)),
                     aes(X, Y, label = rank_theta_78Sites), size = 10, fontface = "bold") +
  ggplot2::ggtitle(expression(paste(theta^"i", ~"(rank)"))) +
  ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                 panel.grid.minor = element_line(colour = "white"),
                 panel.background = element_blank(),
                 strip.background = element_rect(fill = NA),
                 axis.line = element_blank(), axis.ticks = element_blank(),
                 axis.title = element_blank(), axis.text = element_blank(),
                 plot.title = element_text(hjust = 0.5, size = 36, face = "bold"),
                 plot.margin=unit(c(t=0.5,r=-1,b=-0,l=-1), "cm"))



# LOOP ACROSS PRICES TO GENERATE Z50 MAP AND SAVE Z0, Z50, GAMMA, AND THETA MAPS BY PRICE (WITH NUMBERS DISPLAYED)
mapList <- list()
mapIndex <- 1

for (price in aux.prices) {
transfer <- price - 7
  # z50 map (vary by model)
  mapList[[mapIndex]] <-
    ggplot2::ggplot(data = prediction.78SitesModel %>%
                      dplyr::filter(time == 50, p_e == price) %>%
                      dplyr::mutate(z_t = round(z_t))) +
    ggplot2::geom_sf(fill = NA) +
    ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1) +
    ggplot2::geom_text(aes(X, Y, label = z_t), size = 10) +
    ggplot2::geom_text(data = prediction.78SitesModel %>%
                         dplyr::filter(time == 50, p_e == price, rank_theta_78Sites == 8, rank_gamma_78Sites == 78) %>%
                         dplyr::mutate(z_t = round(z_t)),
                       aes(X, Y, label = z_t), size = 10, fontface = "bold") +
    ggplot2::ggtitle(expr(paste("Z"[2067]^"i", ~"(%), ", "b", "=", !!transfer))) +
    ggplot2::theme(panel.grid.major = element_line(colour = "White"),
                   panel.grid.minor = element_line(colour = "white"),
                   panel.background = element_blank(),
                   strip.background = element_rect(fill = NA),
                   axis.line = element_blank(), axis.ticks = element_blank(),
                   axis.title = element_blank(), axis.text = element_blank(),
                   plot.title = element_text(hjust = 0.5, size = 36, face = "bold"),
                   plot.margin=unit(c(t=0.5,r=-1,b=-0,l=-1), "cm"))

  # SAVE MAP Z0, Z50, GAMMA, THETA
  ggpubr::ggarrange(zNumber_2017_78Sites, mapList[[mapIndex]], gammaNumber_78Sites, thetaNumber_78Sites, ncol = 2, nrow = 2) %>%
    ggpubr::ggexport(filename = here::here(glue::glue("results/prediction/78-MPC/map_z0z50GammaThetaNumber_78Sites_pe{price}.png")),
                     width = 2400, height = 1500)

  mapIndex <- mapIndex + 1

}




########## Map in 30 years


zNumber_2017_78Sites <-
  ggplot2::ggplot(data = prediction.78SitesModel %>% dplyr::filter(time == 1, p_e == 7) %>%
                    dplyr::mutate(z_2017_78Sites = round(z_2017_78Sites))) +
  ggplot2::geom_sf(fill = NA) +
  ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1) +
  ggplot2::geom_text(aes(X, Y, label = z_2017_78Sites, color = ifelse(rank_theta_78Sites == 2, "red", ifelse(rank_theta_78Sites == 1, "#3333FF", "black"))), size = 10,fontface = "bold") +
  ggplot2::scale_color_identity()+
  ggplot2::geom_text(data = prediction.78SitesModel %>%
                       dplyr::filter(time == 50, p_e == 7, rank_theta_78Sites == 8, rank_gamma_78Sites == 78) %>%
                       dplyr::mutate(z_2017_78Sites = round(z_2017_78Sites)),
                     aes(X, Y, label = z_2017_78Sites), size = 10, fontface = "bold") +
  ggplot2::ggtitle(expression(paste("Z"[2017]^"i", ~"(%)"))) +
  ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                 panel.grid.minor = element_line(colour = "white"),
                 panel.background = element_blank(),
                 strip.background = element_rect(fill = NA),
                 axis.line = element_blank(), axis.ticks = element_blank(),
                 axis.title = element_blank(), axis.text = element_blank(),
                 plot.title = element_text(hjust = 0.5, size = 36, face = "bold"),
                 plot.margin=unit(c(t=0.5,r=-1,b=-0,l=-1), "cm"))

gammaNumber_78Sites <-
  ggplot2::ggplot(data = prediction.78SitesModel %>% dplyr::filter(time == 1, p_e == 7) %>%
                    dplyr::mutate(rank_gamma_78Sites = round(rank_gamma_78Sites))) +
  ggplot2::geom_sf(fill = NA) +
  ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1) +
  ggplot2::geom_text(aes(X, Y, label = rank_gamma_78Sites, color = ifelse(rank_theta_78Sites == 2, "red", ifelse(rank_theta_78Sites == 1, "#3333FF", "black"))), size = 10,fontface = "bold") +
  ggplot2::scale_color_identity()+
  ggplot2::geom_text(data = prediction.78SitesModel %>%
                       dplyr::filter(time == 50, p_e ==7, rank_theta_78Sites == 8, rank_gamma_78Sites == 78) %>%
                       dplyr::mutate(rank_gamma_78Sites = round(rank_gamma_78Sites)),
                     aes(X, Y, label = rank_gamma_78Sites), size = 10, fontface = "bold") +
  ggplot2::ggtitle(expression(paste(gamma^"i", ~"(rank)"))) +
  ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                 panel.grid.minor = element_line(colour = "white"),
                 panel.background = element_blank(),
                 strip.background = element_rect(fill = NA),
                 axis.line = element_blank(), axis.ticks = element_blank(),
                 axis.title = element_blank(), axis.text = element_blank(),
                 plot.title = element_text(hjust = 0.5, size = 36, face = "bold"),
                 plot.margin=unit(c(t=0.5,r=-1,b=-0,l=-1), "cm"))


thetaNumber_78Sites <-
  ggplot2::ggplot(data = prediction.78SitesModel %>% dplyr::filter(time == 1, p_e == 7) %>%
                    dplyr::mutate(rank_theta_78Sites = round(rank_theta_78Sites))) +
  ggplot2::geom_sf(fill = NA) +
  ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1) +
  ggplot2::geom_text(aes(X, Y, label = rank_theta_78Sites, color = ifelse(rank_theta_78Sites == 2, "red", ifelse(rank_theta_78Sites == 1, "#3333FF", "black"))), size = 10,fontface = "bold") +
  ggplot2::scale_color_identity()+
  ggplot2::geom_text(data = prediction.78SitesModel %>%
                       dplyr::filter(time == 50, p_e == 7, rank_theta_78Sites == 8, rank_gamma_78Sites == 78) %>%
                       dplyr::mutate(rank_theta_78Sites = round(rank_theta_78Sites)),
                     aes(X, Y, label = rank_theta_78Sites), size = 10, fontface = "bold") +
  ggplot2::ggtitle(expression(paste(theta^"i", ~"(rank)"))) +
  ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                 panel.grid.minor = element_line(colour = "white"),
                 panel.background = element_blank(),
                 strip.background = element_rect(fill = NA),
                 axis.line = element_blank(), axis.ticks = element_blank(),
                 axis.title = element_blank(), axis.text = element_blank(),
                 plot.title = element_text(hjust = 0.5, size = 36, face = "bold"),
                 plot.margin=unit(c(t=0.5,r=-1,b=-0,l=-1), "cm"))


mapList <- list()
mapIndex <- 1

for (price in aux.prices) {
transfer <- price - 7
  # z50 map (vary by model)
  mapList[[mapIndex]] <-
    ggplot2::ggplot(data = prediction.78SitesModel %>%
                      dplyr::filter(time == 30, p_e == price) %>%
                      dplyr::mutate(z_t = round(z_t))) +
    ggplot2::geom_sf(fill = NA) +
    ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1) +
    ggplot2::geom_text(aes(X, Y, label = z_t, color = ifelse(rank_theta_78Sites == 2, "red", ifelse(rank_theta_78Sites == 1, "#3333FF", "black"))), size = 10,fontface = "bold") +
    ggplot2::scale_color_identity()+
    ggplot2::geom_text(data = prediction.78SitesModel %>%
                         dplyr::filter(time == 30, p_e == price, rank_theta_78Sites == 8, rank_gamma_78Sites == 78) %>%
                         dplyr::mutate(z_t = round(z_t)),
                       aes(X, Y, label = z_t), size = 10, fontface = "bold") +
    ggplot2::ggtitle(expr(paste("Z"[2047]^"i", ~"(%), ", "b", "=", !!transfer))) +
    ggplot2::theme(panel.grid.major = element_line(colour = "White"),
                   panel.grid.minor = element_line(colour = "white"),
                   panel.background = element_blank(),
                   strip.background = element_rect(fill = NA),
                   axis.line = element_blank(), axis.ticks = element_blank(),
                   axis.title = element_blank(), axis.text = element_blank(),
                   plot.title = element_text(hjust = 0.5, size = 36, face = "bold"),
                   plot.margin=unit(c(t=0.5,r=-1,b=-0,l=-1), "cm"))

  # SAVE MAP Z0, Z50, GAMMA, THETA
  ggpubr::ggarrange(zNumber_2017_78Sites, mapList[[mapIndex]], gammaNumber_78Sites, thetaNumber_78Sites, ncol = 2, nrow = 2) %>%
    ggpubr::ggexport(filename = here::here(glue::glue("results/prediction/78-MPC/map_z0z30GammaThetaNumber_78Sites_pe{price}.png")),
                     width = 2400, height = 1500)

  mapIndex <- mapIndex + 1

}

######### End map here

# LOOP ACROSS PRICES TO GENERATE Z50 MAP AND SAVE Z0, Z50, GAMMA, AND THETA MAPS BY PRICE
mapList <- list()
mapIndex <- 1
for (price in aux.prices) {
  transfer<-price-7
  # z50 map (vary by model)
  mapList[[mapIndex]] <-
    ggplot2::ggplot(data = prediction.78SitesModel %>%
                      dplyr::filter(time == 50, p_e == price) %>%
                      dplyr::mutate(z_t = cut(round(z_t,2),
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
    ggpubr::ggexport(filename = here::here(glue::glue("results/prediction/78-MPC/map_z50_78Sites_pe{price}.png")),
                     width = 2400, height = 1500)

  mapIndex <- mapIndex + 1

}

# FIGURE WTIH ALL MAPS
# SAVE MAP Z0, GAMMA, THETA, with varying prices
ggpubr::ggarrange(gamma_78Sites, theta_78Sites, z_2017_78Sites, "",
                  mapList[[1]], mapList[[2]], mapList[[3]], mapList[[4]],
                  ncol = 4, nrow = 2) %>%
  ggpubr::ggexport(filename = here::here(glue::glue("results/prediction/78-MPC/map_z0z50GammaTheta_78Sites_allPrices_mpc.png")),
                   width = 2700, height = 1500)



# LOOP ACROSS PRICE AND YEARS TO GENERATE Z DYNAMICS FOR EACH PRICE IN A SINGLE FIGURE
# years
aux.years <- c(2017, 2022, 2027, 2032, 2037, 2067)

# generate empty map list
for (p in seq_along(aux.prices)) {

  aux.mapList <- list()
  transfer <- aux.prices[p]-7
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
    ggpubr::ggexport(filename = here::here(glue::glue("results/prediction/78-MPC/map_zDecades_78Sites_pe{aux.prices[p]}_mpc.png")),
                     width = 2700, height = 1500)

}



# LOOP ACROSS YEARS TO GENERATE Z DYNAMICS IN MULTIPLE FIGURES WITH PRICE = 20.76
# years by decade
aux.years <- c(2017, 2022, 2027, 2032, 2037, 2042, 2047, 2052, 2057, 2062, 2067)

for (y in seq_along(aux.years)) {

  # z50 map (vary by model)
  ggplot2::ggplot(data = prediction.78SitesModel %>%
                    dplyr::filter(time == aux.years[y]-2017, p_e == 22) %>%
                    dplyr::mutate(z_t = cut(round(z_t,5),
                                            breaks = c(0,0.5,20,40,60,80,105),
                                            include.lowest = T,
                                            dig.lab = 3,
                                            labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")))) +
    ggplot2::geom_sf(aes(fill = z_t)) +
    ggplot2::scale_fill_manual(name = expr(paste("Z"[!!aux.years[y]]^"i", ~"(%), ", "P"^"e", "=", 22)), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
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
  ggplot2::ggsave(here::here(glue::glue("results/prediction/78-MPC/map_z{aux.years[y]}_78Sites_pe20.76.png")), width = 12, height = 9, dpi = 300)

}




# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("code/analysis")






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------