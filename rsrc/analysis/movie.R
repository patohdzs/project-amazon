
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


install.packages("animation")
library(animation)
install.packages("png")
library(png)
install.packages("patchwork")
library(patchwork)
install.packages("cowplot")
library(cowplot)
conflicts_prefer(patchwork::area)

# SETUP ----------------------------------------------------------------------------------------------------------------------------------------------

# RUN 'setup.R' TO CONFIGURE INITIAL SETUP (mostly installing/loading packages)
source("rsrc/setup.R")


# START TIMER
tictoc::tic(msg = "mapsPrediction_1043SitesModel.R script", log = T)


# OPTIONS
options(scipen = 999)





# DATA INPUT ----------------------------------------------------------------------------------------------------------------------------------------

# # GLOBAL MODEL CALIBRATION VARIABLES
# load(here::here("data/calibration/globalModel", "calibration_globalModel.Rdata"))

# extract high and low price
p_high <- 42.03


# # clear unnecessary objects
# rm(matrixTransition.2prices, calibration.globalModel)


# 1043 SITES MODEL CALIBRATION VARIABLES
load(here::here("data/hmc", "hmc_1043SitesModel.Rdata"))




# 1043 SITES AGGREGATE PREDICTION
aux.prices <- c(7.6, 17.6, 22.6, 27.6,32.6)



# read the site-level predictions with high cattle price to calculate aggregated agricultural output
prediction.1043SitesModel <- purrr::map_df(.x = aux.prices,
                                           .f = ~readr::read_delim(here::here("data/prediction/1043-det",
                                                                              glue::glue("p_e_{.x}/amazon_data_z.dat")),
                                                                   col_types = cols(.default = "n")) %>%
                                             dplyr::mutate(p_e = .x, p_a = p_high)) 





# # read the site-level predictions with low cattle price to calculate aggregated agricultural output
# prediction.1043SitesModel <- purrr::map_df(.x = aux.prices,
#                                            .f = ~readr::read_delim(here::here("data/prediction/1043SitesModel",
#                                                                               glue::glue("pl_p_e_{.x}/amazon_data_z.dat")),
#                                                                    col_types = cols(.default = "n")) %>%
#                                              dplyr::mutate(p_e = .x, p_a = p_low)) %>%
#   bind_rows(prediction.1043SitesModel)


# AMAZON BIOME VECTOR DATA
load(here::here("data/raw2clean/amazonBiome_ibge/output/clean_amazonBiome.Rdata"))




# DATASET CLEANUP AND PREP ---------------------------------------------------------------------------------------------------------------------------

# select relevant calibrated variables
calibration.1043SitesModel <-
  calibration.1043SitesModel %>%
  dplyr::mutate(rank_theta_1043Sites = dense_rank(desc(theta_1043Sites)), # rank values
                rank_gamma_1043Sites = dense_rank(desc(gamma_1043Sites))) %>%  # rank values
  dplyr::select(id, z_2017_1043Sites, theta_1043Sites,  rank_theta_1043Sites, rank_gamma_1043Sites, x_2017_1043Sites, zbar_2017_1043Sites)

# adjust predicted data to site-year panel + add calibrated variables
prediction.1043SitesModel <-
  prediction.1043SitesModel %>%
  tidyr::pivot_longer(cols = -c("T/R ", "p_a", "p_e"), names_to = "id", values_to = "z_t") %>%
  dplyr::mutate(id = as.numeric(stringr::str_trim(id))) %>%
  dplyr::rename(time = "T/R ") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
  dplyr::right_join(calibration.1043SitesModel, by = c("id" = "id")) %>% # match by id, guarantee that prediction data uses the same id than the calibrated data
  sf::st_as_sf() %>%
  dplyr::mutate(z_t = 100*1e11*z_t/zbar_2017_1043Sites,
                z_2017_1043Sites = 100*z_2017_1043Sites/zbar_2017_1043Sites) # transform share to %


# clean environment
rm(calibration.1043SitesModel)

# adjust projection
prediction.1043SitesModel <- sf::st_as_sf(prediction.1043SitesModel)
prediction.1043SitesModel <- sf::st_transform(prediction.1043SitesModel, sf::st_crs(4326))
clean.amazonBiome <- sf::st_transform(clean.amazonBiome, sf::st_crs(prediction.1043SitesModel))





# GENERATE MAPS --------------------------------------------------------------------------------------------------------------------------------------

# PLOT z_2017, GAMMA, NAD THETA CALIBRATED VALUES


#### movie

####test 

aux.years <- seq(from = 2017, to = 2050, by = 1)

saveGIF({
  for (y in seq_along(aux.years)) {
plot_b15  <-
  ggplot2::ggplot(data = prediction.1043SitesModel %>%
                    dplyr::filter(time == aux.years[y]-2017, p_e == 22.6) %>%
                    dplyr::mutate(z_t = cut(round(z_t,5),
                                            breaks = c(0,0.5,20,40,60,80,105),
                                            include.lowest = T,
                                            dig.lab = 3,
                                            labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")))) +
  ggplot2::geom_sf(aes(fill = z_t)) +
  ggplot2::scale_fill_manual(
    name = bquote("Year" ~ .(aux.years[y]) ~ "Percent of land allocated to agriculture"),
    values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), 
    drop = FALSE
  ) +
  ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::labs(subtitle = bquote("(b) Transfer payment of $15/ton CO2 captured"))+
  ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                panel.background = element_blank(),
                strip.background = element_rect(fill = NA),
                axis.line = element_blank(), axis.ticks = element_blank(),
                axis.title = element_blank(), axis.text = element_blank(),
                legend.title = element_text(hjust = 0.5, size = 21, face = "bold"),
                legend.position = "bottom", legend.margin=margin(t=-1, r=-0, b=0.3, l=-0, unit="cm"),
                legend.text = element_text(size = 21),
                plot.margin = margin(t = 2, r = 10, b = 1.5, l = -2, unit = "cm"),
                  plot.subtitle = element_text(size = 21, face = "bold", hjust = 0.5))

# Additional code to display or save plot_b15
# ...



  plot_b0  <-
    ggplot2::ggplot(data = prediction.1043SitesModel %>%
                      dplyr::filter(time == aux.years[y]-2017, p_e == 7.6) %>%
                      dplyr::mutate(z_t = cut(round(z_t,5),
                                              breaks = c(0,0.5,20,40,60,80,105),
                                              include.lowest = T,
                                              dig.lab = 3,
                                              labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")))) +
    ggplot2::geom_sf(aes(fill = z_t)) +
  ggplot2::scale_fill_manual(
    name = bquote("Year" ~ .(aux.years[y]) ~ "Percent of land allocated to agriculture"),
    values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), 
    drop = FALSE
  ) +
    ggplot2::geom_sf(data = clean.amazonBiome, fill = NA, color = "darkgreen", size = 1.2) +
    ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
    ggplot2::labs(subtitle = bquote("(a) Business as usual")) +
    ggplot2::theme(panel.grid.major = element_line(colour = "white"),
                  panel.grid.minor = element_line(colour = "white"),
                  panel.background = element_blank(),
                  strip.background = element_rect(fill = NA),
                  axis.line = element_blank(), axis.ticks = element_blank(),
                  axis.title = element_blank(), axis.text = element_blank(),
                  legend.title = element_text(hjust = 0.5, size = 21, face = "bold"),
                  legend.position = "bottom", legend.margin=margin(t=-1, r=-0, b=0.3, l=-0, unit="cm"),
                  legend.text = element_text(size = 21),
                  plot.margin = margin(t = 2, r = -2, b = 1.5, l = 10, unit = "cm"),
                  plot.subtitle = element_text(size = 21, face = "bold", hjust = 0.5))

  # Additional code to display or save plot_b0
  # ...   

  #  combined_plot <- plot_grid(plot_b0, plot_b15, ncol = 2, align = "v", axis = "tblr", rel_widths = c(1, 1))
#     layout <- c(
#   area(t = 1, l = 1, b = 3, r = 3.4),
#   area(t = 1, l = 3.6, b = 3, r = 5)
# )
    combined_plot <- plot_b0 + plot_b15 + 
      plot_layout(ncol = 2,guides = "collect") & theme(legend.position='bottom')
  # combined_plot <- plot_b0 +plot_spacer()+ plot_b15 + plot_layout(ncol = 3,widths=c(1,-0.3,1),guides = "collect") & theme(legend.position='bottom')
    print(combined_plot)
  }
}, movie.name = "movie_1043site.gif", interval = 0.5, ani.width = 1650, ani.height = 650)


# END TIMER
tictoc::toc(log = T)

# export time to csv table
ExportTimeProcessing("code/analysis")






# END OF SCRIPT --------------------------------------------------------------------------------------------------------------------------------------