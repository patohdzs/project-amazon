library(stargazer)
library(ggplot2)
library(terra)
library(tidyverse)
library(sf)
library(conflicted)
library(units)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)


# Load 1043 calibration data set
load("data/calibration/gamma_calibration_1043_sites.Rdata")
high_res_df <- calib_df

load("data/calibration/gamma_calibration_78_sites.Rdata")

# Validate and clean geometries
high_res_df <- st_make_valid(high_res_df)
calib_df <- st_make_valid(calib_df)

df <- high_res_df %>%
  select(gamma_1043_site_reg = site_reg_gamma) %>%
  st_join(calib_df) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(gamma_1043_site_reg = mean(gamma_1043_site_reg, na.rm = TRUE))

# Add 1043 site-level CO2e to spatial variables
calib_df <- left_join(calib_df, df) %>%
  rename(gamma_78_site_reg = site_reg_gamma)


fig_1 <- calib_df %>%
  ggplot(aes(x = gamma_78_site_reg, y = gamma_1043_site_reg)) +
  geom_point() +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    x = "Gamma estimates from 78 sites regression",
    y = "Gamma estimates from 1043 sites regression, aggregated to 78 sites"
  )

ggsave(filename = "plots/gamma_calib/scatter_site_regs.pdf", plot = fig_1)


model_1 <- lm(gamma_1043_site_reg ~ 0 + gamma_78_site_reg, data = calib_df)
model_2 <- lm(gamma_1043_site_reg ~ gamma_78_site_reg, data = calib_df)

# Compute SSE (Sum of Squared Errors)
SSE <- sum((calib_df$gamma_1043_site_reg - predict(model_1))^2)

# Compute TSS (Total Sum of Squares) relative to the mean of y
TSS <- sum((calib_df$gamma_1043_site_reg - mean(calib_df$gamma_1043_site_reg))^2)

# Compute R^2 for model w/o intercept
R2 <- 1 - SSE / TSS

stargazer(model_1, model_2, out = "plots/gamma_calib/site_regs_comparison.tex")
stargazer(model_1, model_2, out = "plots/gamma_calib/site_regs_comparison.txt")
