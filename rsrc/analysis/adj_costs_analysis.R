library(sf)
library(terra)
library(stargazer)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

# Load raster data (amazon biome share, pixel areas, and land uses)
amazon_rsts <- rast(
  list.files(
    "data/processed/",
    pattern = ".tif",
    full.names = TRUE
  )
)

# Load percipitation data
precip_rsts <- rast("data/clean/precipitation.tif")

# Resample percipitation data into sites
mean_precip_rst <- precip_rsts |>
  resample(amazon_rsts) |>
  mean()

names(mean_precip_rst) <- "mean_yearly_percipitation"

# Stack with rest of amazon rasters
amazon_rsts <- c(amazon_rsts, mean_precip_rst)

# Extract variables as polygons and project data
calib_df <- as.polygons(amazon_rsts, dissolve = FALSE) %>%
  st_as_sf() %>%
  st_transform(5880)

# Remove sites with less than 3% overlap with the amazon biome
calib_df <- calib_df %>%
  filter(share_amazon_biome >= 0.03)

# Add id variable
calib_df$id <- seq_len(nrow(calib_df))

# Transform share variables into area (ha)
calib_df <- calib_df %>%
  mutate(across(
    starts_with("share_"),
    .names = "area_{.col}",
    ~ . * pixel_area_ha
  )) %>%
  rename_with(
    ~ str_replace(., "share_", ""),
    starts_with("area_")
  )

# Check if areas with poor pasture have less secondary vegetation
# jose model

model_1 <- calib_df %>% lm(
  area_secondary_vegetation_2017 ~
    share_pasture_quality_1_2017 +
    share_pasture_quality_3_2017 +
    area_agricultural_use_2017 +
    mean_yearly_percipitation,
  data = .
)

stargazer(model_1, out = "plots/adj_costs_calib/table_1.tex")
stargazer(model_1, out = "plots/adj_costs_calib/table_1.txt")


# Levels
model_1 <- calib_df %>% lm(
  area_secondary_vegetation_2017 ~ area_pasture_quality_1_2017,
  data = .
)

model_2 <- calib_df %>%
  lm(
    area_secondary_vegetation_2017 ~
      area_pasture_quality_1_2017 +
      area_agricultural_use_2017,
    data = .
  )

# Logs
model_3 <- calib_df %>%
  filter(area_secondary_vegetation_2017 > 0) %>%
  filter(area_pasture_quality_1_2017 > 0) %>%
  lm(
    log(area_secondary_vegetation_2017) ~
      log(area_pasture_quality_1_2017),
    data = .
  )

model_4 <- calib_df %>%
  filter(area_secondary_vegetation_2017 > 0) %>%
  filter(area_pasture_quality_1_2017 > 0) %>%
  lm(
    log(area_secondary_vegetation_2017) ~
      log(area_pasture_quality_1_2017) +
      log(area_agricultural_use_2017),
    data = .
  )


stargazer(model_1, model_2, model_3, model_4, out = "plots/adj_costs_calib/table_2.tex")
stargazer(model_1, model_2, model_3, model_4, out = "plots/adj_costs_calib/table_2.txt")


# Now control for percipitation
# Levels
model_5 <- calib_df %>%
  filter(area_secondary_vegetation_2017 > 0) %>%
  lm(
    area_secondary_vegetation_2017 ~
      area_pasture_quality_1_2017 +
      mean_yearly_percipitation,
    data = .
  )

model_6 <- calib_df %>%
  filter(area_secondary_vegetation_2017 > 0) %>%
  lm(
    area_secondary_vegetation_2017 ~
      area_pasture_quality_1_2017 +
      area_agricultural_use_2017 +
      mean_yearly_percipitation,
    data = .
  )

# Logs
model_7 <- calib_df %>%
  filter(area_secondary_vegetation_2017 > 0) %>%
  filter(area_pasture_quality_1_2017 > 0) %>%
  lm(
    log(area_secondary_vegetation_2017) ~
      log(area_pasture_quality_1_2017) +
      mean_yearly_percipitation,
    data = .
  )

model_8 <- calib_df %>%
  filter(area_secondary_vegetation_2017 > 0) %>%
  filter(area_pasture_quality_1_2017 > 0) %>%
  lm(
    log(area_secondary_vegetation_2017) ~
      log(area_pasture_quality_1_2017) +
      log(area_agricultural_use_2017) +
      mean_yearly_percipitation,
    data = .
  )


stargazer(model_5, model_6, model_7, model_8, out = "plots/adj_costs_calib/table_3.tex")
stargazer(model_5, model_6, model_7, model_8, out = "plots/adj_costs_calib/table_3.txt")




# Now on lag
# Levels
model_9 <- calib_df %>%
  filter(area_secondary_vegetation_2017 > 0) %>%
  lm(
    area_secondary_vegetation_2017 ~
      area_pasture_quality_1_2011 +
      mean_yearly_percipitation,
    data = .
  )

model_10 <- calib_df %>%
  filter(area_secondary_vegetation_2017 > 0) %>%
  lm(
    area_secondary_vegetation_2017 ~
      area_pasture_quality_1_2011 +
      area_agricultural_use_2017 +
      mean_yearly_percipitation,
    data = .
  )

# Logs
model_11 <- calib_df %>%
  filter(area_secondary_vegetation_2017 > 0) %>%
  filter(area_pasture_quality_1_2011 > 0) %>%
  lm(
    log(area_secondary_vegetation_2017) ~
      log(area_pasture_quality_1_2011) +
      mean_yearly_percipitation,
    data = .
  )

model_12 <- calib_df %>%
  filter(area_secondary_vegetation_2017 > 0) %>%
  filter(area_pasture_quality_1_2011 > 0) %>%
  lm(
    log(area_secondary_vegetation_2017) ~
      log(area_pasture_quality_1_2011) +
      log(area_agricultural_use_2017) +
      mean_yearly_percipitation,
    data = .
  )


stargazer(model_9, model_10, model_11, model_12, out = "plots/adj_costs_calib/table_4.tex")
stargazer(model_9, model_10, model_11, model_12, out = "plots/adj_costs_calib/table_4.txt")



fig_1 <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = share_secondary_vegetation_2017 * 100), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    fill = "Secondary vegetation area(%)",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/adj_costs_calib/sec_veg_share.pdf", plot = fig_1)

fig_2 <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = share_pasture_quality_1_2017 * 100), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    fill = "Low pasture quality area (%)",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/adj_costs_calib/low_pasture_quality_share.pdf", plot = fig_2)


fig_3 <- calib_df %>%
  filter(area_secondary_vegetation_2017 > 0) %>%
  ggplot(aes(x = area_secondary_vegetation_2017, y = area_pasture_quality_1_2017)) +
  geom_point() +
  labs(
    x = "Secondary vegetation area (hectares)",
    y = "Low pasture quality area (hectares)"
  )

ggsave(filename = "plots/adj_costs_calib/scatter_sec_veg_low_pasture.pdf", plot = fig_3)


fig_4 <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = mean_yearly_percipitation), color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey") +
  labs(
    fill = "Average yearly percipitation (mm)",
    x = "Longitude",
    y = "Latitude"
  )

ggsave(filename = "plots/adj_costs_calib/avg_yearly_percipitation.pdf", plot = fig_4)


df <- calib_df %>%
  mutate(diff = share_pasture_quality_3_2017 - share_pasture_quality_1_2017) %>%
  select(
    diff,
    mean_yearly_percipitation
  ) %>%
  st_drop_geometry() %>%
  data.frame()

cor(df)


calib_df %>%
  mutate(
    total_pasture_area = area_pasture_quality_1_2017 +
      area_pasture_quality_2_2017 +
      area_pasture_quality_3_2017
  ) %>%
  ggplot(aes(x = total_pasture_area, y = area_agricultural_use_2017)) +
  geom_point() +
  labs(
    x = "total pasture area (hectares)",
    y = "Z (hectares)"
  )

ggsave(filename = "plots/adj_costs_calib/Z_vs_pasture_area.pdf", plot = fig_3)
