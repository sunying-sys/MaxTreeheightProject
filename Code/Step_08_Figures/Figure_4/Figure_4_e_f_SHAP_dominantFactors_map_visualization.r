library(terra)
library(dplyr)
library(ggplot2)
library(sf)
library(scales)
library(grid)
library(cowplot)

terraOptions(memfrac = 0.85, progress = 1)

# 0. Projection
equalEarthProj <- "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# 1. Read dominant variable raster
dominance_layer <- rast("D:/MaxTreeHeightProject/ProcessData/SHAPresult_map/var8/Dominant_variable_SHAP.tif")

predictor_files <- c(
  WaterAvailability = "D:/MaxTreeHeightProject/ProcessData/SHAPresult_map/var8/SHAP_WaterAvailability.tif",
  ForestAge = "D:/MaxTreeHeightProject/ProcessData/SHAPresult_map/var8/SHAP_ForestAge.tif",
  TemperatureRange = "D:/MaxTreeHeightProject/ProcessData/SHAPresult_map/var8/SHAP_CHELSA_Temperature_Annual_Range.tif",
  SoilpH = "D:/MaxTreeHeightProject/ProcessData/SHAPresult_map/var8/SHAP_SG_Soil_pH_H2O_0_100cm.tif",
  VPD = "D:/MaxTreeHeightProject/ProcessData/SHAPresult_map/var8/SHAP_CHELSA_vpd.tif",
  Slope = "D:/MaxTreeHeightProject/ProcessData/SHAPresult_map/var8/SHAP_EarthEnvTopoMed_Slope.tif",
  HumanDisturbance = "D:/MaxTreeHeightProject/ProcessData/SHAPresult_map/var8/SHAP_Human_Disturbance.tif",
  WindSpeed = "D:/MaxTreeHeightProject/ProcessData/SHAPresult_map/var8/SHAP_WorldClim2_WindSpeed_AnnualMean.tif"
)

shap_stack <- rast(predictor_files)
names(shap_stack) <- names(predictor_files)

# 2. Check alignment; resample only if needed
is_aligned <- compareGeom(dominance_layer, shap_stack[[1]], stopOnError = FALSE)

if (!isTRUE(is_aligned)) {
  dominance_layer <- resample(dominance_layer, shap_stack[[1]], method = "near")
}

# 3. Extract SHAP value of dominant factor
dominant_shap <- app(c(dominance_layer, shap_stack), fun = function(x) {
  id <- x[1]
  if (is.na(id) || id < 1 || id > 8) return(NA)
  x[id + 1]
})

dominant_strength <- abs(dominant_shap)
dominant_direction <- ifel(dominant_shap > 0, 1, 0)

names(dominance_layer) <- "dominant_id"
names(dominant_strength) <- "strength"
names(dominant_direction) <- "positive"

# 4. Crop to map extent
dominance_layer <- crop(dominance_layer, ext(-180, 180, -60, 85))
dominant_strength <- crop(dominant_strength, ext(-180, 180, -60, 85))
dominant_direction <- crop(dominant_direction, ext(-180, 180, -60, 85))

# 5. No aggregation: use original resolution
plot_stack <- c(dominance_layer, dominant_strength, dominant_direction)

# 6. Project to Equal Earth
template_equal <- project(plot_stack[[1]], equalEarthProj, method = "near")

dominance_equal <- project(dominance_layer, template_equal, method = "near")
strength_equal <- project(dominant_strength, template_equal, method = "bilinear")
direction_equal <- project(dominant_direction, template_equal, method = "near")

names(dominance_equal) <- "dominant_id"
names(strength_equal) <- "strength"
names(direction_equal) <- "positive"

plot_stack_equal <- c(dominance_equal, strength_equal, direction_equal)

# 7. Continent boundary
continent_vect <- terra::vect("D:/MaxTreeHeightProject/Data/world_continent_shp/continent.shp")
continent_vect <- crop(continent_vect, ext(-180, 180, -60, 85))
continent_vect <- project(continent_vect, equalEarthProj)
continent_poly <- st_as_sf(continent_vect)

# 8. Mask by land
plot_stack_equal_masked <- mask(plot_stack_equal, continent_vect)

# 9. Convert to dataframe
df <- as.data.frame(plot_stack_equal_masked, xy = TRUE, na.rm = TRUE)

factor_names <- c(
  "Water availability",
  "Forest age",
  "Seasonal temperature variation",
  "Soil pH",
  "Vapour pressure deficit",
  "Slope",
  "Human modification",
  "Wind speed"
)

legend_order <- c(
  "Water availability",
  "Forest age",
  "Seasonal temperature variation",
  "Soil pH",
  "Vapour pressure deficit",
  "Slope",
  "Human modification",
  "Wind speed"
)

df <- df |>
  filter(!is.na(dominant_id), !is.na(strength), dominant_id %in% 1:8) |>
  mutate(
    dominant_factor = factor(factor_names[dominant_id], levels = legend_order),
    direction_label = factor(
      ifelse(positive == 1, "Positive effect", "Negative effect"),
      levels = c("Positive effect", "Negative effect")
    ),
    strength_scaled = rescale(
      pmin(strength, quantile(strength, 0.99, na.rm = TRUE)),
      to = c(0.20, 0.75)
    )
  )

# 10. Colors
factor_cols <- c(
  "Water availability" = "#5B8CC0",
  "Forest age" = "#77B07A",
  "Seasonal temperature variation" = "#D4A35A",
  "Soil pH" = "#D78A8A",
  "Vapour pressure deficit" = "#6CB5B0",
  "Slope" = "#9B85BD",
  "Human modification" = "#C7B15A",
  "Wind speed" = "#A9C7DD"
)

direction_cols <- c(
  "Positive effect" = "#D98C8C",
  "Negative effect" = "#8DA8C9"
)

# 11. Dominant factor map
p_main <- ggplot() +
  geom_sf(data = continent_poly, fill = "#F3F1EC", color = NA) +
  geom_raster(
    data = df,
    aes(x = x, y = y, fill = dominant_factor, alpha = strength_scaled)
  ) +
  geom_sf(data = continent_poly, fill = NA, color = "grey75", linewidth = 0.08) +
  scale_fill_manual(values = factor_cols, name = NULL, na.translate = FALSE) +
  scale_alpha(range = c(0.20, 0.75), guide = "none") +
  coord_sf(crs = equalEarthProj, datum = NA, expand = FALSE) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2))

legend_factor_df <- data.frame(
  dominant_factor = factor(legend_order, levels = legend_order),
  x = c(1, 2.4, 3.8, 5.4, 1, 2.4, 3.8, 5.4),
  y = c(2, 2, 2, 2, 1, 1, 1, 1),
  label = legend_order
)

p_legend_factor <- ggplot(legend_factor_df) +
  geom_tile(aes(x = x, y = y, fill = dominant_factor), width = 0.18, height = 0.25) +
  geom_text(aes(x = x + 0.20, y = y, label = label), hjust = 0, size = 3.1) +
  scale_fill_manual(values = factor_cols, guide = "none") +
  annotate("text", x = 0.65, y = 2.65, label = "Dominant factor", hjust = 0, fontface = "bold", size = 3.4) +
  xlim(0.6, 6.8) +
  ylim(0.45, 2.85) +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, color = "grey85", linewidth = 0.35),
        plot.margin = margin(1, 1, 1, 1))

legend_strength_df <- data.frame(
  x = seq(1, 8, length.out = 8),
  y = 1,
  strength = seq(0.20, 0.75, length.out = 8)
)

p_legend_strength <- ggplot(legend_strength_df) +
  geom_tile(aes(x = x, y = y, alpha = strength), fill = "black", width = 0.9, height = 0.32, color = NA) +
  scale_alpha(range = c(0.20, 0.75), guide = "none") +
  annotate("text", x = 1, y = 1.65, label = "SHAP strength", hjust = 0, fontface = "bold", size = 3.4) +
  annotate("text", x = 1, y = 0.45, label = "Weak", hjust = 0, size = 3.1) +
  annotate("text", x = 8, y = 0.45, label = "Strong", hjust = 1, size = 3.1) +
  xlim(0.8, 8.2) +
  ylim(0.25, 1.85) +
  theme_void()

p_legend <- ggdraw() +
  draw_plot(p_legend_factor, x = 0.02, y = 0.05, width = 0.70, height = 0.90) +
  draw_plot(p_legend_strength, x = 0.77, y = 0.13, width = 0.20, height = 0.74)

p_final <- plot_grid(
  p_main,
  p_legend,
  ncol = 1,
  rel_heights = c(1, 0.22)
)

ggsave(
  "D:/MaxTreeHeightProject/Plots/Global_SHAP_dominant_factor_strength_EqualEarth.pdf",
  p_final, width = 12, height = 7
)

# 12. Direction map
p_dir <- ggplot() +
  geom_sf(data = continent_poly, fill = "#F3F1EC", color = NA) +
  geom_raster(
    data = df,
    aes(x = x, y = y, fill = direction_label, alpha = strength_scaled)
  ) +
  geom_sf(data = continent_poly, fill = NA, color = "grey75", linewidth = 0.08) +
  scale_fill_manual(values = direction_cols, name = NULL, na.translate = FALSE) +
  scale_alpha(range = c(0.20, 0.75), guide = "none") +
  coord_sf(crs = equalEarthProj, datum = NA, expand = FALSE) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2))

legend_df <- expand.grid(
  strength = seq(0, 1, length.out = 100),
  direction = c("Negative effect", "Positive effect")
) |>
  mutate(direction = factor(direction, levels = c("Negative effect", "Positive effect")))

p_dir_legend <- ggplot(legend_df, aes(x = strength, y = direction)) +
  geom_tile(aes(fill = direction, alpha = strength)) +
  scale_fill_manual(values = direction_cols, guide = "none") +
  scale_alpha(range = c(0.20, 0.75), guide = "none") +
  scale_x_continuous(breaks = c(0, 1), labels = c("Weak", "Strong"), expand = c(0, 0)) +
  scale_y_discrete(labels = c("Negative", "Positive"), expand = c(0, 0)) +
  labs(x = "SHAP strength", y = NULL) +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    plot.margin = margin(1, 1, 1, 1),
    panel.border = element_rect(fill = NA, color = "grey40", linewidth = 0.3)
  )

p_dir_final <- ggdraw(p_dir) +
  draw_plot(p_dir_legend, x = 0.06, y = 0.12, width = 0.14, height = 0.14)

ggsave(
  "D:/MaxTreeHeightProject/Plots/Global_SHAP_dominant_direction_strength_EqualEarth.pdf",
  p_dir_final, width = 12, height = 7
)