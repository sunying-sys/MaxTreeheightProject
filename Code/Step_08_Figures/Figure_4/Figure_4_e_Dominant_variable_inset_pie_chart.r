library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("D:/MaxTreeHeightProject")

df <- fread("ProcessData/SHAPresult_map/var8/Dominant_variable_area_stat_by_pixels.csv")

factor_cols <- c(
  "Water availability" = "#3B6FB6",
  "Forest age" = "#5AA469",
  "Seasonal temperature variation" = "#E39C37",
  "Soil pH" = "#D65F5F",
  "Vapour pressure deficit" = "#4FA7A0",
  "Slope" = "#8C6BB1",
  "Human modification" = "#C9B037",
  "Wind speed" = "#a6cee3")

df <- df %>%
  mutate(
    dominant_var = factor(dominant_var, levels = names(factor_cols)),
    label = paste0(round(percent, 0), "%"),
    ymax = cumsum(proportion),
    ymin = lag(ymax, default = 0),
    ymid = (ymax + ymin) / 2
  )

p_pie <- ggplot(df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5, fill = dominant_var)) +
  geom_rect(color = "white", linewidth = 0.6) +
  coord_polar(theta = "y") +
  geom_text_repel(
    aes(x = 4, y = ymid, label = label),
    inherit.aes = FALSE,
    size = 5,
    segment.color = "grey40",
    segment.linewidth = 0.4,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = factor_cols) +
  xlim(0, 5.2) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 30, 10, 30)
  )

p_pie

ggsave(
  "Plots/Figure_04_02_Dominant_variable_pie_chart.pdf",
  p_pie,
  width = 7,
  height = 7
)
