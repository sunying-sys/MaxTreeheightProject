library(ggplot2)
library(dplyr)
library(viridis)
library(data.table)
library(patchwork)
library(tibble)

setwd("D:/MaxTreeHeightProject")


tableAggred1 = fread("D:/MaxTreeHeightProject/RF_importance/var_importance_BiomeQuantile975.csv")
tableAggred2 = fread("D:/MaxTreeHeightProject/PartialResults/Partial_Regression_coef_BiomeQuantile975.csv")
colnames(tableAggred1)
colnames(tableAggred2)[colnames(tableAggred2) == "Variables"] <- "variable"

df_merge <- tableAggred2 %>%
  left_join(tableAggred1, by = "variable")

var_info <- tribble(
  ~variable, ~varType,
  "CHELSA_vpd",                      "Climate",
  "ForestAge",                       "Vegetation",
  "CHELSA_Temperature_Annual_Range", "Climate",
  "SG_Soil_pH_H2O_0_100cm",          "Soil",
  "EarthEnvTopoMed_Slope",           "Topographty",
  "Human_Disturbance",               "Human",
  "WaterAvailability",                 "Climate",
  "WorldClim2_WindSpeed_AnnualMean", "Climate"
)
tableAggred <- df_merge %>%
  left_join(var_info, by = "variable")

tableAggred <- as.data.table(tableAggred)

my_order <- c("CHELSA_Temperature_Annual_Range","WaterAvailability","ForestAge","SG_Soil_pH_H2O_0_100cm",
              "WorldClim2_WindSpeed_AnnualMean","EarthEnvTopoMed_Slope","CHELSA_vpd","Human_Disturbance")
tableAggred[, variable := factor(variable, levels = rev(my_order))]

var_labs <- c(
  "Human_Disturbance" = "Human modification",
  "ForestAge" = "Forest age",
  "SG_Soil_pH_H2O_0_100cm" = "Soil pH",
  "CHELSA_Temperature_Annual_Range" = "Seasonal temperature variation",
  "EarthEnvTopoMed_Slope" = "Slope",
  "WaterAvailability" = "Water availability",
  "CHELSA_vpd" = "Vapour pressure deficit",
  "WorldClim2_WindSpeed_AnnualMean" = "Wind speed"
)


p1 <- ggplot(tableAggred, aes(x = variable, y = percentage, fill = percentage)) +
  geom_col(width = 0.72, colour = "grey25", linewidth = 0.2) +
  coord_flip() +
  scale_x_discrete(labels = var_labs) +
  scale_y_continuous(limits = c(0,0.4),expand = c(0, 0)) +
  # scale_fill_gradient(
  #   low  = "#DCEFEA",
  #   high = "#1F6F78",
  #   guide = "none"
  # ) +
  scale_fill_gradient(
    low  = "#E5EEF5",
    high = "#2B7BBC",
    guide = "none"
  ) +
  labs(x = NULL, y = "Importance (Percentage)") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, colour = "black"),
    legend.position = "none"
  )
p1

p2 <- ggplot(tableAggred, aes(x = variable, y = Fit)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper),
                  colour = "#2B7BBC", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, colour = "#B04A4A", linewidth = 0.5) +
  coord_flip() +
  scale_x_discrete(labels = var_labs) +
  labs(x = NULL, y = "Partial regression coefficients") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text = element_text(size = 12, colour = "black")
  )

combined_plot <- p1 + p2 +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 18),
    plot.tag.position = c(0.05, 0.97)
  )

combined_plot

ggsave(filename = "Plots/Figure_03_Varimportance_Partial_all.pdf",plot = combined_plot,
  width = 20, height = 12, units = "cm")
ggsave(filename = "Plots/Figure_03_Varimportance_Partial_all.jpeg",plot = combined_plot,dpi=300,
       width = 20, height = 12, units = "cm")
