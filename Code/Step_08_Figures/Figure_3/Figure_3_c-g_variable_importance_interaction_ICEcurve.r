library(data.table)
library(dplyr)
library(randomForest)
library(vivid)
library(intergraph)
library(igraph)
library(raster)
library(ggplot2)
library(patchwork)

# -----------------------------
# 1. Read and prepare real data
# -----------------------------

setwd("D:/MaxTreeHeightProject")

rawData <- fread(
  "Data/Tree_height_covariates_extraction_merged_20260415_filtered_BIENdata_biome975.csv"
)[, -1] %>% 
  na.omit()

# Forest type classification
rawData <- rawData %>%
  mutate(
    ForestType = case_when(
      WWF_Biome %in% c(1, 2, 3, 7, 9, 14) ~ "Tropical",
      WWF_Biome %in% c(4, 5, 8, 10) ~ "Temperate",
      WWF_Biome %in% c(6, 11) ~ "Boreal",
      WWF_Biome %in% c(12, 13) ~ "Dryland",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(ForestType))


# -----------------------------
# 2. Define response and predictors
# -----------------------------

response <- "MaxHeight"

vars_new <- c(
  "ForestAge",
  "CHELSA_Temperature_Annual_Range",
  "SG_Soil_pH_H2O_0_100cm",
  "EarthEnvTopoMed_Slope",
  "Human_Disturbance",
  "WaterAvailability",
  "CHELSA_vpd",
  "WorldClim2_WindSpeed_AnnualMean"
)

myData <- rawData %>%
  dplyr::select(all_of(c(response, c(vars_new,"x","y")))) %>%
  na.omit()

# Optional: unit conversion, consistent with your SHAP plots
cleanedTrainTable <- myData %>%
  mutate(
    CHELSA_Temperature_Annual_Range = CHELSA_Temperature_Annual_Range / 10,
    SG_Soil_pH_H2O_0_100cm = SG_Soil_pH_H2O_0_100cm / 10,
    CHELSA_vpd = CHELSA_vpd / 10000
  )

# -----------------------------
# 3. Fit random forest model
# -----------------------------

source("D:/MaxTreeHeightProject/Rcode/sample.grid_from_GSIF_package.r")
duplicateTable = cleanedTrainTable
# tranform the data frame format lat lon into spatial lat long as spatial points
coordinates(cleanedTrainTable) = ~ x + y
# allocate the projection
proj4string(cleanedTrainTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
cleanedTrainTable@data = duplicateTable
# # set the 
set.seed(1701)
gridSubsampledPoints = sample.grid(cleanedTrainTable, cell.size = c(0.25,0.25), n = 1)
# print(dim(gridSubsampledPoints$subset@data))
gridSubsampledPoints = gridSubsampledPoints[[1]]
# get the grid subsample data frame
gridSubsampledTable = gridSubsampledPoints@data %>%
  dplyr::select(-x, -y)

set.seed(1701)
rf <- randomForest(
  MaxHeight ~ .,
  data = gridSubsampledTable,
  ntree = 250,
  mtry = 3,
  nodesize = 20,
  importance = TRUE
)

print(rf)
varImpPlot(rf)

# -----------------------------
# 4. Calculate variable importance and interaction matrix
# -----------------------------

set.seed(1701)
viviRf <- vivi(
  fit = rf,
  data = gridSubsampledTable,
  response = "MaxHeight",
  gridSize = 50,
  importanceType = "agnostic",
  nmax = 1000,
  reorder = TRUE,
  predictFun = NULL,
  numPerm = 4,
  showVimpError = FALSE
)

# Interaction network
viviNetwork(viviRf)

var_labels <- c(
  "Human_Disturbance" = "Human modification",
  "ForestAge" = "Forest age",
  "SG_Soil_pH_H2O_0_100cm" = "Soil pH",
  "CHELSA_Temperature_Annual_Range" = "Seasonal temperature variation",
  "EarthEnvTopoMed_Slope" = "Slope",
  "WaterAvailability" = "Water availability",
  "CHELSA_vpd" = "Vapour pressure deficit",
  "WorldClim2_WindSpeed_AnnualMean" = "Wind speed"
)

viviRf_plot <- viviRf

rownames(viviRf_plot) <- var_labels[rownames(viviRf_plot)]
colnames(viviRf_plot) <- var_labels[colnames(viviRf_plot)]
int_col <- colorRampPalette(c("#F2F4F7","#C6DBEF","#6BAED6","#2171B5"))(100)
#imp_col <- colorRampPalette(c("#F0F7F4", "#BFD8C2","#74A892","#2F6F5E"))(100)
imp_col <- rep("#C75B4E", 100)

source("Rcode/viviNetwork_plot_variable_interaction_network.R")


# Figure 3 panel c
pdf("Plots/Figure_03_vars_viviNetwork.pdf",width = 8,height = 6)
par(mar = c(0, 0, 0, 0))
set.seed(2000)

p<-viviNetwork_new(viviRf_plot,intPal = int_col,impPal = imp_col,edgeWidths = 1:4,nudge_x = 0.05,nudge_y = 0.03)
p +
  theme(
    text = element_text(size = 16),
    axis.text = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
dev.off()


# -----------------------------
# 5. Partial dependence plots for single variables

top4 <- colnames(viviRf)[1:4]
topvar<- c("ForestAge","WaterAvailability","CHELSA_Temperature_Annual_Range","SG_Soil_pH_H2O_0_100cm")
args(pdpVars)
# -----------------------------
set.seed(2000)
p_pdp <- pdpVars(
  data = gridSubsampledTable,
  fit = rf,
  response = "MaxHeight",
  vars = topvar,
  nIce = 100
)
str(p_pdp)

# extract data plot ----------
p_pdp_gg <- p_pdp[sapply(p_pdp, function(x) inherits(x, "ggplot"))]

length(p_pdp_gg)

pdp_data <- lapply(seq_along(p_pdp_gg), function(i) {
  dat <- p_pdp_gg[[i]]$data   # 推荐用 $data，不用 @data
  dat$var_plot <- topvar[i]
  dat
}) %>%
  bind_rows()

y_lim <- c(20,70)

plot_order <- c(
  "CHELSA_Temperature_Annual_Range",
  "WaterAvailability",
  "ForestAge",
  "SG_Soil_pH_H2O_0_100cm"
)


x_labels <- c(
  "CHELSA_Temperature_Annual_Range" = "Seasonal temperature variation",
  "WaterAvailability" = "Water availability",
  "ForestAge" = "Forest age",
  "SG_Soil_pH_H2O_0_100cm" = "Soil pH"
)

# Figure 3 panel d-f
panel_labels <- c("d", "e", "f", "g")
names(panel_labels) <- plot_order
make_pdp_plot <- function(var_name, show_y = TRUE, panel_label = NULL) {
  dat <- pdp_data %>% 
    filter(var_plot == var_name)
  pdp_mean <- dat %>%
    group_by(x = .data[[var_name]]) %>%
    summarise(fit_mean = mean(fit, na.rm = TRUE), .groups = "drop")
  
  ggplot(dat, aes(x = .data[[var_name]], y = fit)) +
    geom_line(aes(group = .id, colour = MaxHeight),linewidth = 0.25,alpha = 0.6) +
    geom_line(data = pdp_mean,aes(x = x, y = fit_mean),inherit.aes = FALSE,colour = "black",linewidth = 1.1) +
    annotate( "text",x = -Inf, y = Inf,label = panel_label,hjust = -0.6,vjust = 1.4,size = 10) +
    scale_colour_gradientn(
      colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
      name = "MaxHeight (m)"
    ) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    labs(
      x = x_labels[var_name],
      y = if (show_y) "Maximum tree height (m)" else NULL
    ) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      axis.text.y = if (show_y) element_text(colour = "black") else element_blank(),
      axis.ticks.y = if (show_y) element_line(colour = "black") else element_blank(),
      legend.position = "right",
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(t = 6, r = 10, b = 6, l = 10)
    )
}
p_list <- lapply(seq_along(plot_order), function(i) {
  make_pdp_plot(
    var_name = plot_order[i],
    show_y = i == 1,
    panel_label = panel_labels[plot_order[i]]
  )
})

p_final <- wrap_plots(p_list, ncol = 4, guides = "collect") & 
  theme(legend.position = "none", legend.title = element_text(size = 14), legend.text = element_text(size = 12))
p_final

ggsave(filename = "Plots/Figure_03_PDP_ICE_final.pdf",plot = p_final,width = 14,height = 4,units = "in",dpi = 300)


