library(data.table)
library(dplyr)
library(ggplot2)
library(xgboost)
library(patchwork)
library(scales)

rawData <- fread("D:/MaxTreeHeightProject/Data/Tree_height_covariates_extraction_merged_20260407_filtered_BIENdata_biome975.csv")[,-1]

y_var <- "MaxHeight"
vars_new <- c("ForestAge", "CHELSA_Temperature_Annual_Range",
              "SG_Soil_pH_H2O_0_100cm", "EarthEnvTopoMed_Slope",
              "Human_Disturbance", "WaterAvailability",
              "CHELSA_vpd", "WorldClim2_WindSpeed_AnnualMean")

dat <- rawData[, c(y_var, vars_new), with = FALSE] %>% na.omit()
dat <- dat %>%
  mutate(
    CHELSA_Temperature_Annual_Range = CHELSA_Temperature_Annual_Range / 10,
    SG_Soil_pH_H2O_0_100cm = SG_Soil_pH_H2O_0_100cm / 10,
    CHELSA_vpd = CHELSA_vpd / 10000
  )
y <- dat[[y_var]]
X <- dat[, vars_new, with = FALSE]
X_mat <- as.matrix(X)

set.seed(123)
train_id <- sample(seq_len(nrow(dat)), size = 0.8 * nrow(dat))

X_train <- X_mat[train_id,]
X_test  <- X_mat[-train_id,]
y_train <- y[train_id]
y_test  <- y[-train_id]

dtrain <- xgb.DMatrix(data = X_train, label = y_train)
dtest  <- xgb.DMatrix(data = X_test, label = y_test)

params <- list(objective = "reg:squarederror",eval_metric = "rmse",eta = 0.05,
               max_depth = 4,min_child_weight = 5,subsample = 0.8,colsample_bytree = 0.8)

set.seed(123)
xgb_fit <- xgb.train(params = params,data = dtrain,nrounds = 500,watchlist = list(train = dtrain, test = dtest),
                     early_stopping_rounds = 30,verbose = 1)

# -----------------------------
# 1. Prepare SHAP data
# -----------------------------
shap_mat <- predict(xgb_fit, X_test, predcontrib = TRUE)

shap_df <- as.data.frame(shap_mat[, -ncol(shap_mat)])
colnames(shap_df) <- colnames(X_test)

x_df <- as.data.frame(X_test)
tree_height_test <- y_test

plot_df <- lapply(colnames(x_df), function(v) {
  data.frame(variable = v,value = x_df[[v]],shap = shap_df[[v]],tree_height = tree_height_test)}) %>%
  bind_rows() %>%
  filter(is.finite(value), is.finite(shap))

# -----------------------------
# 2. Fit LOESS curves and identify SHAP = 0 thresholds
# -----------------------------
get_loess_threshold <- function(df, span = 0.75, n = 300) {
  
  df <- df %>% arrange(value)
  
  if (nrow(df) < 20 || length(unique(df$value)) < 10) {
    return(data.frame())
  }
  
  fit <- loess(shap ~ value, data = df, span = span, degree = 2)
  newx <- seq(min(df$value, na.rm = TRUE), max(df$value, na.rm = TRUE),
  length.out = n)
  pred <- predict(fit, newdata = data.frame(value = newx), se = TRUE)
  smooth_df <- data.frame(
    variable = unique(df$variable),value = newx,
    fit = pred$fit,se = pred$se.fit
  ) %>% 
    mutate(ymin = fit - 1.96 * se,ymax = fit + 1.96 * se,sign = ifelse(fit >= 0, "Positive", "Negative"))
  
  # Identify zero-crossing points where the LOESS curve intersects SHAP = 0
  idx <- which(diff(sign(smooth_df$fit)) != 0)
  threshold_df <- data.frame()
  
  if (length(idx) > 0) {
    threshold_df <- lapply(idx, function(i) {
      x1 <- smooth_df$value[i]
      x2 <- smooth_df$value[i + 1]
      y1 <- smooth_df$fit[i]
      y2 <- smooth_df$fit[i + 1]
      
      x0 <- x1 - y1 * (x2 - x1) / (y2 - y1)
      
      data.frame(variable = unique(df$variable),threshold = x0)
    }) %>%
      bind_rows()
  }
  list(smooth = smooth_df, threshold = threshold_df)
}

res_list <- plot_df %>%
  group_split(variable) %>%
  lapply(get_loess_threshold)

smooth_df <- bind_rows(lapply(res_list, `[[`, "smooth"))
threshold_df <- bind_rows(lapply(res_list, `[[`, "threshold"))

# -----------------------------
# 3. Define a plotting function for individual SHAP dependence plots
# -----------------------------
# Define variable labels, plotting order, and common y-axis range
var_labels <- c(
  "Human_Disturbance" = "Human modification",
  "ForestAge" = "Forest age (yr)",
  "SG_Soil_pH_H2O_0_100cm" = "Soil pH",
  "CHELSA_Temperature_Annual_Range" = "Seasonal temperature variation (°C)",
  "EarthEnvTopoMed_Slope" = "Slope (°)",
  "WaterAvailability" = "Water availability (mm)",
  "CHELSA_vpd" = "Vapour pressure deficit (kPa)",
  "WorldClim2_WindSpeed_AnnualMean" = "Wind speed (m/s)")

my_order <- c(
  "CHELSA_Temperature_Annual_Range",
  "WaterAvailability",
  "ForestAge",
  "SG_Soil_pH_H2O_0_100cm",
  "WorldClim2_WindSpeed_AnnualMean",
  "EarthEnvTopoMed_Slope",
  "CHELSA_vpd",
  "Human_Disturbance"
)

# Set a common y-axis range for all panels
y_lim <- range(plot_df$shap, na.rm = TRUE)
y_pad <- diff(y_lim) * 0.08
y_lim <- c(y_lim[1] - y_pad, y_lim[2] + y_pad)

# Use a consistent position for threshold labels across all panels
label_df <- threshold_df %>%
  mutate(
    y_lab = y_lim[2] - diff(y_lim) * 0.08,
    label = round(threshold, 1))

# -----------------------------
# 4.  SHAP dependence plot
# -----------------------------

plot_one_shap <- function(v, panel_lab, show_legend = FALSE) {
  
  df <- plot_df %>% filter(variable == v)
  
  sm <- smooth_df %>%
    filter(variable == v) %>%
    arrange(value) %>%
    mutate(xmax = lead(value)) %>%
    filter(!is.na(xmax))
  
  th <- threshold_df %>% filter(variable == v)
  lab <- label_df %>% filter(variable == v)
  
  ggplot(df, aes(x = value, y = shap)) +
    geom_rect(data = sm %>% filter(sign == "Negative"),aes(xmin = value, xmax = xmax, ymin = -Inf, ymax = 0),
              inherit.aes = FALSE,fill = "#DCEAF7",alpha = 0.55) +       # "#d9d2ea" "#E7F1EC"
    geom_rect(data = sm %>% filter(sign == "Positive"),aes(xmin = value, xmax = xmax, ymin = 0, ymax = Inf),
              inherit.aes = FALSE,fill = "#F8DDDD",alpha = 0.55) +       # "#e6e8b8" "#F6E7E5"
    geom_hline(yintercept = 0, linetype = "dashed",color = "grey65", linewidth = 0.5) +
    geom_point(aes(color = tree_height), alpha = 0.58, size = 1.35) +
    geom_ribbon(data = smooth_df %>% filter(variable == v),aes(x = value, ymin = ymin, ymax = ymax),inherit.aes = FALSE,
                fill = "#80B9B8", alpha = 0.25) +  # "#8e79c6"
    geom_line(
      data = smooth_df %>% filter(variable == v),aes(x = value, y = fit),
      inherit.aes = FALSE, color = "#157A7F",linewidth = 1.05) +         # "#6b4fb3"
    geom_vline(data = th,aes(xintercept = threshold),color = "#4A90C2",  # "#d84b4b"  "#E67E22"
               linetype = "longdash",linewidth = 0.65) +
    geom_text(data = lab,aes(x = threshold, y = y_lab, label = label),
              inherit.aes = FALSE,color = "#4A90C2", size = 6,vjust = 0.2,hjust = -0.2) +     # "#E67E22"
    scale_color_gradientn(colors = c("#3B2C73", "#2B8CBE", "#41AB5D", "#FDE725"),name = "MaxHeight (m)") +  # low = "#7452a3", high = "#b6bd45"
    annotate( "text", x = -Inf,y = Inf,label = panel_lab, hjust = -0.45,vjust = 1.35, size = 12 ) +
    guides(color = guide_colorbar(direction = "horizontal",title.position = "left",barwidth = unit(2.5, "cm"),
                                  barheight = unit(0.25, "cm"))) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    labs(x = var_labels[v],y = "SHAP value") +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 16, color = "black"),
      legend.position = if (show_legend) c(0.62, 0.10) else "none",
      legend.background = element_rect(fill = alpha("white", 0.25), color = NA),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      panel.border = element_rect(color = "grey30", linewidth = 0.6),
      plot.margin = ggplot2::margin(4, 10, 4, 10)
    )
}

panel_labs <- letters[1:length(my_order)]

p_list <- Map(function(v, lab, i) {
    plot_one_shap(v = v,panel_lab = lab,show_legend = i == length(my_order))
    },my_order,panel_labs,seq_along(my_order))


# delete right panel y axis
right_ids <- c(2:4, 6:8)

p_list[right_ids] <- lapply(p_list[right_ids],
                            function(p) {p + theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())})
p_all <- wrap_plots(p_list, ncol = 4)
p_all

ggsave(filename = "Plots/Figure_04_SHAP_dependence_MaxHeight.pdf",plot = p_all,width = 16,height = 8,units = "in")
ggsave(filename = "Plots/Figure_04_SHAP_dependence_MaxHeight.jpeg",plot = p_all,width = 16,height = 8,units = "in",dpi = 300)
