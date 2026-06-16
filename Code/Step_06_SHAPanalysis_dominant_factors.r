# =========================================================
# Calculate pixel-wise SHAP values and derive:
# 1) Dominant variable map
# 2) Dominance strength map
# 3) Area statistics of dominant variables
# =========================================================

# =========================
# 1. Load required packages
# =========================
library(data.table)
library(dplyr)
library(terra)
library(xgboost)
library(sp)
library(raster)
# =========================
# 2. Set working directory
# =========================
setwd("~/Desktop/ YingSun/MaxTreeHeightProject")
source("Code/sample.grid_from_GSIF_package.r")
# =========================
# 3. Read training data
# =========================
train_df_raw <- fread("Data/Tree_height_covariates_extraction_merged_20260407_filtered_BIENdata_biome975.csv")

cleanedTrainTable <- train_df_raw %>%
  dplyr::select(
    MaxHeight,
    ForestAge,
    CHELSA_Temperature_Annual_Range,
    SG_Soil_pH_H2O_0_100cm,
    WaterAvailability,
    CHELSA_vpd,
    EarthEnvTopoMed_Slope,
    Human_Disturbance,
    WorldClim2_WindSpeed_AnnualMean,
    x, y
  ) %>%
  na.omit()

duplicateTable <- cleanedTrainTable

coordinates(cleanedTrainTable) <- ~ x + y
proj4string(cleanedTrainTable) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
cleanedTrainTable@data <- duplicateTable

set.seed(123)
gridSubsampledPoints <- sample.grid(cleanedTrainTable, cell.size = c(0.25, 0.25), n = 1)
gridSubsampledPoints <- gridSubsampledPoints[[1]]
gridSubsampledTable <- gridSubsampledPoints@data

# Define response variable
response_var <- "MaxHeight"

# Define predictor variables used in the SHAP analysis
predictor_vars <- c(
  "WaterAvailability",
  "ForestAge",
  "CHELSA_Temperature_Annual_Range",
  "SG_Soil_pH_H2O_0_100cm",
  "CHELSA_vpd",
  "EarthEnvTopoMed_Slope",
  "Human_Disturbance",
  "WorldClim2_WindSpeed_AnnualMean"
)

# Keep only response + predictor columns and remove missing values
train_df <- gridSubsampledTable %>%
  dplyr::select(all_of(c(response_var, predictor_vars))) %>%
  na.omit()

# ---------------------------------------------------------
# If unit conversion is needed, uncomment and make sure the
# same transformation is also applied to raster layers below
# ---------------------------------------------------------
train_df <- train_df %>%
  mutate(
    CHELSA_Temperature_Annual_Range = CHELSA_Temperature_Annual_Range / 10,
    SG_Soil_pH_H2O_0_100cm = SG_Soil_pH_H2O_0_100cm / 10,
    CHELSA_vpd = CHELSA_vpd / 10000
  )

# Create predictor matrix and response vector
X_train <- as.data.frame(train_df[, ..predictor_vars])
y_train <- train_df[[response_var]]

# =========================
# 4. Fit XGBoost model
# =========================
set.seed(123)

dtrain <- xgb.DMatrix(
  data = as.matrix(X_train),
  label = y_train
)

# xgboost ver 3.2.2.1
params <- xgb.params(
  max_depth = 6,
  learning_rate = 0.05,
  subsample = 0.8,
  colsample_bytree = 0.8,
  objective = "reg:squarederror",
  eval_metric = "rmse"
)

xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 300,
  verbose = 1
)


# =========================
# 5. Read raster predictors
# =========================
predictor_files <- c(
  WaterAvailability = "Data/layer/water_availability_CHELSE_merged.tif",
  ForestAge = "Data/layer/forestAge_merged.tif",
  CHELSA_Temperature_Annual_Range = "Data/layer/seasonTempVariation_CHELSE_merged.tif",
  SG_Soil_pH_H2O_0_100cm = "Data/layer/SG_Soil_pH_H2O_0_100cm_merged.tif",
  CHELSA_vpd = "Data/layer/CHELSA_vpd.tif",
  EarthEnvTopoMed_Slope = "Data/layer/EarthEnvSlope_merged.tif",
  Human_Disturbance = "Data/layer/humanModification_merged.tif",
  WorldClim2_WindSpeed_AnnualMean = "Data/layer/windSpeed_merged.tif"
)

r_list <- lapply(predictor_files, rast)
pred_stack0 <- rast(r_list)
names(pred_stack0) <- names(predictor_files)

# ---------------------------------------------------------
# If unit conversion is needed, uncomment and make sure the
# same transformation is also applied to training data above
# ---------------------------------------------------------
pred_stack0[["CHELSA_Temperature_Annual_Range"]] <- pred_stack0[["CHELSA_Temperature_Annual_Range"]] / 10
pred_stack0[["SG_Soil_pH_H2O_0_100cm"]] <- pred_stack0[["SG_Soil_pH_H2O_0_100cm"]] / 10
pred_stack0[["CHELSA_vpd"]] <- pred_stack0[["CHELSA_vpd"]] / 10000

# =========================
# 6. Apply forest mask
# =========================
presentCover <- rast("Data/layer/presentTreeCover_fromComposite_merged.tif")

forest_mask <- presentCover
forest_mask[forest_mask == 0] <- NA
forest_mask[forest_mask > 0] <- 1

pred_stack <- mask(pred_stack0, forest_mask)

# =========================
# 7. Optional: check consistency
# =========================
print(summary(train_df[, ..predictor_vars]))
print(global(pred_stack, fun = "range", na.rm = TRUE))

# =========================
# 8. Convert raster stack to data frame
# =========================
pred_df <- as.data.frame(pred_stack, xy = TRUE, cells = TRUE, na.rm = FALSE)

complete_idx <- complete.cases(pred_df[, predictor_vars])
pred_df_complete <- pred_df[complete_idx, ]

X_global <- as.data.frame(pred_df_complete[, predictor_vars, drop = FALSE])

# =========================
# 9. Compute SHAP values
# =========================
# IMPORTANT:
# predcontrib = TRUE must use matrix/data.frame, not SpatRaster directly

shap_mat <- predict(
  xgb_model,
  newdata = as.matrix(X_global),
  predcontrib = TRUE
)

shap_mat <- as.data.frame(shap_mat)

print(colnames(shap_mat))

# Remove baseline / bias column if present
bias_col <- grep("BIAS|Baseline|base|(Intercept)", colnames(shap_mat), ignore.case = TRUE, value = TRUE)

if (length(bias_col) > 0) {
  shap_features <- shap_mat[, setdiff(colnames(shap_mat), bias_col), drop = FALSE]
} else {
  shap_features <- shap_mat
}

# Keep predictor columns in the same order
shap_features <- shap_features[, predictor_vars, drop = FALSE]

# =========================
# 10. Derive dominant variable
# =========================
abs_shap <- abs(shap_features)

dominant_idx <- apply(abs_shap, 1, which.max)
dominant_var <- predictor_vars[dominant_idx]

pred_df_complete$dominant_var <- dominant_var

# =========================
# 11. Derive dominance strength
# =========================
sorted_vals <- t(apply(abs_shap, 1, function(z) sort(z, decreasing = TRUE)))
pred_df_complete$dominance_strength <- sorted_vals[, 1] - sorted_vals[, 2]

# =========================
# 12. Convert results back to rasters
# =========================
var_levels <- data.frame(
  value = seq_along(predictor_vars),
  variable = predictor_vars
)

dominant_code <- match(pred_df_complete$dominant_var, var_levels$variable)

dominant_r <- rast(pred_stack[[1]])
values(dominant_r) <- NA_integer_
values(dominant_r)[pred_df_complete$cell] <- dominant_code
names(dominant_r) <- "Dominant_variable"
levels(dominant_r) <- var_levels

dom_strength_r <- rast(pred_stack[[1]])
values(dom_strength_r) <- NA_real_
values(dom_strength_r)[pred_df_complete$cell] <- pred_df_complete$dominance_strength
names(dom_strength_r) <- "Dominance_strength"

# =========================
# 13. Area statistics of dominant variables
# =========================
# 13.1 By pixel count
area_stat_pixels <- pred_df_complete %>%
  dplyr::count(dominant_var, name = "n_pixels") %>%
  dplyr::mutate(
    proportion = n_pixels / sum(n_pixels),
    percent = proportion * 100
  )

print(area_stat_pixels)

# 13.2 By actual area (km2)
cell_area_r <- cellSize(pred_stack[[1]], unit = "km")
pred_df_complete$cell_area_km2 <- values(cell_area_r)[pred_df_complete$cell]

area_stat_km2 <- pred_df_complete %>%
  dplyr::group_by(dominant_var) %>%
  dplyr::summarise(
    area_km2 = sum(cell_area_km2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    proportion = area_km2 / sum(area_km2),
    percent = proportion * 100
  )

print(area_stat_km2)

write.csv(area_stat_pixels, "Plots/Dominant_variable_area_stat_by_pixels.csv", row.names = FALSE)
write.csv(area_stat_km2, "Plots/Dominant_variable_area_stat_by_km2.csv", row.names = FALSE)

# =========================
# 14. Visualization
# =========================
dom_cols <- c("#4C78A8", "#54A24B", "#E45756", "#8C6D31")
dom_cols <- c("#4C78A8", "#54A24B", "#E45756", "#8C6D31", "#D98C2E", "#7BB874", "#4C8CA6", "#B26DA3")

names(dom_cols) <- predictor_vars

plot(
  dominant_r,
  col = dom_cols,
  type = "classes",
  main = "Dominant variable"
)

plot(
  dom_strength_r,
  col = hcl.colors(100, "YlOrRd", rev = FALSE),
  main = "Dominance strength"
)

barplot(
  height = area_stat_km2$percent,
  names.arg = area_stat_km2$dominant_var,
  col = dom_cols[area_stat_km2$dominant_var],
  las = 2,
  ylab = "Dominated area (%)",
  main = "Area fraction dominated by each variable"
)

# =========================
# 15. Write outputs
# =========================
writeRaster(dominant_r, "Plots/Dominant_variable_SHAP.tif", overwrite = TRUE)
writeRaster(dom_strength_r, "Plots/Dominance_strength_SHAP.tif", overwrite = TRUE)

png("Plots/Dominant_variable_SHAP.png", width = 2400, height = 1400, res = 200)
plot(
  dominant_r,
  col = dom_cols,
  type = "classes",
  main = "Dominant variable"
)
dev.off()

png("Plots/Dominance_strength_SHAP.png", width = 2400, height = 1400, res = 200)
plot(
  dom_strength_r,
  col = hcl.colors(100, "YlOrRd", rev = FALSE),
  main = "Dominance strength"
)
dev.off()

png("Plots/Dominant_variable_area_fraction.png", width = 1800, height = 1400, res = 200)
barplot(
  height = area_stat_km2$percent,
  names.arg = area_stat_km2$dominant_var,
  col = dom_cols[area_stat_km2$dominant_var],
  las = 2,
  ylab = "Dominated area (%)",
  main = "Area fraction dominated by each variable"
)
dev.off()


# =========================
# Save SHAP rasters for each predictor
# =========================
shap_rasters <- list()

for (v in predictor_vars) {
  r_tmp <- rast(pred_stack[[1]])
  values(r_tmp) <- NA_real_
  values(r_tmp)[pred_df_complete$cell] <- shap_features[[v]]
  names(r_tmp) <- paste0("SHAP_", v)
  
  shap_rasters[[v]] <- r_tmp
  
  writeRaster(
    r_tmp,
    filename = paste0("Plots/SHAP_", v, ".tif"),
    overwrite = TRUE
  )
}


# ===layer# =========================
# Save absolute SHAP rasters for each predictor
# =========================
asb_shap_rasters <- list()

for (v in predictor_vars) {
  r_tmp <- rast(pred_stack[[1]])
  values(r_tmp) <- NA_real_
  values(r_tmp)[pred_df_complete$cell] <- abs_shap[[v]]
  names(r_tmp) <- paste0("AbsSHAP_", v)
  abs_shap_rasters[[v]] <- r_tmp
  writeRaster(
    r_tmp,
    filename = paste0("Plots/AbsSHAP_", v, ".tif"),
    overwrite = TRUE
  )
}

