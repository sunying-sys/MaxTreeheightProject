# ============================================================
# Step 02: Prepare species-level maximum tree height records
# ============================================================

library(data.table)
library(dplyr)
library(stringr)
library(rWCVP)
library(rWCVPdata)
library(plantlist)

setwd("D:/MaxTreeHeightProject")

data("wcvp_names", package = "rWCVPdata")
data("orders_dat", package = "plantlist")

# 1. Load and harmonize tree height records from TRY, Tallo and BIEN data.
tree_height_try <- fread("Data/OriginalData/Rawdata/Cleaned_tree_height_data_from_TRY.csv")[, -1, with = FALSE]
tree_height_merged <- fread("Data/OriginalData/Rawdata/Merged_tree_height_data_for_individuals.csv")

load("D:/DataCollection/BIEN_traits/BIEN_PlantHeight.RData")
tree_height_bien <- as.data.table(tr_height_sub[, c("longitude", "latitude", "trait_value", "scrubbed_species_binomial")])
setnames(tree_height_bien, c("x", "y", "Height", "Binomial"))

required_cols <- c("x", "y", "Height", "Binomial")
tree_height <- rbindlist(list(tree_height_try[, ..required_cols], tree_height_merged[, ..required_cols], tree_height_bien), use.names = TRUE)
tree_height <- unique(na.omit(tree_height, cols = required_cols))
tree_height[, `:=`(x = as.numeric(x), y = as.numeric(y), Height = as.numeric(Height), Binomial = as.character(Binomial))]

# 2. Parse species names and match them to WCVP accepted names.
species_df <- data.table(full_name = unique(tree_height$Binomial))
species_df[, scientific_name := str_extract(full_name, "^[A-Z][a-zA-Z-]+\\s+[a-zA-Z-]+")]
species_df[, authority := str_trim(str_remove(full_name, "^[A-Z][a-zA-Z-]+\\s+[a-zA-Z-]+"))]
species_df <- species_df[!is.na(scientific_name)]

# This step can be time-consuming for large species lists.
matched_wcvp <- as.data.table(wcvp_match_names(names_df = as.data.frame(species_df[, .(full_name, scientific_name)]), name_col = "scientific_name"))

# 3. Add family and genus information from WCVP and remove hybrid or unmatched names.
matched_wcvp[, Family := wcvp_names$family[match(wcvp_name, wcvp_names$taxon_name)]]
matched_wcvp[, Genus := wcvp_names$genus[match(wcvp_name, wcvp_names$taxon_name)]]

species_wcvp <- matched_wcvp[!grepl("×", wcvp_name) & !is.na(wcvp_name) & !is.na(Family), .(full_name, wcvp_name, Family, Genus)]
species_wcvp <- unique(species_wcvp)[order(full_name, wcvp_name), .SD[1], by = full_name]
fwrite(species_wcvp, "Data/SpeciesNames_WCVP_Matched.csv")

# 4. Merge WCVP-standardized species names with individual tree height records.
tree_height_matched <- merge(tree_height, species_wcvp, by.x = "Binomial", by.y = "full_name", all.x = FALSE)
tree_height_matched <- unique(na.omit(tree_height_matched, cols = c("x", "y", "Height", "Binomial", "wcvp_name", "Family", "Genus")))
save(tree_height_matched, file = "Data/TreeHeight_IndividualRecords_WCVPmatched_20260410.RData")

# 5. Calculate species-level height statistics and geographic ranges.
load("Data/Tree_height_merged_20260410_for_single_trees_spWCVPmatched.RData")
treeHeightTable_matched[, Height := as.numeric(Height)]
tree_height_matched <- treeHeightTable_matched[Height >= 2 & Height <= 130]

species_height <- tree_height_matched[, .(
  n_records = .N,
  maxHeight = max(Height, na.rm = TRUE),
  q975Height = quantile(Height, 0.975, na.rm = TRUE),
  lon_range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE),
  lat_range = max(y, na.rm = TRUE) - min(y, na.rm = TRUE),
  meanHeight = mean(Height, na.rm = TRUE),
  x = x[which.max(Height)],
  y = y[which.max(Height)],
  Family = Family[which.max(Height)],
  Genus = Genus[which.max(Height)]
), by = wcvp_name]

setnames(species_height, "wcvp_name", "species")

# 6. Retain species with at least three records and either high upper-tail height values or broad geographic ranges.
mean_max_height <- mean(species_height$maxHeight, na.rm = TRUE)
species_height_filtered <- species_height[n_records >= 3 & (maxHeight > mean_max_height | lon_range > 9 | lat_range > 9)]

# 7. Add plant group information, including Angiosperms and Gymnosperms.
orders_table <- as.data.table(orders_dat)[, .(FAMILY, ORDER, GROUP)]
species_height_final <- merge(species_height_filtered, orders_table, by.x = "Family", by.y = "FAMILY", all.x = TRUE)
species_height_final[GROUP == "", GROUP := "Angiosperms"]
species_height_final[is.na(GROUP), GROUP := "Unknown"]
species_height_final <- unique(species_height_final)

# 8. Export the final species-level maximum tree height table.
fwrite(species_height_final, "Data/MaxTreeheight_20260410_for_single_trees_filter_by_meanmaxHeight_lonlatrange_gt3.csv")
