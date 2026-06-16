# ============================================================
# Step 01: Prepare pixel-level maximum tree height records
# ============================================================

library(terra)
library(data.table)
library(raster)
library(dplyr)
library(sf)
setwd("D:/MaxTreeHeightProject")

# 1. Create a 30-arcsec global grid template and derive grid IDs from cell numbers.
StartRaster = raster(nrows=17400, ncols=43200, xmn=-180.0001, xmx=179.9999, ymn=-90.00014,ymx=83.99986,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",resolution=c(0.008333333, 0.008333333))
StartRaster[] <- 1:ncell(StartRaster)
grid <- StartRaster

# 2. Load tree height records from TRY, Tallo, and BIEN  data.
treeHeightTable1 = fread("Data/OriginalData/Rawdata/Cleaned_tree_height_data_from_TRY.csv")[,-c(1,5)]
treeHeightTable2 = fread("Data/OriginalData/Rawdata/Merged_tree_height_data_for_individuals.csv")[,-4]
load("D:/DataCollection/BIEN_traits/BIEN_PlantHeight.RData")
treeHeightTable3 = tr_height_sub[,c("longitude","latitude","trait_value")]
colnames(treeHeightTable3)<-c("x", "y", "Height")

# 3. Clean invalid coordinates and biologically implausible height records.
treeHeightTable = rbind(treeHeightTable1,treeHeightTable2,treeHeightTable3)
colnames(treeHeightTable)
treeHeightTable_clean <- treeHeightTable %>%
  mutate(Height = as.numeric(Height)) %>%
  filter(!is.na(x), !is.na(y), !is.na(Height), Height > 0,Height <= 130)
summary(treeHeightTable_clean$Height)

# 4. Assign each tree height record to a 30-arcsec grid cell.
df = treeHeightTable_clean
coordinates(df) <- ~x+y
proj4string(df) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
gridID <- extract(StartRaster, df)
df@data$plotID <- gridID
str(df)

# 5. Extract WWF biome for each tree height record.
wwf <- st_read("Data/Terrestrial_Ecoregions_of_the_World/wwf_terr_ecos.shp")
# convert sf to sp
wwf_sp <- as(wwf[, "BIOME"], "Spatial")
# if CRS is different, transform wwf_sp to match df
if (!identical(sp::proj4string(wwf_sp), sp::proj4string(df))) {
  wwf_sp <- spTransform(wwf_sp, CRSobj = CRS(sp::proj4string(df)))
}
# extract BIOME (point and layer overlap )
df@data$BIOME <- over(df, wwf_sp)$BIOME

# creat data frame
DF_treeHeightTable <- cbind(data.frame(x = coordinates(df)[, 1],y = coordinates(df)[, 2]),df@data)
# transform data.frame to data.table, get speed and memory advantages
setDT(DF_treeHeightTable)
table(DF_treeHeightTable$BIOME)
# 1       2       3       4       5       6       7       8       9      10      11      12      13      14      98      99 
# 234307   43225    9006 4767001 3523874   70115   36594  233676    6659   14705   22757  105418  241861    1613    4385     155 
colnames(DF_treeHeightTable) <- c("x","y","MaxHeight","plotID","WWF_Biome")

save(DF_treeHeightTable, file = "Data/treeHeightTable_WWFextract.RData")
# 
load("Data/treeHeightTable_WWFextract.RData")

# 6. filter MaxHeight based on 97.5% value of same biome and MaxHeight in a plotID 
# (1) keep treeHeight value more than 97.5% value of same biome
biome_q975 <- DF_treeHeightTable[, .(
  q975_treeHeight = quantile(MaxHeight, 0.975, na.rm = TRUE)), by = WWF_Biome]
DF_treeHeightTable_q975 <- merge(DF_treeHeightTable, biome_q975,by = "WWF_Biome",all.x = TRUE)
DF_treeHeightTable_filter <- DF_treeHeightTable_q975[MaxHeight > q975_treeHeight]
# (2) Keep only the record with the maximum Height within each gridID
DF_treeHeightTable_final <- DF_treeHeightTable_filter[order(plotID, -MaxHeight)][, .SD[1], by = plotID]

# 7. output the final data
fwrite(DF_treeHeightTable_final, "D:/MaxTreeHeightProject/Data/Pixel_Level_tree_height_max_20260416_filtered_BIEN_biome975.csv")

# 8. Export the pixel-level maximum tree height records as a shapefile.
library(sf)
plot_stat  = fread("D:/MaxTreeHeightProject/Data/Pixel_Level_tree_height_max_20260416_filtered_BIEN_biome975.csv") %>% 
  #na.omit()%>%
  select(x, y, MaxHeight)
plot_stat_sf <- st_as_sf(plot_stat,coords = c("x", "y"),crs = 4326,remove = FALSE)
st_write(plot_stat_sf,"Data/SHPfile/Pixel_Level_tree_height_max_20260407_filtered.shp",delete_layer = TRUE)

