library(raster)
library(RColorBrewer)
library(viridis)
library(sp)
library(scico)

# maps coould be downloade from zenodo composite
setwd("D:/MaxTreeHeightProject")
# potential map---------
maxHeightMerged_var = raster("D:/MaxTreeHeightProject/Data/layer_data/maxHeight_var_Modelling_20260407_merged.tif")
# Calculate quantiles (default: 0%, 25%, 50%, 75%, 100%)
quantile(maxHeightMerged_var, probs = c(0.025, 0.975), na.rm = TRUE)
# 2.5%    97.5% 
#   0.01156827 0.06363697
# modify the values in the raster for further plot making 
maxHeightMerged_var [maxHeightMerged_var < 0] = 0
maxHeightMerged_var [maxHeightMerged_var >0.05] = 0.05

# define the equal earth projection
equalEarthProj = "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# project raster to equal earth
extentInfo = raster::projectExtent(maxHeightMerged_var, crs = equalEarthProj)
# reproject the rasters
maxHeightprojected = raster::projectRaster(maxHeightMerged_var, to=extentInfo,over=T)

worldBorderInit = shapefile("Data/world_continent_shp/continent.shp")
worldBorderRaw = crop(worldBorderInit,c(-180, 180, -60, 84))
# reproject the world border
worldEqualEarth = spTransform(worldBorderRaw,to=extentInfo, CRS(equalEarthProj),over=T)

maxHeightprojectedFinal = mask(maxHeightprojected,worldEqualEarth)
# make the plot and write into local PDF
pdf("Plots/Figure_X_variation_Maximum_tree_height_ensemble_mean.pdf",width =12, height=18)
par(mfrow=c(3,1))
rbPal1= colorRampPalette(c("#FDE725FF","#3CBB75FF","#2D708EFF","#481567FF"))
rbPal10 = scico(100, palette="vik")     
plot(worldEqualEarth,border="gray70",col="gray70",axes=T,main="Cofficients of variation",lwd=0.2)
plot(maxHeightprojectedFinal,col=rbPal10,add=T,breaks=c(seq(0,0.05,0.0005)),legend=T,maxpixels = 50000000)
dev.off() 