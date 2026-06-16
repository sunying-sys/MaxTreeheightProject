# this is the code for Figure 5 making

library(raster)
library(RColorBrewer)
library(viridis)
library(sp)

# maps coould be downloade from zenodo composite
setwd("D:/MaxTreeHeightProject")

# Maxheight map----------
# load the maximum tree height map
maxHeightMerged = raster("D:/MaxTreeHeightProject/Data/layer_data/maxHeight_mean_Modelling_20260407_merged.tif")
# Calculate quantiles (default: 0%, 25%, 50%, 75%, 100%)
# quantile(maxHeightMerged, probs = c(0.025, 0.975), na.rm = TRUE)
# 2.5%    97.5% 
#   26.94963 50.84094 
# modify the values in the raster for further plot making 
maxHeightMerged [maxHeightMerged <25] = 25
maxHeightMerged [maxHeightMerged >50] = 50
# define the equal earth projection
equalEarthProj = "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# project raster to equal earth
extentInfo = raster::projectExtent(maxHeightMerged, crs = equalEarthProj)
# reproject the rasters
maxHeightprojected = raster::projectRaster(maxHeightMerged, to=extentInfo,over=T)

worldBorderInit = shapefile("Data/world_continent_shp/continent.shp")
worldBorderRaw = crop(worldBorderInit,c(-180, 180, -60, 84))
# reproject the world border
worldEqualEarth = spTransform(worldBorderRaw,to=extentInfo, CRS(equalEarthProj),over=T)

maxHeightprojectedFinal = mask(maxHeightprojected,worldEqualEarth)
# make the plot and write into local PDF
pdf("Plots/Figure_05_01_Maximum_tree_height_ensemble_mean_3.pdf",width =12, height=18)
par(mfrow=c(3,1))
rbPal1= colorRampPalette(c("#FDE725FF","#3CBB75FF","#2D708EFF","#481567FF"))
plot(worldEqualEarth,border="gray70",col="gray70",axes=T,main="Maximum tree height",lwd=0.2)
plot(maxHeightprojectedFinal,col=rbPal1(100),add=T,breaks=c(seq(25,50,0.25)),legend=T,maxpixels = 50000000)
dev.off() 

# rbPal2= colorRampPalette(c("#F7F7F7","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))

# plot(worldEqualEarth,border="gray70",col="gray25",axes=T,main="Uncertainty",lwd=0.2)
# plot(modelUncertainWD_Projected,col=rbPal2(100),add=T,breaks=c(seq(0,0.05,0.0005)),legend=T,maxpixels = 40000000)


