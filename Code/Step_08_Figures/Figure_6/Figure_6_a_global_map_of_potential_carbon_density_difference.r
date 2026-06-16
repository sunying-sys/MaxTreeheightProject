library(raster)
library(sp)
library(dplyr)
library(data.table)
library(ggplot2)

setwd("/Volumes/EcoComp01/sunying/Carbon_biomass_compare/Compare_Carbon_potential")

# load the map-----
potentialCarbon_ESA <- raster("/Volumes/EcoComp01/sunying/Carbon_biomass_compare/RasterFiles/Global_Carbon_Potential_ESA_spatSample_smallvalue_to_originalvalue_Cfactor0.5_linear_BIENdata20260415.tif")
potentialCarbon_WD <- raster("/Volumes/EcoComp01/sunying/Carbon_biomass_compare/RasterFiles/Global_Carbon_Potential_WD_spatSample_smallvalue_to_originalvalue_Cfactor0.5_linear_BIENdata20260415.tif")
ESA_C <- raster("/Volumes/EcoComp01/sunying/Carbon_biomass_compare/RasterFiles/Global_Carbon_ESA_Cfactor0.5.tif")
WD_C <- raster("/Volumes/EcoComp01/sunying/Carbon_biomass_compare/RasterFiles/Global_Carbon_WD_Cfactor0.5.tif")
diffESA <- potentialCarbon_ESA - ESA_C
diffWD <- potentialCarbon_WD - WD_C
plot(diffESA)
plot(diffWD)
diff_mean <- mean(stack(diffESA, diffWD), na.rm = TRUE)
qs <- quantile(diff_mean,probs = c(0.025, 0.05,0.25,0.75,0.95,0.975),na.rm = TRUE)
# 2.5%       5%      25%      75%      95%    97.5% 
#   0.00000  0.00000  0.00000 29.98173 41.69479 45.17385 
# X2.5.   X97.5.
# mean     0 32.15016
diff_mean1 <- clamp(diff_mean, upper = 40)
diff_mean1[diff_mean1 == 0] <- NA
plot(diff_mean1)


diff_mean1 <-raster("/Volumes/EcoComp01/sunying/Carbon_biomass_compare/RasterFiles/Global_Carbon_meandiff_treeheightPropotion_Cfactor0.5_valuerange_1_5.tif")
# define the equal earth projection
equalEarthProj = "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# project raster to equal earth
extentInfo = raster::projectExtent(diff_mean1, crs = equalEarthProj)
# reproject the rasters
diff_meanprojected = raster::projectRaster(diff_mean1, to=extentInfo,over=T)

worldBorderInit = shapefile("/Volumes/EcoComp01/sunying/Fig5_Global_map_maximum_tree_height/world_continent_shp/continent.shp")
worldBorderRaw = crop(worldBorderInit,c(-180, 180, -60, 84))
# reproject the world border
worldEqualEarth = spTransform(worldBorderRaw,to=extentInfo, CRS(equalEarthProj),over=T)

diff_meanprojectedFinal = mask(diff_meanprojected,worldEqualEarth)
# make the plot and write into local PDF
pdf("Figure_06_1_compare_carbon_potential_BIENdata_linearModel.pdf",width =12, height=8)
#par(mfrow=c(3,1))
rbPal1= colorRampPalette(c("#FDE725FF","#3CBB75FF","#2D708EFF","#481567FF"))
plot(worldEqualEarth,border="gray70",col="gray70",axes=T,main="Carbon potential",lwd=0.2)
plot(diff_meanprojectedFinal,col=rbPal1(100),add=T,breaks=c(seq(0,40,0.40)),legend=T,maxpixels = 50000000)
dev.off() 



