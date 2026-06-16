library(raster)
library(RColorBrewer)
library(data.table)
library(ggpubr)
library(maps)
library(gridExtra)
library(dplyr)
library(cowplot)
library(ggrastr)
library(ggrepel)

setwd("D:/MaxTreeHeightProject")

trainDataTable = fread("Data/Tree_height_covariates_extraction_merged_20260415_filtered_BIENdata_biome975.csv")[,-1]

# make a copy of the table 
rbPal1= colorRampPalette(c("#aafb14ff","#13a655ff","#2D708EFF","#481168ff"))
# load the world border polygon
worldBorderRaw = shapefile("Data/WORLD_BORDERS/TM_WORLD_BORDERS-0.3.shp")
worldBorder = crop(worldBorderRaw,c(-180, 180, -60, 84)) 
# define the equal earth projection
equalEarthProj =  "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# project raster
extentInfo = projectExtent(worldBorder,equalEarthProj)
# reproject the world border
worldEqualEarth = spTransform(worldBorder ,to=extentInfo, CRS(equalEarthProj))
# transform the data into spatial points
duplicateTable = trainDataTable
coordinates(duplicateTable) = ~x+y
proj4string(duplicateTable) = CRS("+proj=longlat")
transformedPoints = spTransform(duplicateTable, to=extentInfo,CRS(equalEarthProj)) 
newTable = data.frame(transformedPoints@coords,MaxH=trainDataTable$MaxHeight) %>% mutate(MaxH=ifelse(MaxH<20, 20,MaxH)) %>%  mutate(MaxH=ifelse(MaxH>80, 80,MaxH)) %>% mutate(x=coords.x1,y=coords.x2,coords.x2=NULL,coords.x1=NULL)
# define the plot code
panelA = ggplot() + geom_polygon(data = worldEqualEarth,aes(x = long, y = lat, group = group),fill = "grey50",color = "grey50",linewidth = 0.1) + 
                    coord_fixed(1)+
                    ggrastr::rasterise(geom_point(data = newTable,aes(x = x, y = y,colour = MaxH),size=2),dpi=1000,dev = "ragg") +
                    scale_colour_gradientn(colours = rbPal1(100),breaks= c(20,50,80),labels=c(expression(""<=20),"50",expression("">=80)))+
                    theme_classic()+
                    labs(tag = "a")+
                    guides(fill=guide_colorbar(frame.colour = "black", ticks.colour = "white"))+
                    theme(legend.position = c(0.14, 0.3),
                        legend.title = element_text(size=20),
                        legend.text = element_text(size=18),
                        legend.key.size = unit(1.6,"line"),
                        plot.margin= unit(c(-1,0.2,-1,2.1), "cm"), 
                        # panel.border = element_rect(colour = "black", fill=NA, size=2),
                        axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
                        axis.ticks.x = element_blank(),axis.text.x = element_blank(),
                        axis.title.x = element_blank(),axis.title.y = element_blank(),
                        axis.line = element_blank(),
                        plot.tag = element_text(size=32),
                        plot.tag.position = c(0.03, 0.95))+
                    labs(colour = expression("Maximum tree height"))  
                    

#########################################################################################################################
################PANEL B and D
#########################################################################################################################
trainDataTable = fread("Data/Tree_height_covariates_extraction_merged_20260415_filtered_BIENdata_biome975.csv")[,-1] %>% filter(WWF_Biome <15)
# make a copy of the table 
aggregatedTableNew = trainDataTable
# generate a new column to save the forest types not the biome
aggregatedTableNew$ForestType = aggregatedTableNew$WWF_Biome
# change the forest type representative number
aggregatedTableNew$ForestType[aggregatedTableNew$ForestType %in% c(1,2,3,7,9,14)] <- 30 #Tropical
aggregatedTableNew$ForestType[aggregatedTableNew$ForestType %in% c(4,5,8,10)] <- 40 #Temperate
aggregatedTableNew$ForestType[aggregatedTableNew$ForestType %in% c(6,11)] <- 50 #Boreal
aggregatedTableNew$ForestType[aggregatedTableNew$ForestType %in% c(12,13)] <- 60 # Dryland

# transfer th
aggregatedTableNew$WWF_Biome = as.factor(aggregatedTableNew$WWF_Biome)
aggregatedTableNew$ForestType = as.factor(aggregatedTableNew$ForestType)

# biome level statistics of wood density 
biomeStat = aggregate(x=aggregatedTableNew[,c("MaxHeight")],by=aggregatedTableNew[,c("WWF_Biome")],FUN=mean)  %>% `colnames<-`(c("Biome","BiomeMaxHeight"))

biomeSD = aggregate(x=aggregatedTableNew[,c("MaxHeight")],by=aggregatedTableNew[,c("WWF_Biome")],FUN=sd)  %>% `colnames<-`(c("Biome","BiomeSD"))
biomeStat$BiomeMaxHeightSD = biomeSD$BiomeSD
# Biome BiomeMaxHeight BiomeMaxHeightSD
# 1      1       49.71033        10.501706
# 2      2       31.17227        10.180254
# 3      3       28.31806         3.468160
# 4      4       33.25789         4.891902
# 5      5       42.06626         9.024566
# 6      6       24.78185         3.250887
# 7      7       34.25331        10.722649
# 8      8       29.41190         4.345652
# 9      9       23.24980         1.440628
# 10    10       47.20294        10.191304
# 11    11       40.23360         6.575110
# 12    12       35.25984        12.684386
# 13    13       24.56944         5.235493
# 14    14       65.52000        23.379735

# forest type level wood density 
forestTypeStat = aggregate(x=aggregatedTableNew[,c("MaxHeight")],by=aggregatedTableNew[,c("ForestType")],FUN=mean) %>% `colnames<-`(c("ForestType","TypeHeightMean")) %>% mutate(ForestType = c("Tropical","Temperate","Boreal","Dryland"))

forestTypeSD = aggregate(x=aggregatedTableNew[,c("MaxHeight")],by=aggregatedTableNew[,c("ForestType")],FUN=sd) %>% `colnames<-`(c("ForestType","HeightSD"))
forestTypeStat$TypeHeightSD = forestTypeSD$HeightSD

# ForestType TypeHeightMean TypeHeightSD
# 1   Tropical       41.36669    14.161773
# 2  Temperate       36.32123     8.036548
# 3     Boreal       31.35040     9.106118
# 4    Dryland       27.85029     9.623906

write.csv(biomeStat,"BiomeStatistics/Biome_level_MaxHeight_statistics_q975.csv")
write.csv(forestTypeStat,"BiomeStatistics/Forest_type_level_MaxHeight_statistics_q975.csv")

rbPal= colorRampPalette(brewer.pal(9,"RdYlBu"))
aggregatedTableNew <- aggregatedTableNew %>%
  filter(WWF_Biome != c(3,14))
boxPlotMean = ggplot(data = aggregatedTableNew[aggregatedTableNew$WWF_Biome %in% c(1:2,4:13),], aes(x = WWF_Biome, y = MaxHeight,fill=WWF_Biome)) +
            #    geom_jitter(aes(colour = WWF_Biome),alpha = 1/40) +
                 geom_violin(scale = "width", trim = TRUE, alpha = 0.8,drop = FALSE)+ 
                      stat_summary(fun = median, geom = "point", size = 2, color = "black") +
            #    scale_fill_brewer(palette="Spectral")+
            scale_fill_manual(values=rbPal(14))+
                # coord_flip()+
                theme_classic()+
                ylim(5,130)+
                xlab("Biome")+
                labs(tag = "c")+
                scale_x_discrete(labels=c("1" = "Tropical moist",
                                          "2" = "Tropical dry",
                                          "3" = "Tropical coniferous",
                                          "4" = "Temperate broadleaf",
                                          "5" = "Temperate coniferous",
                                          "6" = "Boreal",
                                          "7" = "Tropical savanna",
                                          "8" = "Temperate savanna",
                                          "9" = "Flooded savanna",
                                          "10" = "Montane grassland",
                                          "11" = "Tundra",
                                          "12" = "Mediterranean forest",
                                          "13" = "Desert"
                                          #"14" = "Mangroves"
                                          ),guide = guide_axis(angle = 45))+
                theme(axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.ticks.length=unit(0.15, "cm"),
                        axis.text.x = element_text(colour = "black", size = 18),
                        axis.text.y = element_text(colour = "black", size = 18),
                        plot.tag = element_text(size=32),
                        plot.tag.position = c(0.11, 0.95))+
                theme(legend.position="none")

boxPlotForestTypeMean = ggplot(data = aggregatedTableNew[aggregatedTableNew$ForestType %in% c(30,40,50,60),], aes(x = ForestType, y = MaxHeight,fill=ForestType)) +
                     #    geom_jitter(aes(colour = WWF_Biome),alpha = 1/40) +
                      geom_violin(scale = "width", trim = TRUE, alpha = 0.8)+ 
                      stat_summary(fun = median, geom = "point", size = 2, color = "black") +
                      #    scale_fill_brewer(palette="Spectral")+
            scale_fill_manual(values=rbPal(4))+
                # coord_flip()+
                theme_classic()+
                ylim(5,130)+
            #    scale_color_brewer(palette = "Spectral") +
                ylab("Maximum tree height (m)")+
                # xlab("Forest type")+
                labs(tag = "b")+
                scale_x_discrete(labels=c("30" = "Tropical forest", "40" = "Temperate forest","50" = "Boreal forest","60" = "Dryland forest"),guide = guide_axis(angle = 45))+
                theme(axis.title.x = element_blank(),
                        axis.title.y = element_text(size=20),
                        axis.ticks.length = unit(0.15, "cm"),
                        axis.text.x = element_text(colour = "black", size = 18),
                        axis.text.y = element_text(colour = "black", size = 18),
                        plot.tag = element_text(size=32),
                        plot.tag.position = c(0.3, 0.95))+
                theme(legend.position="none")

# make and save the plot into the local as PDF
pdf("Plots/Fig_1_MaxHeight_Plots_Box+Plot_Panel_a_b_c_q975.pdf",width = 14,height=12)
top_row = plot_grid(panelA, align = "hv")
bottom_row = plot_grid(boxPlotForestTypeMean,boxPlotMean, align = "hv", rel_widths = c(5,14.5))
plot_grid(top_row,bottom_row,ncol=1,rel_heights = c(9,6),align = "hv")
dev.off()
