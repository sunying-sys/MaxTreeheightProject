library(raster)
library(terra)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("D:/MaxTreeHeightProject")

# load the predicted mean map 
maxHeightMap = rast("Data/layer_data/maxHeight_mean_Modelling_20260407_merged.tif")
# tranform into data frame and keep x and y
dataframeHeightMap = as.data.frame(maxHeightMap,xy=T,na.rm=T)
names(dataframeHeightMap) = c("MaxHeight","x","y")
fwrite(dataframeHeightMap,"Plots/CSVTransformed/raster_to_table_for_max_height.csv")

dataframeHeightMap = fread("Plots/CSVTransformed/raster_to_table_for_max_height.csv")
# round the latitue
aggregatedTable = dataframeHeightMap %>%
  mutate(y = round(y, 1)) %>%
  group_by(y) %>%
  summarise(
    MaxHeightMean = mean(MaxHeight, na.rm = TRUE),
    MaxCI = max(MaxHeight, na.rm = TRUE),
    highCI = quantile(MaxHeight, probs = 0.975, na.rm = TRUE),
    lowCI = quantile(MaxHeight, probs = 0.025, na.rm = TRUE))

write.csv(aggregatedTable,"Plots/Aggregated_MaxHeight_along_latitude.csv")

MaxHeightMap_potential = rast("Data/layer_data/adjustedPotentialMaxheight_Mean_Modelling_20260407_merged.tif")
# tranform into data frame and keep x and y
dataframeHeightMap_potential = as.data.frame(MaxHeightMap_potential,xy=T,na.rm=T)
names(dataframeHeightMap_potential) = c("MaxHeight","x","y")
fwrite(dataframeHeightMap_potential,"Plots/raster_to_table_for_max_height.csv")

dataframeHeightMap = fread("Plots/raster_to_table_for_max_height_potential.csv")
# round the latitue
aggregatedTable_potential = dataframeHeightMap %>%
  mutate(y = round(y, 1)) %>%
  group_by(y) %>%
  summarise(
    MaxHeight_potential = mean(MaxHeight, na.rm = TRUE),
    MaxCI = max(MaxHeight, na.rm = TRUE),
    highCI = quantile(MaxHeight, probs = 0.975, na.rm = TRUE),
    lowCI = quantile(MaxHeight, probs = 0.025, na.rm = TRUE))
fwrite.csv(aggregatedTable_potential,"Plots/Aggregated_MaxHeight_potential_along_latitude.csv")



latitudeTable_MaxH = fread("ProcessData/Aggregated_MaxHeight_along_latitude.csv")[,-1] %>%mutate_at(vars(c('lowCI')),~ifelse(lowCI <=0, 0, .))
latitudeTable_MaxQ = fread("ProcessData/Aggregated_MaxHeight_potential_along_latitude_potential.csv")[,-1] %>%mutate_at(vars(c('lowCI')),~ifelse(lowCI <=0, 0, .))
latitudeTable_MaxH$MaxCI = latitudeTable_MaxQ$MaxHeight_potential
colnames(latitudeTable_MaxH)
range(latitudeTable_MaxH$lowCI)
range(latitudeTable_MaxH$highCI)
# display.brewer.pal(n, name)
p1 = latitudeTable_MaxH %>% 
  ggplot(aes(x = y, y= MaxHeightMean)) + 
  geom_ribbon(aes(ymin=lowCI, ymax=highCI), fill="darkcyan",alpha=0.3,color=NA)+
  # scale_fill_manual(values = c("GS_Max1" = "#009e73","GS_Max2" = "#009e73","GS_Mean1" = "#cc79a7","GS_Mean2" = "#cc79a7","HM1"="#56b4e9","HM2"= "#56b4e9","SD1"="#e69f00","SD2"="#e69f00","WK1"="#d55e00","WK2"= "#d55e00"))+
  # geom_hline(yintercept=0.6,linetype="dashed", color="gray30", size=0.5)+
  geom_line(aes(y = MaxHeightMean, color = "Line 1"))+
  geom_line(aes(y = MaxCI, color = "Line 2"),alpha =0.5)+
  scale_color_manual(values = c("darkcyan","#FF00FF"))+
  scale_linetype_manual(values=c(rep(c("solid","solid"), 5)))+
  scale_x_continuous(expand = c(0,0),limits = c(-60, 80), breaks = c(-60,-40,-20,0,20,40,60,80))+
  scale_y_continuous(expand = c(0,0),limits = c(0, 75), breaks = c(0,25,50,75))+
  # scale_alpha_discrete(range = c(0.2, 0.7))+
  ylab("") +
  theme_bw()+
  coord_flip() + 
  #xlim(c(24,61))+
  # ylim(0.3,0.9)+
  xlab("Latitude")+
  # guides( col = guide_legend(ncol=2),linetype = guide_legend())+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.title = element_text(color = "black",size = 20),
        axis.text = element_text(color = "black",size = 20),
        legend.position="none")
p1 

ggsave("Plots/Figure_05_02_treeMaxheight_alongLatitude.pdf",plot = p1,width = 6, height = 10, units = "in")


