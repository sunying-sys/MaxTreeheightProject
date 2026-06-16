library(terra)
library(sp)
library(dplyr)
library(data.table)
library(ggplot2)

setwd("D:/MaxTreeHeightProject")

data <- fread("ProcessData/Summary_potentialCarbon_radiopropotion.csv")
colnames(data)
cor(data$mean_Science,data$mean_WD)
#[1] 0.8509745
cor(data$mean_Science,data$mean_ESA)
# [1] 0.8174079
my_cols <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
  "#66a61e", "#e6ab02", "#a6761d", "#666666",
  "#1f78b4", "#b15928", "#33a02c", "#fb9a99",
  "#cab2d6", "#6a3d9a")

data$subforestType <- factor(data$subforestType,levels = c("Tropical moist","Tropical dry","Tropical coniferous","Temperate broadleaf","Temperate coniferous","Boreal","Tropical savannas","Temperate savannas","Flooded savannas","Montane grassland","Tundra","Mediterranean forests","Desert","Mangroves"))
plotCodeVdataplotCodeVer01 = ggplot(aes(x=mean_WD,y=mean_Science,color=factor(subforestType)),data=data)+
  geom_point(size =3)+
  # geom_errorbarh(aes(xmin = mean_WD-sd_WD,xmax = mean_WD+sd_WD),height = 0.01)+
  # geom_errorbar(aes(ymin = mean_Science-sd_Science,ymax = mean_Science+sd_Science),width = 0.01) +
  geom_errorbarh(aes(xmin = q25_WD,xmax = q75_WD),height = 0.01)+
  geom_errorbar(aes(ymin = q25_Science,ymax = q75_Science),width = 0.01) +
  geom_abline(slope=1, intercept=0,show.legend = NA, linetype="dashed",linewidth = 0.5)+
  coord_cartesian(xlim =c(0,125), ylim = c(0,125))+
  #scale_x_continuous(limits = c(0, 130), breaks = c(0, 25, 50, 75,100,125), expand = expansion(mult = c(0, 0.03))) +
  #scale_y_continuous(limits = c(0, 130), breaks = c(0, 25, 50, 75,100,125), expand = expansion(mult = c(0, 0.03)))+
  theme_classic()+
  scale_color_manual(
    values = my_cols,
    name = NULL,
    labels = c(
      "Tropical moist","Tropical dry","Tropical coniferous",
      "Temperate broadleaf","Temperate coniferous","Boreal",
      "Tropical savanna","Temperate savanna","Flooded savanna",
      "Montane grassland","Tundra","Mediterranean forest",
      "Desert","Mangroves"
    )
  ) +
  #scale_colour_discrete(name=NULL,labels= c("Tropical moist","Tropical dry","Tropical coniferous","Temperate broadleaf","Temperate coniferous","Boreal","Tropical savanna","Temperate savanna","Flooded savanna","Montane grassland","Tundra","Mediterranean forest","Desert","Mangroves"))+
  labs(x = "Carbon density from Mo et al. (t C/ha)", y = "Carbon density from Roebroek et al. (t C/ha)") +
  guides(color = guide_legend(title = NULL,nrow=7),alpha = "none")+
  theme(legend.position = "none",
        #legend.position = c(0.7, 0.15),
        legend.text= element_text(size=16),
        legend.title= element_text(size=16),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        legend.key.size = unit(0.8, 'cm'),
        axis.title = element_text(size = 20),
        axis.ticks.length = unit(0.3, "cm") ,
        axis.text = element_text(size = 20))
plotCodeVdataplotCodeVer01
ggsave("Plots/Figure_06_02_Carbon_WD_vs_Science_radio_meadCI25_75.pdf",
       plot = plotCodeVdataplotCodeVer01,
       #device = cairo_pdf,
       width = 10,
       height = 10,
       units = "in")

plotCodeVdataplotCodeVer02 = ggplot(aes(x=mean_ESA,y=mean_Science,color=factor(subforestType)),data=data)+
  geom_point(size =3)+
  # geom_errorbarh(aes(xmin = mean_ESA-sd_ESA,xmax = mean_ESA+sd_ESA),height = 0.01)+
  # geom_errorbar(aes(ymin = mean_Science-sd_Science,ymax = mean_Science+sd_Science),width = 0.01) + 
  geom_errorbarh(aes(xmin = q25_ESA,xmax = q75_ESA),height = 0.01)+
  geom_errorbar(aes(ymin = q25_Science,ymax = q75_Science),width = 0.01) +
  geom_abline(slope=1, intercept=0,show.legend = NA, linetype="dashed",linewidth = 0.5)+
  coord_cartesian(xlim =c(0,125), ylim = c(0,125))+
  #scale_x_continuous(limits = c(0, 145), breaks = c(0, 25, 50, 75,100,125), expand = expansion(mult = c(0, 0.03))) +
  #scale_y_continuous(limits = c(0, 145), breaks = c(0, 25, 50, 75,100,125), expand = expansion(mult = c(0, 0.03)))+
  theme_classic()+
  scale_color_manual(
    values = my_cols,
    name = NULL,
    labels = c(
      "Tropical moist","Tropical dry","Tropical coniferous",
      "Temperate broadleaf","Temperate coniferous","Boreal",
      "Tropical savanna","Temperate savanna","Flooded savanna",
      "Montane grassland","Tundra","Mediterranean forest",
      "Desert","Mangroves")) +
  #scale_colour_discrete(name=NULL,labels= c("Tropical moist","Tropical dry","Tropical coniferous","Temperate broadleaf","Temperate coniferous","Boreal","Tropical savanna","Temperate savanna","Flooded savanna","Montane grassland","Tundra","Mediterranean forest","Desert","Mangroves"))+
  labs(x = "Carbon density from ESA CCI (t C/ha)", y = "Carbon density from Roebroek et al. (t C/ha)") +
  guides(color = guide_legend(title = NULL,nrow=14),alpha = "none")+
  theme(#legend.position = c(0.7, 0.15),
    legend.position = c(0.82, 0.26),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.text= element_text(size=18),
    legend.title= element_text(size=16),
    panel.border = element_rect(color = "black",fill = NA,size = 1),
    legend.key.size = unit(0.8, 'cm'),
    axis.title = element_text(size = 20),
    axis.ticks.length = unit(0.3, "cm") ,
    axis.text = element_text(size = 20))
plotCodeVdataplotCodeVer02
ggsave("Plots/Figure_06_02_Carbon_ESA_vs_Science_radio_meadCI25_75.pdf",plot = plotCodeVdataplotCodeVer02,device = cairo_pdf,width = 10,height = 10,units = "in")

library(patchwork)
p_all <- plotCodeVdataplotCodeVer01 + plotCodeVdataplotCodeVer02 + plot_layout(ncol = 2) 
ggsave(filename = "Plots/Figure_06_02_carbon_density_comparison.pdf",plot = p_all,width = 20,height = 10,units = "in")