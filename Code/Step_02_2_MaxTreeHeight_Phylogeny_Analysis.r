# load the package which could provide the phylogeny information for all the species
# here is the link : https://github.com/jinyizju/V.PhyloMaker
library(V.PhyloMaker)
library(data.table)
library(ape)
library(ggtree)
library(tidytree)
library(phytools)
library(geiger) #for matching tree and traits
library(diversitree)
library(RColorBrewer)
library(phylocomr)
library(phylosignal)
library(phylobase)
library(plantlist)
library(stringr)
library(castor)
library(dplyr)
library(purrr)
# set the working directory
setwd("D:/MaxTreeHeightProject") 

###################################################################################################################
# STEP1 Phylogeny analysis and plotting for all the species
###################################################################################################################
# load the table of wood density database with species name uniformed in TNRS
woodHeightDataRaw = fread("Data/MaxTreeheight_20260410_for_single_trees_filter_by_meanmaxHeight_lonlatrange_gt3.csv") %>% 
  dplyr::select(species,maxHeight,Family,Genus,ORDER) %>%
  rename(family = Family,genus =Genus) %>%
  na.omit()
# get the phylogeny information
woodHeightData = woodHeightDataRaw
phyloResult = phylo.maker(woodHeightData,scenarios = "S3",output.tree=T)
# phyloResult$scenario.3
phyloTree = phyloResult$tree.scenario.3
# phyloTree = chronos(phyloTree, lambda = 0, model = "relaxed") 
# write the tree to the local folder
write.tree(phyloTree, "PhylogenyAnalysis/Phylogeny_t ee_of_vascular_plants_20260410.tre")
phyloTree <- read.tree("PhylogenyAnalysis/Phylogeny_tree_of_vascular_plants_20260410.tre")
# match tree species to tip.labble

woodHeightData <- as.data.frame(woodHeightData)
rownames(woodHeightData)=gsub(" ", "_", woodHeightData$species)
dat1 = geiger::treedata(phyloTree,as.matrix(woodHeightData))
cleanedWoodHeightData = as.data.frame(dat1$data)
newTree = dat1$phy

# write.nexus(newTree,file="nexus.trees")
# newTree = read.nexus(file="nexus.trees")
write.tree(newTree, "PhylogenyAnalysis/Phylogeny_tree_of_woody_plants_20260410.tre")
newTree<- read.tree("PhylogenyAnalysis/Phylogeny_tree_of_woody_plants_20260410.tre")

# newTree = chronos(chronoTree, lambda = 0, model = "relaxed") 
newWoodHeightData <- woodHeightData[
  match(newTree$tip.label, gsub(" ", "_", woodHeightData$species)),
] %>%
  mutate(maxHeight = as.numeric(as.character(maxHeight))) %>%
  na.omit()

traitVector = newWoodHeightData$maxHeight
names(traitVector)<-newTree$tip.label

traitTable = data.frame(name= as.vector(newTree$tip.label), maxHeight = newWoodHeightData$maxHeight,maxHeight1 = newWoodHeightData$maxHeight) %>% purrr::modify_if(is.factor, as.character) 
# run the aot analysis
aotResult = ph_aot(traits=traitTable,phylo=newTree,trait_contrasts = 1, randomizations = 999,ebl_unstconst = F)
# get the conservatism result
conservTable = aotResult$trait_conservatism %>% filter(trait.name == "maxHeight") # %>% filter(ntaxa<574&ntaxa>50)
# table(conservTable$ntaxa)
# sum(conservTable$ntaxa > 50)
calculationTable = fread("PhylogenyAnalysis/Order_Level_Lambda_BlombergsK_Table_20260410.csv")[,-1]

orderNames = calculationTable$Order

# oder level tree subset and phylogeny analysis

# tip.mn == orderLevelTable$WoodHeightMax  ???不等于

orderLevelNodeFindingFunc = function(ord)
{
  # subet the data frame
  orderLevelTable = calculationTable %>% filter(Order == ord)
  # subset the conservTable bu tip.mn and ntaxa
  orderLevelConservTable = conservTable %>% filter(tip.mn == orderLevelTable$WoodHeightMax& ntaxa==orderLevelTable$SpeciesNumber) %>% dplyr::select(ntaxa,tip.mn,tmn.ranklow,tmn.rankhi,node.mn)
  # get the node position 
  nodePosition = find_root_of_monophyletic_tips(newTree, newWoodHeightData %>% filter(ORDER ==ord) %>% rownames(), as_MRCA=TRUE, is_rooted=FALSE)
  # 
  if (nrow(orderLevelConservTable) ==0)
  {
    ValCompare = NA
  }else
  {
    if (orderLevelConservTable$tmn.ranklow<25|orderLevelConservTable$tmn.rankhi>975)
    {
      ValCompare = "Lower"
    }else 
      if (orderLevelConservTable$tmn.ranklow>975|orderLevelConservTable$tmn.rankhi<25)
      {
        ValCompare = "Higher"
      }else 
      {
        ValCompare = NA
      }
  }
  
  return(data.frame(cbind(orderLevelTable,orderLevelConservTable),NodePosition = nodePosition,ValueCompare = ValCompare))
}

# use the lapply to sun the calculation
calculationList = lapply(orderNames,orderLevelNodeFindingFunc)
# rbind the result
orderPositionTable = rbindlist(calculationList)
write.csv(orderPositionTable,"PhylogenyAnalysis/Order_Level_Phy_AOT_analysis_result_Table_20260410.csv")
# lets get the information for each 
# get the higher or lower wood density clades
# lowerTable = conservTable %>% filter(nmn.ranklow<25 &nmn.rankhi>975) %>% filter(ntaxa<600&ntaxa>50)
# higherTable = conservTable %>% filter(nmn.ranklow<975 &nmn.rankhi>25) %>% filter(ntaxa<600&ntaxa>50)

# Lambda
phylosig(newTree,traitVector,method="lambda", test=T)

# Phylogenetic signal lambda : 0.768924 
# logL(lambda) : -11397.5 
# LR(lambda=0) : 606.08 
# P-value (based on LR test) : 7.96885e-134 

# define the vector which contains the order names 
habitat = setNames(newWoodHeightData$ORDER,rownames(newWoodHeightData))
n = Ntip(newTree)
colourOrders = calculationTable %>% filter(SpeciesNumber>=5) %>% dplyr::select(Order) %>% unlist()
grayOrders = calculationTable %>% filter(SpeciesNumber<5) %>% dplyr::select(Order)%>% unlist()
# ranked vector of order names 
fullOrderVector = unique(newWoodHeightData$ORDER)
colourOrdersUp = fullOrderVector[fullOrderVector %in% colourOrders]
# allocate the colour by order names
#col.hab = setNames(c(colorRampPalette(brewer.pal(8, "Set1"))(length(colourOrdersUp)),rep("gray15",length(grayOrders))),c(colourOrdersUp,grayOrders))
col.hab = setNames(c(colorRampPalette(brewer.pal(8, "Set1"))(length(colourOrdersUp)),rep("white",length(grayOrders))),c(colourOrdersUp,grayOrders))

# construct the plot object
obj <- contMap(newTree, traitVector, plot = FALSE, fsize = 0.2, lwd = 1)
#plot(obj, fsize = 0.2, lwd = 3)
save(obj,file = "D:/MaxTreeHeightProject/PhylogenyAnalysis/contMap_obj.RData")
load("D:/MaxTreeHeightProject/PhylogenyAnalysis/contMap_obj.RData")

#obj = contMap(newTree, traitVector, plot=F, fsize=0.2, tip.labels=T,lwd=1)#type="fan", 
# define a colour pallete
# rbPal = colorRampPalette(c("black","darkred",brewer.pal(9, "RdYlBu")))(1001) %>% rev()
rbPal = colorRampPalette(c("black",'#67001F', '#B2182B', '#D6604D', '#F4A582','#D1E5F0', '#92C5DE', '#2166AC'))(1001) %>% rev()

#This adds a column of color values
# based on the wood density values
obj$cols = setNames(c(rbPal),0:1000)
plottingOrder = orderPositionTable %>% filter(!is.na(ValueCompare)) 
compareCol = setNames(c(ifelse(plottingOrder$ValueCompare=="Higher", "#cb1b16", "#1368aa")),plottingOrder$ValueCompare)
# ttt= setNames(rbPal(10)[as.numeric(cut(x,breaks = 10))],newTree$tip.label)
pdf("Plots/Figure_SX_Phylogeny_and_traits_plot_latest_20260407_1.pdf",width = 50,height=50)
par(#bg="gray15",
  xpd=TRUE)
plotTree.wBars(obj$tree, traitVector, fsize=0.2, scale=1,width=0.5,lwd=2, tip.labels=F,
               method="plotSimmap", colors=obj$cols,type="fan",outline=F,border="transparent",mar=c(2,2,2,25))#type="fan",

#nodelabels(text=str_pad(round(plottingOrder$tip.mn,2), 4, pad = "0",side = c("right")),node=plottingOrder$NodePosition,frame="circle",col="white",cex=2.8,bg=compareCol) #+Ntip(obj$tree)
objTT<-get("last_plot.phylo",envir=.PlotPhyloEnv)
for(i in 1:n)
{cc = if(objTT$xx[i]>0) 120 else -120
th = atan(objTT$yy[i]/objTT$xx[i])
segments(#objTT$xx[i],objTT$yy[i],objTT$xx[i]+cc*cos(th),objTT$yy[i]+cc*sin(th),
  objTT$xx[i]+cc*cos(th),objTT$yy[i]+cc*sin(th),objTT$xx[i]+1.15*cc*cos(th),objTT$yy[i]+1.15*cc*sin(th),
  lwd=10,lend=2,col=col.hab[habitat[newTree$tip.label[i]]])
}

usr <- par("usr")  
add.color.bar(160,cols=obj$cols,title=expression("Wood height (m)"),lims = obj$lims,fcol= "white",digits=3,lwd=60,fsize=6,prompt=FALSE,x=340,y=400,subtitle="",outline=T)
legend(x=450,y=370,colourOrdersUp,pch=15,col=c(col.hab[1:length(colourOrdersUp)]),pt.cex=6,cex=4,bty="o",ncol = 1,text.col="black")

dev.off()

###################################################################################################################
# STEP 2 Phylogeny analysis and plotting for each order
###################################################################################################################
# use the order names as the identifier to subset the data frame and do anaylsis
orderNames = unique(woodHeightData$ORDER)
# oder level tree subset and phylogeny analysis
table(woodHeightData$ORDER)
ord<-orderNames[1]
orderLevelFunc = function(ord)
{
  # subet the data frame
  orderLevelTable = woodHeightData %>% filter(ORDER == ord)
  if(nrow(orderLevelTable)>=5)
  {
    # add the row names 
    rownames(orderLevelTable)=gsub(" ", "_", orderLevelTable$species)
    # use the data frame to get the phylogeny information by the treedata function
    orderLevelTreeData = geiger::treedata(phyloTree,orderLevelTable)
    # get the tree data paired with traits information
    orderWoodHeightData = as.data.frame(orderLevelTreeData$data)
    # get the paired phylogeny tree
    orderTree = orderLevelTreeData$phy
    # write the tree and data frame to the local folder
    write.csv(orderWoodHeightData,paste("PhylogenyAnalysis/OrderLevelWDTables/Order_",ord,"_Subset_wood_Height_data_20260410.csv",sep=""))
    # and write the tree
    write.tree(orderTree, paste("PhylogenyAnalysis/OrderLevelTrees/Order_",ord,"_Subset_wood_Height_data_20260410.tre",sep=""))
    # do the test for phylogeny
    traitVector = as.numeric(as.character(orderWoodHeightData$maxHeight))
    names(traitVector) =  orderTree$tip.label
    # apply the test
    lambdaResult = phylosig(orderTree,traitVector,method="lambda", test=T)
    blombergsResult = phylosig(orderTree,traitVector,method="K", test=T)
    maxWH = mean(as.numeric(orderLevelTable$maxHeight))
    #round(mean(orderLevelTable$WoodDensity),3)
    # define the output table 
    orderLevelOutput = data.frame(Order = ord,
                                  SpeciesNumber = nrow(orderLevelTable),
                                  Lambda = lambdaResult$lambda,
                                  LambdaP_Val = lambdaResult$P,
                                  BlombergsK = blombergsResult$K,
                                  BlombergsP_Val = blombergsResult$P,
                                  WoodHeightMax = maxWH)
    
    
  }else
  {
    orderLevelOutput = data.frame(Order = ord,
                                  SpeciesNumber = nrow(orderLevelTable),
                                  Lambda = NA,
                                  LambdaP_Val = NA,
                                  BlombergsK = NA,
                                  BlombergsP_Val = NA,
                                  WoodHeightMax = mean(as.numeric(orderLevelTable$maxHeight)))
    
  }
  return(orderLevelOutput)
  
}
# use the lapply to sun the calculation
calculationList = lapply(orderNames,orderLevelFunc)
# rbind the result
calculationTable = rbindlist(calculationList)
# write to local folder
write.csv(calculationTable,"PhylogenyAnalysis/Order_Level_Lambda_BlombergsK_Table_20260410.csv")

