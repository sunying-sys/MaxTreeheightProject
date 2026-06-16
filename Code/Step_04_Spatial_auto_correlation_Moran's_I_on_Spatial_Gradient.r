library(sp)
library(spatstat)
library(ncf)
library(spdep)
library(data.table)
library(parallel)
library(mgcv)
# set the working directoy
setwd("D:/MaxTreeHeightProject")

######################################################################################
# read the train data which doesn't have outliers
trainDataTable  = fread("Data/Tree_height_covariates_extraction_merged_20260416_filtered_biome975.csv")[,-1] %>% na.omit()

# set.seed(1000)
# randomSampledTable = trainDataTable[sample(nrow(trainDataTable ), 50000),]
# # make a copy of the random subsampled table
# duplicateTable = randomSampledTable

duplicateTable = trainDataTable
#  define the train vector
retainedVars = c("Aridity_Index",
                 "CHELSA_Annual_Mean_Temperature",
                 "CHELSA_Annual_Precipitation",
                 "CHELSA_Isothermality",
                 "CHELSA_Max_Temperature_of_Warmest_Month",
                 "CHELSA_Mean_Diurnal_Range",
                 "CHELSA_Mean_Temperature_of_Coldest_Quarter",
                 "CHELSA_Mean_Temperature_of_Driest_Quarter",
                 "CHELSA_Mean_Temperature_of_Warmest_Quarter",
                 "CHELSA_Mean_Temperature_of_Wettest_Quarter",
                 "CHELSA_Min_Temperature_of_Coldest_Month",
                 "CHELSA_Precipitation_Seasonality",
                 "CHELSA_Precipitation_of_Coldest_Quarter",
                 "CHELSA_Precipitation_of_Driest_Month",
                 "CHELSA_Precipitation_of_Driest_Quarter",
                 "CHELSA_Precipitation_of_Warmest_Quarter",
                 "CHELSA_Precipitation_of_Wettest_Month",
                 "CHELSA_Precipitation_of_Wettest_Quarter",
                 "CHELSA_Temperature_Annual_Range",
                 "CHELSA_Temperature_Seasonality",
                 "Depth_to_Water_Table",
                 "EarthEnvCloudCover_MODCF_interannualSD",
                 "EarthEnvCloudCover_MODCF_intraannualSD",
                 "EarthEnvCloudCover_MODCF_meanannual",
                 "EarthEnvTopoMed_Eastness",
                 "EarthEnvTopoMed_Elevation",
                 "EarthEnvTopoMed_Northness",
                 "EarthEnvTopoMed_ProfileCurvature",
                 "EarthEnvTopoMed_Roughness",
                 "EarthEnvTopoMed_Slope",
                 "EarthEnvTopoMed_TopoPositionIndex",
                 "SG_Absolute_depth_to_bedrock",
                 "WorldClim2_SolarRadiation_AnnualMean",
                 "WorldClim2_WindSpeed_AnnualMean",
                 "WorldClim2_H2OVaporPressure_AnnualMean",
                 "NDVI",
                 "EVI",
                 "Lai",
                 "Fpar",
                 "Npp",
                 "Tree_Density",
                 "PET",
                 "SG_Clay_Content_0_100cm",
                 "SG_Coarse_fragments_0_100cm",
                 "SG_Sand_Content_0_100cm",
                 "SG_Silt_Content_0_100cm",
                 "SG_Soil_pH_H2O_0_100cm",
                 "LandCoverClass_Cultivated_and_Managed_Vegetation",
                 "LandCoverClass_Urban_Builtup",
                 "Human_Disturbance",
                 "PresentTreeCover",
                 "Nitrogen",
                 "cropland",
                 "grazing",
                 "pasture",
                 "rangeland",
                 "cnRatio",
                 "Cation",
                 "SoilMoisture",
                 "ForestAge",
                 "Organic_Carbon",
                 "GlobBiomass_GrowingStockVolume",
                 "Human_Development_Percentage",
                 "HumanFootprint",
                 "Rainfall_Erosivity",
                 #"WDPA",
                 #"Fire_Frequency",
                 "WaterAvailability",
                 'CHELSA_vpd') #"FireFrequency" due to the data type will not be included in this analysis.

# all(retainedVars %in% colnames(trainDataTable))
# sapply(randomSampledTable[, ..retainedVars], function(x) length(unique(x)))
# sort(sapply(randomSampledTable[, ..retainedVars], function(x) length(unique(x))))

# generate the formula for the game model
formulaString = as.formula(paste('MaxHeight','~', paste("s(",retainedVars,")", collapse="+",sep=""), sep=""))
# traint the gam model
gamMulti = gam(data = duplicateTable, formula = formulaString) 
duplicateTable$gamResiduals = gamMulti$residuals
coordinates(duplicateTable) = ~x+y
proj4string(duplicateTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# # calculate the Moran I changes along the spatial distance
spatialMoranDynamics = spline.correlog(x=coordinates(duplicateTable)[,1], y=coordinates(duplicateTable)[,2],z=duplicateTable$gamResiduals, resamp=50, quiet=TRUE,latlon=T,xmax=1000,df=500)
# make the plot
pdf(paste("Plots/Figure_SX_MoranI_along_distance_of_gam_model_residuals_for_MaxHeight_twoDataset_50km.pdf",sep=""),width = 5, height=4)
plot(spatialMoranDynamics,xlim=c(0,1000),ylim=c(-1,1))
abline(v=50,col="red", lty=2)
text(150,0.5,"50km")
dev.off()