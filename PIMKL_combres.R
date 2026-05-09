######################
### LOAD LIBRARIES ###
######################

library(reshape2)
library(ggplot2)
library(mltools)
library(MLmetrics)


################################
### DEFINE GLOBAL PARAMETERS ###
################################


setwd("/sps/bioaster/Projects/ISYBIO/")

dataType = "simulations" #"realData" # "simulations_Noise"
nRep = 40
if(dataType=="simulations"){load("Data_Simulations/Data/MOFA/final/MultiOmics_Simulations_final.RData")}
if(dataType=="realData"){
	load("Public_MultiOmics_datasets/Data/All_Datasets_Predictive_RmConfEffect.RData")
	MO = MO[-c(1,3,4)]
}

######### PIMKL 
MCC_F1_PIMKL<-c()
for (sim in names(MO)){
	if(dataType=="realData"){
		directory = paste0("Benchmark/Data/Predictive/PIMKL/",dataType,"/",sim,"/")
		PIMKLres_files = dir(directory)
	}else{
		directory = paste0("Benchmark/Data/Predictive/PIMKL/",dataType,"/output/")
		pattern_scenario=paste0("Y_df_",gsub("_","",sim),"view")
		PIMKLres_files = dir(directory,pattern=pattern_scenario)
	}
	PIMKL_mcc<-c()
	PIMKL_F1<-c()
	for (i in PIMKLres_files)	{
		tdf<-read.csv(paste0(directory,i))
		#print(table(tdf[,2:3]))
		PIMKL_mcc<-c(PIMKL_mcc,mcc(tdf[,2],tdf[,3]))
		PIMKL_F1<-c(PIMKL_F1,F1_Score(tdf[,2],tdf[,3]))
	}
	MCC_F1_PIMKL<-rbind(MCC_F1_PIMKL,data.frame(Method="PIMKL",Scenario=sim, Rep=1:40, MCC=PIMKL_mcc,F1=PIMKL_F1))
}

save(MCC_F1_PIMKL,file=paste0("Benchmark/Data/Predictive/MCC_F1_PIMKL_",dataType,".RData"))
