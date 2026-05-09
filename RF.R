######################
### LOAD LIBRARIES ###
######################


library(randomForest)
library(mltools)
library(matrixStats)
library(parallel)
library(MLmetrics)
library(caret)


################################
### DEFINE GLOBAL PARAMETERS ###
################################

setwd("/sps/bioaster/Projects/ISYBIO/")
source("/sps/bioaster/Projects/ISYBIO/Benchmark/src/Shared/nestedCV.R")

dataType ="simulations" #"realData" # 
realData2Keep = c("rv144_trt","COVID19","TCGA")
nbcores = 40
nRep = 40
nFolds = 5
nodesizes = c(1,3,6)
seedval<-c(895404, 392954, 953004, 876087, 263713, 193086, 675373, 536070, 189527, 793618,
		833665, 941433, 652399, 704627, 527325, 767149, 873279, 317183, 845440, 722964,
		103341, 214969, 259822, 918061, 589755, 661554, 934839, 838399, 459699, 757984,
		591911, 429887, 927056, 618164, 784879, 641952, 324179, 486233, 475693, 757759)

MCC_RF_sep = list()
MCC_RF_concat = list()

if(dataType=="simulations"){load("Data_Simulations/Data/MOFA/final/MultiOmics_Simulations_final.RData")}
if(dataType=="realData"){
	load("Public_MultiOmics_datasets/Data/All_Datasets_Predictive_RmConfEffect.RData")
	MO = MO[realData2Keep]
	names(MO) = gsub("rv144_trt","RV144",names(MO))
}


##############################################
################# RF_concat ##################
##############################################


MCC_F1_RF_concat = do.call("rbind",lapply(names(MO), function(sim){
	tmp=do.call("rbind",mclapply(1:nRep, function(rep) {

		set.seed(seedval[rep])
		Y = factor(MO[[sim]][[rep]]$Z_covar[,1])
		X = data.frame(do.call(cbind,lapply(MO[[sim]][[rep]]$data,t)))
		if(dataType=="simulations"){colnames(X)<-paste0("Feature_",1:ncol(X))}
		colnames(X) =  rmSpecialChar(colnames(X))
		
		### Run nestedCV (tune hyperparameters and generate predictions) 
		mtry = getMtryGrid(ncol(X))
		nestedcv = nestedCV(X=X, Y=Y, funName="randomForest", nFolds=nFolds, nodesizes=nodesizes, mtry=mtry)

		# ### Compute MCC/F1 and return results
		# MCCF1 = c(MCC=mcc(preds=nestedcv$predictions,actuals=Y), F1=F1_Score(nestedcv$predictions,Y))
		# return(data.frame(Method="RF_concat",Scenario=sim, Rep=rep, t(MCCF1)))

	}, mc.cores=nbcores))
}))

save(MCC_F1_RF_concat,file=paste0("Benchmark/Data/Predictive/MCC_F1_RF_concat_",dataType,"2.RData"))


##############################################
################### RF_sep ###################
##############################################


MCC_F1_RF_sep = do.call("rbind",lapply(names(MO), function(sim){
	do.call("rbind",mclapply(1:nRep, function(rep){
		
		Y <- factor(MO[[sim]][[rep]]$Z_covar[,1])
		
		MCCF1_View = t(sapply(MO[[sim]][[rep]]$data, function(X){
			
			set.seed(seedval[rep])
			rownames(X) = rmSpecialChar(rownames(X))

			### Run nestedCVRandomForest
			mtry = getMtryGrid(nrow(X))
			nestedcv = nestedCV(X=t(X), Y=Y, funName="randomForest", nFolds=nFolds, nodesizes=nodesizes, mtry=mtry)

			### Compute MCC/F1 and return results
			MCCF1 = c(MCC=mcc(preds=nestedcv$predictions,actuals=Y), F1=F1_Score(nestedcv$predictions,Y))
			return(MCCF1)
		}))

		### Keep max MCC/F1 across all views
		return(data.frame(Method="RF_max_single_view",Scenario=sim, Rep=rep, t(MCCF1_View[which.max(rowSums(MCCF1_View)),])))

	}, mc.cores=nbcores))
}))

save(MCC_F1_RF_sep,file=paste0("Benchmark/Data/Predictive/MCC_F1_RF_sep_",dataType,"2.RData"))

