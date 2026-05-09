######################
### LOAD LIBRARIES ###
######################


library(RSpectra)
library(CVXR)
library(mclust)
library(randomForest)
library(doParallel)
library(mixOmics)
library(dplyr)
library(mltools)
library(SIDA)
library(MLmetrics)
library(caret)


################################
### DEFINE GLOBAL PARAMETERS ###
################################

setwd("/sps/bioaster/Projects/ISYBIO/")
source("/sps/bioaster/Projects/ISYBIO/Benchmark/src/Shared/cv_SIDA.R") #### need to additionally source SIDA help functions
source("/sps/bioaster/Projects/ISYBIO/Benchmark/src/Shared/addSIDAfunc.R")

dataType = "simulations"# "realData" #
nbcores = 40
nRep = 40
wht=1
lambdavec<-0:8/10
nFolds=5

seedvalvec<-c(895414, 392974, 953004, 876087, 263713, 193086, 675373, 536070, 189527, 793618,
		833665, 941433, 652399, 704627, 527325, 767149, 873279, 317183, 845450, 722974,
		103341, 214969, 259822, 918061, 589755, 661554, 934839, 838399, 459699, 757984,
		591911, 429887, 927056, 618164, 784879, 641952, 324179, 486233, 475693, 757759,
		210741, 967381, 603458, 178028, 285096, 962465, 435403, 133645, 288010, 794855)

if(dataType=="simulations"){
	load("Data_Simulations/Data/MOFA/final/MultiOmics_Simulations_final.RData")
	MO = MO[-c(7,8,12,14,16)]
}
if(dataType=="realData"){
	load("Public_MultiOmics_datasets/Data/All_Datasets_Predictive_RmConfEffect.RData")
	MO = MO[-c(1,3,4)]
}

for (sim in names(MO)){
	cat("Sim ",sim,"\n")
	SIDA_res <- mclapply(1:nRep, function(rep) {

		Xdata <- lapply(MO[[sim]][[rep]]$data,function(x) scale(t(x)))
		if(dataType=="realData" & sim=="MMRN"){Xdata$MetaboPos <- Xdata$MetaboPos[,-c(12)]}
		Y <- MO[[sim]][[rep]]$Z_covar[,1]+1
		seedval<-seedvalvec[rep]
		iFolds<-createFolds(Y, k = nFolds)
		predY<-c()
		actY<-c()
		taulist<-c()

		for (ifold in 1:nFolds){
			Xtrain<-lapply(Xdata,function(x) x[-iFolds[[ifold]],])
			Ytrain<-Y[-iFolds[[ifold]]]
			Xtest<-lapply(Xdata,function(x) x[iFolds[[ifold]],])
			Ytest<-Y[iFolds[[ifold]]]
			criter<-0
			while (criter<10){
				mycv2 <- tryCatch({
						cv_SIDA(Xtrain,Ytrain,withCov=FALSE,plotIt=FALSE, Xtestdata=Xtest,Ytest=Ytest,
						isParallel=FALSE,ncores=nbcores,gridMethod='GridSearch',
						AssignClassMethod='Joint',nfolds=5,ngrid=10,standardize=FALSE,
						maxiteration=20, weight=wht,thresh=1e-03,seedval=seedval,lambda=lambdavec)
				},error=function(e) NA)
				if (!is.na(mycv2)){
						criter = criter + 10
					} else {
						criter=criter+1
						seedval = seedval+1
					}
			}
			if (!is.na(mycv2)){
				predY<-c(predY,mycv2$mysida$PredictedClass)
				actY<-c(actY,Ytest)
				taulist<-c(taulist,unlist(mycv2$optTau))
			}
			rm(mycv2)
		}
		return(list(prediction=predY,actuals=actY,taulist=taulist))
	},mc.cores=nbcores)
	save(SIDA_res,file=paste0("Benchmark/Data/Predictive/SIDA/SIDA_",dataType,"_",sim,".RData"))
	rm(SIDA_res)
}

##################################
### COMBINE SCENARIOS/DATASETS ###
##################################

SIDA_res_final<-list()
for (sim in names(MO)){
	load(paste0("Benchmark/Data/Predictive/SIDA/SIDA_",dataType,"_",sim,".RData"))
	SIDA_res_final[[sim]] = SIDA_res
	rm(SIDA_res)
}

SIDA_res<-SIDA_res_final
save(SIDA_res,file=paste0("Benchmark/Data/Predictive/SIDA_",dataType,".RData"))

MCC_SIDA<-lapply(SIDA_res, function(x) lapply(x,function(y) tryCatch({mcc(y$actuals,y$prediction)},error=function(e) NA)))
F1_SIDA<-lapply(SIDA_res, function(x) lapply(x,function(y) tryCatch({F1_Score(y$actuals,y$prediction)},error=function(e) NA)))
MCC_F1_SIDA<-data.frame(Method="SIDA",Scenario=rep(names(MO),each=40), Rep=1:40, MCC=unlist(MCC_SIDA),F1=unlist(F1_SIDA))
save(MCC_F1_SIDA,file=paste0("Benchmark/Data/Predictive/MCC_F1_SIDA_",dataType,".RData"))






