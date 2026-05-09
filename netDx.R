######################
### LOAD LIBRARIES ###
######################


library(mltools)
library(MultiAssayExperiment)
library(doParallel)
library(BiocFileCache)
library(ROCR)
library(parallel)
library(glmnet)
library(MLmetrics)



################################
### DEFINE GLOBAL PARAMETERS ###
################################

setwd("/sps/bioaster/Projects/ISYBIO/")
source("/sps/bioaster/Projects/ISYBIO/Benchmark/src/Shared/netDxfunc.R")
dataType = "simulations" # "realData" #
nbcores = 40
nRep = 40
seedval<-c(895404, 392954, 953004, 876087, 263713, 193086, 675373, 536070, 189527, 793628,
		833665, 941433, 652399, 704627, 527325, 767149, 873279, 317183, 845440, 722964,
		103341, 214969, 259822, 918061, 589755, 661554, 934839, 838399, 459699, 757984,
		591911, 429887, 927056, 618164, 784879, 641952, 324179, 486233, 475693, 757759)


if(dataType=="simulations"){load("Data_Simulations/Data/MOFA/final/MultiOmics_Simulations_final.RData")}
if(dataType=="realData"){
	load("Public_MultiOmics_datasets/Data/All_Datasets_Predictive_RmConfEffect.RData")
	MO = MO[-c(1,3,4)]
}

MCC_netDx = netDx_res = list()

for (sim in names(MO)){
	netDx_res = list()
	for(rep in 1:nRep){
	
		cat("Sim ",sim,"Nrep= ",rep)
		#netDx_res[[sim]][[rep]]<-NA
		
		Xdata <- MO[[sim]][[rep]]$data
		Xdata<-lapply(Xdata,function(x) t(scale(t(x),T,T)))
		if(dataType=="realData"){for (i in 1:length(Xdata)) {colnames(Xdata[[i]])<-paste0("Sample",colnames(Xdata[[i]]))}}   ####important for netDx
		
		Y <- as.character(MO[[sim]][[rep]]$Z_covar[,1])
		patient.data <- data.frame(patiendID=colnames(Xdata[[1]]),ID=colnames(Xdata[[1]]),STATUS=Y, row.names=colnames(Xdata[[1]]))
		
		simMOmultiassay <- MultiAssayExperiment(experiments = Xdata,colData = patient.data)

		groupliste<-list()
		for (i in 1:length(Xdata))
		{
		tmp <- list(rownames(experiments(simMOmultiassay)[[i]]));
		names(tmp) <- names(simMOmultiassay)[i]
		groupliste[[names(simMOmultiassay)[[i]]]] <- tmp
		}
		
		sims <- rep(list("pearsonCorr"),length(Xdata))
		
		# map layer names to sims
		names(sims) <- names(groupliste)
		
		criter=0
		while (criter<10)
		{
			outDir <- paste(tempdir(),paste0("resMO/",sim,"_Rep_",rep),sep=getFileSep()) # use absolute path
			if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
			
			t0 <- Sys.time()
			resrep <- tryCatch(suppressMessages(
			buildPredictor_modif(
				dataList=simMOmultiassay,          ## your data
				groupList=groupliste,    ## grouping strategy
				sims=sims,
				outDir=outDir,          ## output directory
				trainProp=0.8,          ## pct of samples to use to train model in
				preFilter=TRUE,
				## each split
				numSplits=5L,            ## number of train/test splits
				featSelCutoff=9L,       ## threshold for calling something
				## feature-selected
				featScoreMax=10L,    ## max score for feature selection
				numCores=1,            ## set higher for parallelizing
				debugMode=FALSE,
				keepAllData=FALSE,     ## set to TRUE for debugging
				logging="default",     ## set to "default" for messages
				seedval=seedval[rep]
			)
			),error=function(e) NA)
			t1 <- Sys.time()
			print(t1-t0)
			if (!is.na(resrep)) 
			{
				criter = criter+10 
			} else {
			criter=criter+1
			seedval = seedval+1}
		}

		if (!is.na(resrep)) {
			dfresrep<-rbind(resrep$Split1$predictions,resrep$Split2$predictions,
							resrep$Split3$predictions,resrep$Split4$predictions,
							resrep$Split5$predictions)
			mcc<-mcc(dfresrep$PRED_CLASS,dfresrep$STATUS)
			netDx_res[[rep]] = list(mcc=mcc, prediction=dfresrep, resrep=resrep)
		}
	} ###mc.cores should be equal to 1 because of parallelization libraries incompatibility
	save(netDx_res,file=paste0("/sps/bioaster/Projects/ISYBIO/Benchmark/Data/Predictive/netDx/netDx_",dataType,"_",sim,".RData"))
}


##################################
### COMBINE SCENARIOS/DATASETS ###
##################################

netDx_res_final<-list()
for (sim in names(MO)){
	load(paste0("Benchmark/Data/Predictive/netDx/netDx_",dataType,"_",sim,".RData"))
	print(sim)
	print(summary(sapply(netDx_res, function(x) x[[1]])))
	netDx_res_final[[sim]] = netDx_res
	rm(netDx_res)
}
netDx_res<-netDx_res_final
save(netDx_res,file=paste0("Benchmark/Data/Predictive/netDx_",dataType,".RData"))

#names(netDx_res)<-names(MO)
MCC_netDx<-lapply(netDx_res, function(x) lapply(x,function(y) ifelse(is.null(y[[1]]),NA,y[[1]])))
F1_netDx<-c()
for (sim in names(MO)){
	for(rep in 1:nRep) {
		F1_netDx<-c(F1_netDx,tryCatch({F1_Score(netDx_res[[sim]][[rep]]$prediction$STATUS,netDx_res[[sim]][[rep]]$prediction$PRED_CLASS)},error=function(e) NA))


	}
}
MCC_F1_netDx<-data.frame(Method="netDx",Scenario=rep(names(MO),each=40), Rep=1:40, MCC=unlist(MCC_netDx),F1=F1_netDx)

save(MCC_F1_netDx,file=paste0("Benchmark/Data/Predictive/MCC_F1_netDx_",dataType,".RData"))


