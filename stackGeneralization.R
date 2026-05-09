######################
### LOAD LIBRARIES ###
######################


library(mltools)
library(reshape2)
library(parallel)
library(SuperLearner)
library(randomForest)
library(MLmetrics)
library(caret)


################################
### DEFINE GLOBAL PARAMETERS ###
################################


dataType ="simulations" #realData" #  
realData2Keep = c("rv144_trt","COVID19","TCGA")
nodesizes = c(1,3,6)
nbcores = 40
nFolds = 5


setwd("/sps/bioaster/Projects/ISYBIO/")
source("/sps/bioaster/Projects/ISYBIO/Benchmark/src/Shared/nestedCV.R")
if(dataType=="simulations"){load("Data_Simulations/Data/MOFA/final/MultiOmics_Simulations_final.RData")}
if(dataType=="realData"){
	load("Public_MultiOmics_datasets/Data/All_Datasets_Predictive_RmConfEffect.RData")
	MO = MO[realData2Keep]
	names(MO) = gsub("rv144_trt","RV144",names(MO))
}
method <- get("method.NNLS", mode = "function")()


#####################################
### FIT AND EVALUATE SUPERLEARNER ###
#####################################


predictions = lapply(names(MO), function(sim){
	
	mclapply(1:length(MO[[sim]]), function(rep){
		
		nView = length(MO[[sim]][[rep]]$data)
		Y = MO[[sim]][[rep]]$Z_covar[,"main_MO"]
		Y = Y - min(Y)
		iFolds = createFolds(Y, k = nFolds)
		level1Predictions = rep(NA,length(Y))
		
		for(iFold in 1:nFolds){
            
			level0Predictions = matrix(NA,nrow=length(Y),ncol=nView)
			ytrain = Y[-iFolds[[iFold]]]
			
            for(iView in 1:nView){

				xtrain = t(MO[[sim]][[rep]]$data[[iView]][,-iFolds[[iFold]]])
                xtest = t(MO[[sim]][[rep]]$data[[iView]][,iFolds[[iFold]]])
				colnames(xtrain) =colnames(xtest) = rmSpecialChar(colnames(xtrain))

                ### Run nestedCV to determine optimal hyperparameters
                mtry = getMtryGrid(ncol(xtrain))
                nestedcv = nestedCV(X=xtrain, Y=as.factor(ytrain), funName="randomForest", nFolds=nFolds, nodesizes=nodesizes, mtry=mtry)
				
				### Estimate RF parameters on nestedcv$finalGridParam to predict test fold
                rf = randomForest(x=xtrain, y=as.factor(ytrain), nodesize=nestedcv$finalGridParam["nodesize"], mtry=nestedcv$finalGridParam["mtry"])

				level0Predictions[-iFolds[[iFold]],iView] = nestedcv$predictions
				level0Predictions[iFolds[[iFold]],iView] = predict(rf, newdata=xtest, type="vote")[,2]
			}
			coef = method$computeCoef(Y=ytrain, Z=level0Predictions[-iFolds[[iFold]],],libraryNames=NA, verbose=F, obsWeights = rep(1, length(unlist(iFolds[-iFold]))))$coef
			level1Predictions[iFolds[[iFold]]] = round(level0Predictions[iFolds[[iFold]],]%*%coef)
		}
		return(level1Predictions)
	}, mc.cores=nbcores)
})

names(predictions) = names(MO)
save(predictions,file=paste0("Benchmark/Data/Predictive/predictions_stack_generalization_SL_",dataType,"2.RData"))


###################
### COMPUTE MCC ###
###################


### COMPUTE MCC AND F1
MCC_F1_stack_generalization_SL = do.call("rbind",lapply(names(MO), function(sim){
	do.call("rbind",lapply(1:length(MO[[sim]]),function(rep){
		data.frame(Method="Stack_Generalization", Scenario=sim, Rep=rep,
		MCC = tryCatch(mcc(preds=predictions[[sim]][[rep]], actuals=MO[[sim]][[rep]]$Z_covar[,"main_MO"]),error=function(e) NA),
		F1 = tryCatch(F1_Score(predictions[[sim]][[rep]], MO[[sim]][[rep]]$Z_covar[,"main_MO"]),error=function(e) NA))
	}))
}))

### SAVE MCC AND F1
save(MCC_F1_stack_generalization_SL,file=paste0("Benchmark/Data/Predictive/MCC_F1_stack_generalization_SL_",dataType,"2.RData"))


############################################
### COMPARE SUPERLEARNER VS RANDOMFOREST ###
############################################


load(paste0("Benchmark/Data/MCC_stack_generalization_SL",dataType,".RData"))
MCC.df = MCC_stackG
MCC.df$Method = sub("Stack_Generalization","RF_glmnet_SVM",MCC.df$Method)
load("/sps/bioaster/pt6/Projets/ISYBIO/Benchmark/Data/MCC_stack_generalization_RF.RData")
MCC.df = rbind(MCC.df,MCC_stackG)
MCC.df$Method = sub("Stack_Generalization","RF",MCC.df$Method)

pdf("Benchmark/MCC_stack_generalization.pdf", width=12)
ggplot(MCC.df, aes(y=MCC,x=Scenario, fill=Method)) + geom_boxplot(outlier.shape = NA) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("")
dev.off()
