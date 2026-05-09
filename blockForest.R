######################
### LOAD LIBRARIES ###
######################


library(mltools)
library(matrixStats)
library(parallel)
library(blockForest)
library(reshape)
library(ggplot2)
library(MLmetrics)
library(caret)


################################
### DEFINE GLOBAL PARAMETERS ###
################################


setwd("/sps/bioaster/Projects/ISYBIO/")
source("/sps/bioaster/Projects/ISYBIO/Benchmark/src/Shared/nestedCV.R")

dataType ="simulations" # "realData" #
nbcores = 40
nRep = 40
nFolds = 5
nodesizes = c(1,3,6)
nRepBlockWeights = 5
seedval<-c(895404, 392954, 953004, 876087, 263713, 193086, 675373, 536070, 189527, 793618,
		833665, 941433, 652399, 704627, 527325, 767149, 873279, 317183, 845440, 722964,
		103341, 214969, 259822, 918061, 589755, 661554, 934839, 838399, 459699, 757984,
		591911, 429887, 927056, 618164, 784879, 641952, 324179, 486233, 475693, 757759)

dataset2keep = list(realData = c("RV144","COVID19","TCGA"),
	simulations = c("Reference","n/5","px5","Case_Control_1:7","High_Main_MO","High_Conf_MO_Overlap","Main_MO_2Smallest_Omics","Main_MO_1Smallest_Omic","High_Fract_Signal_Feat","nx5","p/5","High_Conf_SO_Overlap","Main_MO_1Largest_Omic","Main_MO_2Largest_Omics","Noise")
)


if(dataType=="simulations"){
	load("Data_Simulations/Data/MOFA/final/MultiOmics_Simulations_final.RData")
}else{
	load("Public_MultiOmics_datasets/Data/All_Datasets_Predictive_RmConfEffect.RData")
	MO = MO[realData2Keep]
	names(MO) = gsub("rv144_trt","RV144",names(MO))
	convertViewName = c("Transcripts (mRNA)","Haplotypes","Cytokines","Cell Types","Metabolites","Lipids","Proteins","Transcripts (mRNA)","Transcripts (mRNA)","Transcripts (miRNA)","Proteins")
	names(convertViewName) = c("esetBaselined","haplotypeSet","luminexSet","phenotypeSet","metabolites","lipids","proteins","transcripts","RNASeq2GeneNorm","miRNASeqGene","RPPAArray")
	boxplotCols = c("grey100", "grey90", "grey50")
}


#################################################
################# BLOCK FOREST ##################
#################################################


MCC_F1_blockWeights = lapply(names(MO), function(sim){
	mclapply(1:nRep, function(rep) {

		set.seed(seedval[rep])
		Y = factor(MO[[sim]][[rep]]$Z_covar[,1])
		X = data.frame(do.call(cbind,lapply(MO[[sim]][[rep]]$data,t)))
		if(dataType=="simulations"){colnames(X)<-paste0("Feature_",1:ncol(X))}
		colnames(X) = rmSpecialChar(colnames(X))

		### Define block variables
		blocks = rep(1:length(MO[[sim]][[rep]]$data), times=sapply(MO[[sim]][[rep]]$data,nrow))
		blocks = lapply(1:length(MO[[sim]][[rep]]$data), function(x) which(blocks==x))
		
		### Run nestedCV (tune hyperparameters and generate predictions) 
		mtry = sapply(blocks, function(x) getMtryGrid(length(x)))
		nestedcv = nestedCV(X=X, Y=Y, blocks=blocks, fun="blockForest", nFolds=nFolds, nodesizes=nodesizes, mtry=mtry)

		### Compute MCC/F1 and return  max MCC/F1
		MCCF1 = c(MCC=mcc(preds=nestedcv$predictions,actuals=Y), F1=F1_Score(nestedcv$predictions,Y))
		return(list(blockWeights=nestedcv$blockWeights, MCC_F1_BlockForest = data.frame(Method="BlockForest",Scenario=sim, Rep=rep, t(MCCF1))))
	}, mc.cores=nbcores)
})

MCC_F1_BlockForest = do.call("rbind",lapply(MCC_F1_blockWeights,function(x) do.call("rbind",lapply(x, function(y) y$MCC_F1_BlockForest))))
save(MCC_F1_BlockForest,file=paste0("Benchmark/Data/Predictive/MCC_F1_BlockForest_",dataType,".RData"))


#################################################
### BOXPLOT blockWeights PER SCENARIO/DATASET ###
#################################################


### EXTRACT WEIGHTS FROM MCC_BF OBJECT
blockWeights.df  = do.call("rbind",lapply(1:length(MCC_F1_blockWeights), function(i) data.frame(Scenario=names(MO)[[i]], melt(t(sapply(MCC_F1_blockWeights[[i]], function(y) y$blockWeights))))))
colnames(blockWeights.df) = c("Scenario","Repeat","View","MCC")

### FORMAT blockWeights.df
if(dataType=="simulations"){
	blockWeights.df$View = as.factor(paste0("view_",blockWeights.df$View))
	blockWeights.df$Scenario = gsub("Intersect","Overlap",gsub("_Multiplied_","x",gsub("_Divided_","/",gsub("Equal_p_","p=",gsub("largest","Largest",gsub("smallest","Smallest",blockWeights.df$Scenario))))))
	blockWeights.df = blockWeights.df[blockWeights.df$Scenario %in% dataset2keep[[dataType]],]
	blockWeights.df$Scenario = factor(blockWeights.df$Scenario, levels=dataset2keep[[dataType]])
}else{
	for(sim in names(MO)){
		blockWeights.df[blockWeights.df$Scenario==sim,"View"] = convertViewName[names(MO[[sim]][[1]]$data)[as.numeric(blockWeights.df[blockWeights.df$Scenario==sim,"View"])]]
	}
	blockWeights.df$Scenario = factor(blockWeights.df$Scenario,levels=names(MO))
}

### PLOT WEIGHTS PER SCENARIO/DATASET
p <- ggplot(blockWeights.df)
if(dataType=="simulations"){
	p <- p + geom_boxplot(aes(y=MCC, fill=View, x=Scenario), outlier.shape=NA)
}else{
	p <- p + geom_boxplot(aes(y=MCC, fill=Scenario, x=View), outlier.shape=NA) + 
	scale_fill_manual(values = boxplotCols) + theme(legend.position = "none") + facet_wrap(~Scenario, scales="free") +
	guides(fill="none")
}
p <- p + theme_bw() + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(paste0("/sps/bioaster/Projects/ISYBIO/Benchmark/Figures/Predictive/blockForest_boxplots_weights_",dataType,"2.pdf"), width=10, height=6)

