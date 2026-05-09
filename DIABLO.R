######################
### LOAD LIBRARIES ###
######################


library(mclust)
library(randomForest)
library(parallel)
library(mixOmics)
library(dplyr)
library(mltools)
library(reshape2)
library(RColorBrewer)
library(MLmetrics)


################################
### DEFINE GLOBAL PARAMETERS ###
################################


dataType = "simulations" #"realData" # 
nRep = 40
maxComp = 1
nFolds = 5
nSparsity = 8
nrepeat = 3
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))

setwd("/sps/bioaster/Projects/ISYBIO/")
source("/sps/bioaster/Applications/Rscripts/stats_functions.R")
source("/sps/bioaster/Projects/ISYBIO/Benchmark/src/Shared/perfdiablomodif.R") #### need to additionally source some mixOmics functions
if(dataType=="simulations"){load("Data_Simulations/Data/MOFA/final/MultiOmics_Simulations_final.RData")}
if(dataType=="realData"){
	load("Public_MultiOmics_datasets/Data/All_Datasets_Predictive_RmConfEffect.RData")
	MO = MO[-c(1,3,4)]
	names(MO) = gsub("rv144_trt","RV144",names(MO))
	convertViewName = c("Transcripts (mRNA)","Haplotypes","Cytokines","Cell Types","Metabolites","Lipids","Proteins","Transcripts (mRNA)","Transcripts (mRNA)","Transcripts (miRNA)","Proteins")
	names(convertViewName) = c("esetBaselined","haplotypeSet","luminexSet","phenotypeSet","metabolites","lipids","proteins","transcripts","RNASeq2GeneNorm","miRNASeqGene","RPPAArray")
}
scenarioNames = list(
	figure2 = c("Reference","n/5","px5","Case_Control_1:7","High_Main_MO","High_Conf_MO_Overlap","Main_MO_2Smallest_Omics","Main_MO_1Smallest_Omic","High_Fract_Signal_Feat"),
	supplementary = c("nx5","p/5","High_Conf_SO_Overlap","Main_MO_1Largest_Omic","Main_MO_2Largest_Omics","Noise")
)


###############################
### FIT AND EVALUATE DIABLO ###
###############################


DIABLOres = lapply(names(MO), function(sim){
	mclapply(1:nRep, function(rep){

		cat("Sim ",sim,"Nrep= ",rep)
		Xdata = lapply(MO[[sim]][[rep]]$data,t)
		Xdata = lapply(Xdata,scale)
		Y = factor(MO[[sim]][[rep]]$Z_covar[,1]+1)

		#design 
		design = matrix(0, ncol = length(Xdata), nrow = length(Xdata), 
						dimnames = list(names(Xdata), names(Xdata)))

		### generate parameter grid 
		test.keepX = lapply(Xdata, function(x) unique(round(exp(seq(1,floor(log(ncol(x))),length.out=nSparsity)))))

		### find best predictive parameter combination
		tune = tune.block.splsda(X = Xdata, Y = Y, ncomp = maxComp, 
									  test.keepX = test.keepX, design = design,
									  validation = 'Mfold', folds = nFolds, nrepeat = nrepeat, dist = "mahalanobis.dist")

		if(is.null(ncomp <- tune$choice.ncomp$ncomp)){ncomp = maxComp}
		list.keepX = lapply(tune$choice.keepX, function(x) x[ncomp])

		### fit model
		final.diablo.model = block.splsda(X = Xdata, Y = Y, ncomp = ncomp, keepX = list.keepX, design = design)
							  
		### compute performances
		perf = perf.sgccdamodif(final.diablo.model, validation = 'Mfold', folds = nFolds, nrepeat = nrepeat, dist = "mahalanobis.dist")
		weights = as.data.frame(perf$weights %>% group_by(block) %>% summarize(avg = mean(comp1)))$avg
		weightedvote = lapply(perf$WeightedVote,unlist)
		prediction = apply(do.call(cbind,weightedvote),1,function(x) names(which.max(table(x))))
		return(list(prediction=prediction, ncomp=ncomp, weights=weights, nVarSelected=unlist(list.keepX)))
	}, mc.cores=nRep)
})

names(DIABLOres) = names(MO)
save(DIABLOres,file=paste0("Benchmark/Data/Predictive/predictions_DIABLO_",dataType,".RData"))


##################################################
### PLOT WEIGHTS DISTRIBUTION ACROSS SCENARIOS ###
##################################################


load(paste0("Benchmark/Data/Predictive/predictions_DIABLO_",dataType,".RData"))

### EXTRACT WEIGHTS AND NUMBER OF SELECTED VARIABLES
weights.nVarSelected = NULL
for(var in c("weights","nVarSelected")){
	tmp = lapply(1:length(DIABLOres), function(iSim){
		data.frame(Scenario=names(MO)[iSim], Variable=var, do.call("rbind",lapply(DIABLOres[[iSim]], function(y) y[[var]])))
	})
	ncol = sapply(tmp,ncol)
	for(i in 1:length(tmp)){
		if(ncol[i]<max(ncol)){
			tmp[[i]] = data.frame(tmp[[i]],rep(NA,nrow(tmp[[i]])))
		}
		colnames(tmp[[i]])[-(1:2)]=paste0("view_",1:(max(ncol)-2))
	}
	weights.nVarSelected = rbind(weights.nVarSelected, melt(do.call("rbind",tmp),variable.name = "View"))
}

### RENAME VARIABLES, FILTER nVarSelected AND Scenario
weights.nVarSelected$View = gsub("X","view_",weights.nVarSelected$View)
weights.nVarSelected$Variable = gsub("weights","Weights",gsub("nVarSelected","Number of Selected Variables",weights.nVarSelected$Variable))

### PLOT WEIGHTS AND NUMBER OF SELECTED VARIABLES
if(dataType=="simulations"){
	weights.nVarSelected$Scenario = gsub("Intersect","Overlap",gsub("_Multiplied_","x",gsub("_Divided_","/",gsub("Equal_p_","p=",gsub("largest","Largest",gsub("smallest","Smallest",weights.nVarSelected$Scenario))))))
	weights.nVarSelected = weights.nVarSelected[weights.nVarSelected$Scenario %in% unlist(scenarioNames),]
	weights.nVarSelected$Scenario = factor(weights.nVarSelected$Scenario, levels=unlist(scenarioNames))
	weights.nVarSelected = weights.nVarSelected[weights.nVarSelected$value<400,]
	p <- ggplot(weights.nVarSelected, aes(y=value, fill=View, x=Scenario)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme() +  xlab("") +
		facet_grid(Variable~., scales="free") + ylab("") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 
	ggsave(paste0("/sps/bioaster/Projects/ISYBIO/Benchmark/Figures/Predictive/DIABLO_boxplots_weights_nVarSelected_",dataType,"2.pdf"), width=10, height=6)

}else{
	for(sim in names(MO)){
		weights.nVarSelected[weights.nVarSelected$Scenario==sim,"View"] = convertViewName[names(MO[[sim]][[1]]$data)[as.numeric(gsub("view_","",weights.nVarSelected[weights.nVarSelected$Scenario==sim,"View"]))]]
	}
	# weights.nVarSelected[,"View"] = factor(weights.nVarSelected[,"View"], levels=convertViewName)
	weights.nVarSelected[,"Scenario"] = factor(weights.nVarSelected[,"Scenario"], levels=names(MO)[c(1,3:2)])
	weights.nVarSelected = weights.nVarSelected[!is.na(weights.nVarSelected$value),]
	p <- ggplot(weights.nVarSelected, aes(y=value, fill=View, x=View)) + geom_boxplot(outlier.shape=NA) + theme_bw() + theme(legend.title = element_blank()) +  xlab("") +
		facet_grid(Variable~Scenario, scales="free") + ylab("") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  theme(legend.position = "none")
	ggsave(paste0("/sps/bioaster/Projects/ISYBIO/Benchmark/Figures/Predictive/DIABLO_boxplots_weights_nVarSelected_",dataType,"2.pdf"), width=10, height=8)
}


##########################
### COMPUTE MCC AND F1 ###
##########################


### COMPUTE MCC AND F1
MCC_F1_DIABLO = do.call("rbind",lapply(names(MO), function(sim){
	do.call("rbind",lapply(1:length(MO[[sim]]),function(rep){
		data.frame(Method="DIABLO", Scenario=sim, Rep=rep,
		MCC = mcc(preds=as.numeric(DIABLOres[[sim]][[rep]]$prediction)-1, actuals=MO[[sim]][[rep]]$Z_covar[,"main_MO"]),
		F1 = F1_Score(as.numeric(DIABLOres[[sim]][[rep]]$prediction)-1, MO[[sim]][[rep]]$Z_covar[,"main_MO"]))
	}))
}))

### SAVE MCC AND F1
save(MCC_F1_DIABLO,file=paste0("Benchmark/Data/Predictive/MCC_F1_DIABLO_",dataType,".RData"))


##########################
### PLOT MCC AND NCOMP ###
##########################


NCOMP = sapply(DIABLO, function(x) sapply(x,function(y) y$ncomp))
MCC.df = data.frame("DIABLO",melt(MCC)[,c(2,1,3)])
NCOMP.df = data.frame("DIABLO",melt(NCOMP)[,c(2,1,3)])
colnames(MCC.df) = c("Method","Scenario","Rep","MCC")
colnames(NCOMP.df) = c("Method","Scenario","Rep","nComp")

pdf(paste0("Benchmark/MCC_DIABLO_",dataType,".pdf"), width=12)
ggplot(MCC.df, aes(y=MCC,x=Scenario, fill=Method)) + geom_boxplot(outlier.shape = NA) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("")
dev.off()

pdf(paste0("Benchmark/Barplot_ncomp_DIABLO.pdf"), width=9)
ggplot(NCOMP.df) + geom_bar(aes(x=nComp, fill=Scenario),position=position_dodge())+ scale_fill_manual(values = getPalette(16)) + theme_bw() + xlab("Number of Components") + ylab("Number of simulations")
dev.off()
