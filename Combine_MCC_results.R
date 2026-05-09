####################
### LOAD LIBRARY ###
####################


library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(randomForest)
library(mltools)
library(parallel)
library(reshape)
library(MLmetrics)
library(caret)


################################
### DEFINE GLOBAL PARAMETERS ###
################################

setwd("/sps/bioaster/Projects/ISYBIO/")
source("/sps/bioaster/Projects/ISYBIO/Benchmark/src/Shared/nestedCV.R")

dataType = "realData" #"simulations" #
supplementary = F
nbcores = 40
nRep = 40
nFolds = 5
nodesizes = c(1,3,6)
methodNames = list(RData = c("RF_sep","RF_concat","BlockForest","PIMKL","netDx","SIDA","DIABLO","stack_generalization_SL"),
					manuscript = c("DIABLO", "SIDA", "PIMKL", "netDx","Stacked_Generalization","BlockForest","RF_Concat", "RF_Max_Single_View"))
convertViewName = c("Transcripts (mRNA)","Haplotypes","Cytokines","Cell Types","Metabolites","Lipids","Proteins","Transcripts (mRNA)","Transcripts (mRNA)","Transcripts (miRNA)","Proteins")
names(convertViewName) = c("esetBaselined","haplotypeSet","luminexSet","phenotypeSet","metabolites","lipids","proteins","transcripts","RNASeq2GeneNorm","miRNASeqGene","RPPAArray")
datasetNames = list()
datasetNames[["simulations"]] = list(
	figure2 = c("Reference","n/5","px5","Case_Control_1:7","High_Main_MO","High_Conf_MO_Overlap","Main_MO_2Smallest_Omics","Main_MO_1Largest_Omic","High_Fract_Signal_Feat"), #"p=240","High_Conf_MO",
	supplementary = c("nx5","p/5","High_Conf_SO_Overlap","Main_MO_1Smallest_Omic","Main_MO_2Largest_Omics","Noise") #"High_Conf_SO",
)
datasetNames[["realData"]] = c("RV144","COVID19","TCGA")

boxplotCols = list(integrative = c("#B2182B","#FF7F00","#6A3D9A","#FFD92F","#1F78B4","#006D2C","#74C476","#E5F5E0"), nonIntegrative = c("grey100", "grey90", "grey50"))

seedval<-c(895404, 392954, 953004, 876087, 263713, 193086, 675373, 536070, 189527, 793618,
		833665, 941433, 652399, 704627, 527325, 767149, 873279, 317183, 845440, 722964,
		103341, 214969, 259822, 918061, 589755, 661554, 934839, 838399, 459699, 757984,
		591911, 429887, 927056, 618164, 784879, 641952, 324179, 486233, 475693, 757759)


#######################################################
### COMPUTE MCC ON EACH DATASET SEPARATELY USING RF ###
#######################################################


if(dataType=="realData"){

	try(load(paste0("Benchmark/Data/Predictive/MCC_F1_RF_single_",dataType,".RData")))
	if(!exists("MCC_F1_RF_single")){

		load("Public_MultiOmics_datasets/Data/All_Datasets_Predictive_RmConfEffect.RData")
		names(MO) = gsub("rv144_trt","RV144",names(MO))
		nVar = unlist(sapply(datasetNames[[dataType]], function(sim) sapply(MO[[sim]][[1]]$data, nrow)))
		names(nVar) = gsub("[^\\.]+\\.(.*)","\\1",names(nVar))

		MCC_F1_RF_single = do.call("rbind",lapply(datasetNames[[dataType]], function(sim){
			
			do.call("rbind",mclapply(1:nRep, function(rep){
				
				Y = factor(MO[[sim]][[rep]]$Z_covar[,1])
				
				do.call("rbind",lapply(names(MO[[sim]][[rep]]$data), function(omicName){
					set.seed(seedval[rep])
					X = t(MO[[sim]][[rep]]$data[[omicName]])
					colnames(X) =  rmSpecialChar(colnames(X))

					### Run nestedCV (tune hyperparameters and generate predictions) 
					mtry = getMtryGrid(ncol(X))
					nestedcv = nestedCV(X=X, Y=Y, funName="randomForest", nFolds=nFolds, nodesizes=nodesizes, mtry=mtry)

					### Compute MCC/F1 and return results
					MCCF1 = c(MCC=mcc(preds=nestedcv$predictions,actuals=Y), F1=F1_Score(nestedcv$predictions,Y))
					return(data.frame(Method="RF_single_view", omicName=omicName, Scenario=sim, Rep=rep, t(MCCF1)))
				}))
			}, mc.cores=nbcores))
		}))

		MCC_F1_RF_single$omicName = sapply(MCC_F1_RF_single$omicName, function(omic) paste0(convertViewName[omic]," (",nVar[omic],")"))
		MCC_F1_RF_single$Scenario = factor(MCC_F1_RF_single$Scenario, levels=datasetNames[[dataType]])

		save(MCC_F1_RF_single,file=paste0("Benchmark/Data/Predictive/MCC_F1_RF_single_",dataType,".RData"))
	}
}


###########################
### COMBINE MCC RESULTS ###
###########################


MCC_all = NULL
for(method in methodNames[["RData"]]){
	try(load(paste0("Benchmark/Data/Predictive/MCC_F1_",method,"_",dataType,".RData")))
	try(MCC_all <- rbind(MCC_all,get(paste0("MCC_F1_",method))))
}

if (dataType=="realData"){
	MCC_all$Scenario = sub("rv144_trt","RV144",MCC_all$Scenario)
	MCC_all$Scenario = factor(MCC_all$Scenario,levels=datasetNames[[dataType]])
}

if (dataType=="simulations"){
	MCC_all$Scenario = gsub("Intersect","Overlap",gsub("_Multiplied_","x",gsub("_Divided_","/",gsub("Equal_p_","p=",MCC_all$Scenario))))
	MCC_all = MCC_all[MCC_all$Scenario %in% datasetNames[[dataType]][[1+supplementary]],]
	MCC_all$Scenario = factor(MCC_all$Scenario, levels=datasetNames[[dataType]][[1+supplementary]])
}

MCC_all$Method=factor(gsub("Stack","Stacked",gsub("concat","Concat",gsub("max_single_view","Max_Single_View",MCC_all$Method))), levels=methodNames[["manuscript"]])


#################################
### PLOT MCC/F1 SCORE RESULTS ###
#################################


for(metric in c("MCC","F1")){

	### BOXPLOTS MCC/F1 FOR INTEGRATIVE METHODS
	p1<-ggplot(MCC_all, aes(y=get(metric),x=Scenario,fill=Method)) + geom_boxplot(outlier.shape = NA) + theme_bw() +
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + ylab(metric) +
		scale_fill_manual(values = boxplotCols[["integrative"]]) + coord_cartesian(ylim=c(ifelse(dataType=="simulations",ifelse(metric=="MCC",-0.25,0.2),0),1))
	if(dataType=="simulations"){
		ggsave(filename=paste0("Benchmark/Figures/Predictive/Boxplots_",metric,"_allMethods_",dataType,ifelse(supplementary,"_Supplementary",""),"2.pdf"),plot=p1,width=ifelse(dataType=="simulations",15,10))
	}

	if (dataType=="realData"){
		### BOXPLOTS MCC/F1 FOR RF ON SEPARATE DATASETS
		p2<-ggplot(MCC_F1_RF_single, aes(y=get(metric),x=omicName,fill=Scenario)) + geom_boxplot(outlier.shape = NA) + theme_bw() + 
		theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + ylab(metric) + ylim(-0.1,1)+
		facet_grid(.~Scenario, scales = "free", space = "free") + scale_fill_manual(values=boxplotCols[["nonIntegrative"]]) + guides(fill="none")		

		### COMBINE THE TWO PLOTS
		pdf(paste0("Benchmark/Figures/Predictive/Boxplots_",metric,"_allMethods_",dataType,".pdf"),width=15,height=6)
		grid.arrange(p2 + labs(tag = "A") + theme(plot.tag.position = c(0, 0.98),plot.tag = element_text(face="bold")),
					p1 + labs(tag = "B") + theme(plot.tag.position = c(0, 0.98),plot.tag = element_text(face="bold"))+
					scale_fill_manual(labels=c("DIABLO", "SIDA", "PIMKL", "netDx","Stacked_Generalization","BlockForest","RF_Concat", "RF_Max_Single_View"),
					values = boxplotCols[["integrative"]]),widths=c(1.2,2))
		dev.off()
	}
}


#################################
### COMPARE METHODS ON MEDIAN ###
#################################


medianScenario = sapply(unlist(datasetNames[[dataType]],use.names=F), function(scenario) aggregate(MCC~Method,FUN=median, MCC_all[MCC_all$Scenario==scenario,])[,2])
deltaMedianScenario = medianScenario[,-1] - medianScenario[,1]
rownames(medianScenario) = rownames(deltaMedianScenario) = names(table(MCC_all$Method))
