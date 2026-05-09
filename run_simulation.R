######################
### LOAD LIBRARIES ###
######################


library(gplots)
library(parallel)
library(doMC)
registerDoMC()
library(randomForest)
library(RColorBrewer)
require(pROC)
library(devtools)
library(limma)
library(ggplot2)
library(metap)
library(circlize)
library(matrixStats)
library(reshape2)
library(mixOmics)
library(Jmisc)


############################################
### SET WD AND DEFINE GENERAL PARAMETERS ###
############################################


setwd("/sps/bioaster/Projects/ISYBIO/Data_Simulations/")
source("src/make_example.R")
source("/sps/bioaster/Applications/Rscripts/stats_functions.R")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

nbcores=100
n_rep = 40
n_views = 3
n_samples = 40
BWRpalette = colorRampPalette(colors=c("blue","yellow"))
palSpectral = colorRampPalette(brewer.pal(11,"Spectral"))


####################################
### DEFINE SIMULATION PARAMETERS ###
####################################


### SCENARII FOR BENCHMARK OF INTEGRATIVE METHODS 
# Reference scenario
alpha = list(high=0.01, medium=0.1, low=1)
scenarii = list()
scenarii[["Reference"]] = list(tau=1,
							n_samples=rep(n_samples,2),
							n_features=c(1000,240,60),
							percentageSharedFeat = 0.5,
							DisableMultiomicsFromView = rep(F,3),
							theta_w=rep(0.1,3),
							alpha_w=rep(alpha[["medium"]],3),
							alpha_z=list(main_MO = unlist(alpha[c("medium","low")]), conf_MO = unlist(alpha[c("medium","low")]),
										DS = list(DS1=unlist(alpha[c("medium","low")]),DS2=unlist(alpha[c("medium","low")]),DS3=unlist(alpha[c("medium","low")]))))


### Define alternative scenarii
alternativeScenariiNames = c("p_Divided_5","p_Multiplied_5","n_Divided_5","n_Multiplied_5","Case_Control_1:7","Equal_p_60","Equal_p_240","Main_MO_2Smallest_Omics",
							"Main_MO_2Largest_Omics","High_Main_MO","High_Conf_MO","High_Conf_MO_Intersect","High_Conf_SO","High_Conf_SO_Intersect",
							"Overlap_Across_Effects","High_Fract_Signal_Feat","Noise","Main_MO_1Smallest_Omic","Main_MO_1Largest_Omic")
for(sim in alternativeScenariiNames){scenarii[[sim]] = scenarii[["Reference"]]}

### Dimensionalty
scenarii[["p_Divided_5"]]$n_features = scenarii[["Reference"]]$n_features/5 #/2
scenarii[["p_Multiplied_5"]]$n_features = scenarii[["Reference"]]$n_features*5 #/2
scenarii[["n_Divided_5"]]$n_samples = scenarii[["Reference"]]$n_samples/5 #/2
scenarii[["n_Multiplied_5"]]$n_samples = scenarii[["Reference"]]$n_samples*5 #/2

### Unbalanced case-control groups
scenarii[["Case_Control_1:7"]]$n_samples = scenarii[["Reference"]]$n_samples*c(1/4,7/4)
scenarii[["Equal_p_60"]]$n_features = rep(60,3)
scenarii[["Equal_p_240"]]$n_features = rep(240,3)

### Remove multi-omics effect from views
scenarii[["Main_MO_2Smallest_Omics"]]$DisableMultiomicsFromView = c(T,F,F)
scenarii[["Main_MO_2Largest_Omics"]]$DisableMultiomicsFromView = c(F,F,T)
scenarii[["Main_MO_1Smallest_Omic"]]$DisableMultiomicsFromView = c(T,T,F)
scenarii[["Main_MO_1Largest_Omic"]]$DisableMultiomicsFromView = c(F,T,T)

### Signal/Noise
scenarii[["High_Main_MO"]]$alpha_z = list(main_MO=unlist(alpha[c("high","low")]), conf_MO=unlist(alpha[c("medium","low")]), DS=list(DS1=unlist(alpha[c("medium","low")]),DS2=unlist(alpha[c("medium","low")]),DS3=unlist(alpha[c("medium","low")])))
scenarii[["High_Conf_MO_Intersect"]]$alpha_z = scenarii[["High_Conf_MO"]]$alpha_z = list(main_MO=unlist(alpha[c("medium","low")]), conf_MO=unlist(alpha[c("high","low")]), DS=list(DS1=unlist(alpha[c("medium","low")]),DS2=unlist(alpha[c("medium","low")]),DS3=unlist(alpha[c("medium","low")])))
scenarii[["High_Conf_SO_Intersect"]]$alpha_z = scenarii[["High_Conf_SO"]]$alpha_z = list(main_MO=unlist(alpha[c("medium","low")]), conf_MO=unlist(alpha[c("medium","low")]), DS=list(DS1=unlist(alpha[c("high","low")]),DS2=unlist(alpha[c("high","low")]),DS3=unlist(alpha[c("high","low")])))
scenarii[["High_Conf_MO_Intersect"]]$percentageSharedFeat = 0.95
scenarii[["High_Conf_SO_Intersect"]]$percentageSharedFeat = 0.95

### Intersection across effects
scenarii[["Overlap_Across_Effects"]]$percentageSharedFeat = 0.95
scenarii[["High_Fract_Signal_Feat"]]$theta_w = rep(0.3,3)

### Noise
scenarii[["Noise"]]$theta_w = rep(0,3)
scenarii[["Noise"]]$alpha_w = rep(10,3)
scenarii[["Noise"]]$alpha_z = list(main_MO=unlist(alpha[c("low","low")]), conf_MO=unlist(alpha[c("low","low")]), DS=list(DS1=unlist(alpha[c("low","low")]),DS2=unlist(alpha[c("low","low")]),DS3=unlist(alpha[c("low","low")])))



########################################
### GENERATE MULTI-OMICS SIMULATIONS ###
########################################


MO = list()

for(sce in names(scenarii)){
	MO[[sce]] = list()
	for(rep in 1:n_rep){
		MO[[sce]][[rep]] = make_example_data(n_views = n_views, scenario = scenarii[[sce]])
	}
}
names(MO) = names(scenarii)
save(MO,file="Data/MOFA/final/MultiOmics_Simulations_new.RData")


###########################
### SIMULATIONS SUMMARY ###
###########################


lapply(MO, function(sce) sapply(sce, function(rep) sapply(rep$data, function(ds) dim(ds))))

lapply(MO, function(mo){
	summary(t(sapply(mo, function(rep){
		sapply(rep$S_w, function(sw){
			sum(rowSums(sw)>1)/nrow(sw)
		})
	})))
})


############################################
### CHARACTERIZE MULTI-OMICS SIMULATIONS ###
###       1) COMPUTE SIGNAL/NOISE        ###
############################################


pdf("Analyses/MOFA/final/S2N_all_scenarii.pdf",height=4,width=12)
par(mfrow=c(1,3))
for(sce in names(MO)){
	S2N = mclapply(MO[[sce]], function(x) compute_signal_to_noise_ratio(x), mc.cores=nbcores)
	S2Ncat = lapply(1:n_views, function(iview) do.call("rbind",lapply(S2N, function(x){
		if(nrow(x[[iview]])==3){t(x[[iview]])}else{NULL}
	})))
	for(iview in 1:n_views){colnames(S2Ncat[[iview]]) = c("Main_MO","Confounding_MO",paste("Data_specific",iview,sep="_"))}
	S2Ncat.df = data.frame(melt(as.data.frame(do.call("rbind",S2Ncat)),variable.name="Effect",value.name="Omega_2"), Modality = factor(rep(1:n_views,sapply(S2Ncat,length))))
	for(iview in 1:n_views){boxplot(S2Ncat[[iview]], main=c(sce,paste("Data modality",iview)),ylim=c(0,0.6), ylab="Omega^2")}
}
dev.off()


######################
### 2) COMPUTE DEG ###
######################


pdf("Analyses/MOFA/final/DEG_all_scenarii.pdf",height=4,width=12)
par(mfrow=c(1,3))
for(sce in names(MO)){
	venncounts = mclapply(MO[[sce]], function(x){
		venncounts = list()
		for(view in names(x$data)){
			fit = lmFit(x$data[[view]], design = x$Z_covar[,c("main_MO","conf_MO",paste0("DS",sub(".*([0-9]+)","\\1",view)))])
			fit = eBayes(fit = fit, trend=T)
#			plotSA(fit, main=view)
			venncounts[[view]] = vennCounts(decideTests(fit))
		}
		return(venncounts)
	}, mc.cores=nbcores)
	venncounts_avg = lapply(1:n_views, function(iview) round(Reduce("+",lapply(venncounts, function(x) x[[iview]]))/n_rep))
	for(iview in 1:n_views){print(vennDiagram(venncounts_avg[[iview]], circle.col =palSpectral(3), mar=c(0,0,2,0),cex=c(0.7,1,1), main=c(sce,paste("Data modality",iview))))}
}
dev.off()


#############################################################################
### 3) WRITE SIMULATIONS AND RUN PERMCCA USING run_permcca.sh ON THE FARM ###
#############################################################################


for(sce in names(MO)){
	for(i in 1:(length(MO[[sce]][[1]]$data)-1)){
		for(j in (i+1):length(MO[[sce]][[1]]$data)){
			for(rep in 1:n_rep){
				### WRITE DATA TABLES TO BE PROCESSED WITH PERMCCA
				indRows = lapply(c(i,j), function(k) order(rowSds(MO[[sce]][[rep]]$data[[k]]), decreasing=T)[1:min(c((2*n_samples)/3,nrow(MO[[sce]][[rep]]$data[[k]])))])
				write.table(round(t(MO[[sce]][[rep]]$data[[i]][indRows[[1]],]),digits=3),paste0("Data/MOFA/PermCCA/",sce,"_",names(MO[[sce]][[rep]]$data)[i],"_",names(MO[[sce]][[rep]]$data)[j],"_Rep",rep,"_1.txt"),quote=F,sep="\t", col.names=F, row.names=F)
				write.table(round(t(MO[[sce]][[rep]]$data[[j]][indRows[[2]],]),digits=3),paste0("Data/MOFA/PermCCA/",sce,"_",names(MO[[sce]][[rep]]$data)[i],"_",names(MO[[sce]][[rep]]$data)[j],"_Rep",rep,"_2.txt"),quote=F,sep="\t", col.names=F, row.names=F)
			}
		}
	}
}


###################################################################################################
### 4) LOAD run_permcca RESULTS, RUN CCA AND PERFORM ANOVA ON 1ST CANONICAL VARIABLE, AS CC ARE ###
###  HIGHLY CORRELATED (R>0.9) CHECK ASSOCIATION BETWEEN CANONICAL AXES AND SIMULATED EFFECTS   ###
###################################################################################################


fwer = NULL
for(sce in names(MO)){
	for(i in 1:(length(MO[[sce]][[rep]]$data)-1)){
		for(j in (i+1):length(MO[[sce]][[rep]]$data)){
			for(rep in 1:n_rep){
				fwer=rbind(fwer,data.frame(obsSim="Simulated", dataset=sce, comparison=paste(i,j,sep="_"),read.table(paste0("Data/MOFA/PermCCA/",sce,"_",names(MO[[sce]][[rep]]$data)[i],"_",names(MO[[sce]][[rep]]$data)[j],"_Rep",rep,"_fwer.txt"),sep=",")[,1:5]))
			}
		}
	}
}
fwer.df = melt(fwer, id.var=c("dataset","comparison","obsSim"), variable.name="Canonical_Variable",value.name="pvalue")
fwer.df$Canonical_Variable = sub("V","",fwer.df$Canonical_Variable)
fwer.df$pvalue = -log10(fwer.df$pvalue)

pdf("Analyses/MOFA/final/permCCA_pvalues_all_scenarii.pdf",height=2*length(MO),width=10)
print(ggplot(fwer.df, aes(x=pvalue, group=Canonical_Variable)) + geom_density(aes(fill=Canonical_Variable), alpha=0.3, bw=0.2) + facet_grid(dataset ~ comparison) + theme_bw() + xlab("-log10(p-value)"))
dev.off()


anova.df =NULL
for(sce in names(MO)){
	for(i in 1:2){
		for(j in (i+1):3){
			anova = do.call("rbind",mclapply(1:n_rep, function(rep){
				cca = rcc(t(MO[[sce]][[rep]]$data[[i]]),t(MO[[sce]][[rep]]$data[[j]]), method="shrinkage", ncomp=1)
				df = data.frame(y=cca$variates[["X"]], MO[[sce]][[rep]]$Z_covar)
				summary(aov(y ~ ., data = df))[[1]][1:5,"Pr(>F)"]
			}, mc.cores=nbcores))
			anova.df = rbind(anova.df, data.frame(obsSim="Simulated", dataset=sce, comparison=paste(names(MO[[sce]][[rep]]$data)[c(i,j)],collapse="_"), variable=rep(colnames(MO[[sce]][[1]]$Z_covar),each=nrow(anova)),pvalue=as.vector(anova)))
		}
	}
}
anova.df$pvalue = -log10(anova.df$pvalue)

pdf("Analyses/MOFA/final/CCA_anova_pvalues_all_scenarii.pdf",height=2*length(MO),width=10)
print(ggplot(anova.df, aes(x=pvalue, group=variable)) + geom_density(aes(fill=variable), alpha=0.3, bw=1) + facet_grid(dataset ~ comparison) + theme_bw() + xlab("-log10(p-value)"))
dev.off()


###########################
### 6) SAVE SIMULATIONS ###
###########################


for(sce in names(MO)){
	for(rep in 1:n_rep){
		write.table(MO[[sce]][[rep]]$Z_covar[,"main_MO"],paste0("Data/MOFA/final/",sce,"_Y_Rep_",rep,".txt"),quote=F,sep="\t", row.names=F, col.names=F)
#		for(i in 1:(length(MO[[sce]][[1]]$data))){
#			write.table(round(t(MO[[sce]][[rep]]$data[[i]]),digits=3),paste0("Data/MOFA/final/",sce,"_",names(MO[[sce]][[rep]]$data)[i],"_Rep_",rep,".txt"),quote=F,sep="\t", col.names=T, row.names=T)
#		}
	}
}


##########################################################################################
### 7) CONCAT AND COMPARE SIMULATED AND OBSERVED RESULTS FROM REAL MULTI-OMIC DATASETS ###
##########################################################################################


#if(names(MO) %in% c("rv144","MMRN","TCGA")){

#	### CONCAT FWER RESULTS
#	fwer.all = rbind(fwer.df, do.call("rbind",lapply(names(MO), function(sce){
#		load(paste0("Data/",sce,"/fwer.RData"))
#		tmp = NULL
#		for(viewname in names(fwer)){
#			tmp = rbind(tmp,data.frame(obsSim="Observed", dataset=sce, view=viewname, index=1:5, pvalue=fwer[[viewname]][1:5]))
#		}
#		return(tmp)
#	})))
#	fwer.all$pvalue = -log10(fwer.all$pvalue)
#	fwer.all$dataset = factor(fwer.all$dataset, levels=c("rv144","MMRN","TCGA"))

#	pdf("Analyses/MOFA/final/observed_vs_simulated_CanonivalCorrelation_pvalue.pdf",width=10, height=4)
#	ggplot(data=fwer.all, aes(x=index, y=pvalue, fill=view)) + geom_bar(stat="identity", position=position_dodge()) + facet_grid(obsSim~dataset) + theme_bw() + xlab("index canonical variable") + ylab("-log10(p-value)")
#	dev.off()

#	### CONCAT ANOVA RESULTS
#	anova.df = rbind(anova.df, do.call("rbind",lapply(names(MO), function(sce){
#		tmp = load(paste0("Data/",sce,"/anova_pvalue.RData"))
#		return(data.frame(obsSim="Observed", dataset=sce,get(tmp)))
#	})))
#	anova.df$dataset = factor(anova.df$dataset, levels=c("rv144","MMRN","TCGA"))

#	pdf("Analyses/MOFA/final/observed_vs_simulated_anova_pvalue.pdf",width=6, height=4)
#	ggplot(anova.df,aes(x=pvalue, group=obsSim)) + geom_density(aes(fill = obsSim), alpha = 0.4, bw=1.2) + theme_bw() + xlab("-log10(p-value)") + facet_wrap(~dataset)
#	dev.off()

#}

#### PLOT HEATMAPS
#pdf("Analyses/MOFA/Heatmap_",sce,".pdf",width=12, height=10)
#for(iview in 1:length(MO[[sce]][[rep]]$data)){
#	# Generate clinical Annotations
#	RowSideColors = ColSideColors = NULL
#	for(i in c(1:2,iview+2)){ColSideColors = cbind(ColSideColors, palSpectral(length(unique(MO[[sce]][[rep]]$Z_covar[,i])))[as.numeric(as.factor(MO[[sce]][[rep]]$Z_covar[,i]))])}
#	colnames(ColSideColors) = colnames(MO[[sce]][[rep]]$Z_covar[,c(1:2,iview+2)])
#	# Generate Variable Annotations 
#	RowSideColors = as.matrix(t(apply(MO[[sce]][[rep]]$S_w[[iview]]*MO[[sce]][[rep]]$W[[iview]],2, function(x) BWRpalette(200)[1+round(100*(x-min(x))/max(abs(x)))])))
#	rownames(RowSideColors) = colnames(MO[[sce]][[rep]]$Z_covar[,c(1:2,iview+2)])
#	# Compute random forest
#	rf.roc = compute_random_forest(t(MO[[sce]][[rep]]$data[[iview]]), as.factor(MO[[sce]][[rep]]$Z_covar$Main_MO))

#	heatmap.3(as.matrix(MO[[sce]][[rep]]$data[[iview]]), trace='none', symbreaks=F, symkey=F, col=BWRpalette(n = 128), margins=c(10,10), main=c(paste("Data modality",iview),"",paste0("AUC = ",round(auc(rf.roc),digits=2))), hclustfun=function(d){hclust(d, method = "ward.D2", members = NULL)}, RowSideColors=RowSideColors, ColSideColors = ColSideColors, RowSideColorsSize=4, ColSideColorsSize=8)
#}
#dev.off()

