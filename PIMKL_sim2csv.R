library(matrixTests)

setwd("/sps/bioaster/Projects/ISYBIO/")
load("Data_Simulations/Data/MOFA/final/MultiOmics_Simulations_final.RData")
n_rep=40
dir.create("Benchmark/Data/Predictive/PIMKL/simulations/preprocess",recursive=T)
dir.create("Benchmark/Data/Predictive/PIMKL/simulations/output")


for(sim in names(MO)){
	cat(sim, "\n")
	for(rep in 1:n_rep){
		interaction<-c()
		featuresetF<-c()
		datatowriteF<-c()
		try(for(i in 1:3){
			datatowrite<-t(scale(t(MO[[sim]][[rep]]$data[[i]]),T,T))
			colnames(datatowrite)<-gsub("_","",colnames(datatowrite))
			rownames(datatowrite)<-paste0("view",i,gsub("_","",rownames(datatowrite)))
			datatowriteF<-rbind(datatowriteF,datatowrite)
			write.table(datatowrite,paste0("Benchmark/Data/Predictive/PIMKL/simulations/",gsub("_","",sim),"",gsub("_","",names(MO[[sim]][[rep]]$data)[i]),"Rep",rep,".csv"),quote=F,sep=",", col.names=T, row.names=T)
			interactFeat<-as.data.frame(t(combn(rownames(datatowrite),2)))
			colnames(interactFeat)<-c("e1","e2")
			interactFeat$intensity<-1
			interaction<-rbind(interaction,interactFeat)
			featuresetF <- rbind(featuresetF,data.frame(V1=gsub("_","",names(MO[[sim]][[rep]]$data)[i]),V2=gsub("_","",names(MO[[sim]][[rep]]$data)[i]),V3=paste(rownames(datatowrite),collapse="\t")))
		})

		write.table(interaction,paste0("Benchmark/Data/Predictive/PIMKL/simulations/interaction",gsub("_","",sim),"Rep",rep,".csv"),quote=F,sep=",", col.names=F, row.names=F)
		
		write.table(featuresetF,paste0("Benchmark/Data/Predictive/PIMKL/simulations/featureset",gsub("_","",sim),"Rep",rep,".gmt"),quote=F,sep="\t", col.names=F, row.names=F)

		labels<-cbind(colnames(datatowrite),MO[[sim]][[rep]]$Z_covar[,1])
		write.table(labels,paste0("Benchmark/Data/Predictive/PIMKL/simulations/labels",gsub("_","",sim),"Rep",rep,".csv"),quote=F,sep=",", col.names=F, row.names=F)

	}
}
