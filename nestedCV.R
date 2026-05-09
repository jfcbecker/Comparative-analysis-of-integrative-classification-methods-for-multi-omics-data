getMtryGrid <- function(p) round(c(log10(p), sqrt(p), p/5, p/2))

rmSpecialChar <- function(x) gsub("\\+","plus",gsub("-","minus",gsub("[\\/ ]+","",gsub("^[0-9]+","",x))))

nestedCV <- function(X=NULL, Y=NULL, blocks=NULL, funName=NULL, nFolds=NULL, nodesizes=NULL, mtry=NULL, num.trees = 500){

	### Initialization
	iMaxMCCAllFolds = NULL
	predictions = factor(rep(NA,length(Y)),levels=c(NA,levels(Y)))
	YX_concat = data.frame(Y=Y,X, check.names = F)

	### Define grid search parameters (mtry and nodesize)
	if(funName=="randomForest"){gridParam = expand.grid(nodesize=nodesizes,mtry=mtry)}
	if(funName=="blockForest"){gridParam = do.call("rbind",lapply(nodesizes, function(nodesize) cbind(nodesize,mtry)))}

	### Create folds for Cross-validation
	folds = createFolds(Y, k = nFolds)

	### Cross validation for performance estimation 
	for(iFold in 1:nFolds){

		### Estimate block weights (median over nRepBlockWeights repetitions)
		if(funName == "blockForest"){
			blockWeights = rowMedians(sapply(1:nRepBlockWeights, function(i) blockfor(X[-folds[[iFold]],], Y[-folds[[iFold]]], blocks=blocks, num.trees = num.trees, block.method = "BlockForest")$paramvalues))
		}else{blockWeights=NULL}

		### Run RandomForest/BlockForest on hyper-parameter grid using train set
		RBForest = apply(gridParam, 1, function(param){
			if(funName == "randomForest"){
				return(randomForest(Y~.,YX_concat[-folds[[iFold]],], nodesize=param[1], mtry=param[2], ntree=num.trees))
			}
			if(funName == "blockForest"){
				return(blockForest(Y~.,YX_concat[-folds[[iFold]],], blocks=blocks, block.weights=blockWeights, min.node.size=param[1], mtry = param[-1], num.trees=num.trees, num.threads=1))
			}
		})

		### Keep best hyper-parameters based on OOB using MCC
		iMaxMCC = which.max(sapply(RBForest, function(rbf) mcc(preds=rbf[[ifelse(funName=="randomForest","predicted","predictions")]],actuals=Y[-folds[[iFold]]])))
		iMaxMCCAllFolds = c(iMaxMCCAllFolds, iMaxMCC)

		### Compute performance on test set
		pred = predict(RBForest[[iMaxMCC]], X[folds[[iFold]],])
		if(funName == "randomForest"){predictions[folds[[iFold]]] = pred}
		if(funName == "blockForest"){predictions[folds[[iFold]]] = pred$predictions}
	}
	return(list(predictions=predictions,blockWeights=blockWeights, finalGridParam=unlist(gridParam[as.numeric(names(which.max(table(iMaxMCCAllFolds)))),])))
}
