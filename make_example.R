#' @title Simulate a data set using the generative model of MOFA
#' @name make_example_data
#' @description Function to simulate an example multi-view multi-group data set according to the generative model of MOFA2.
#' @param n_views number of views
#' @param n_features number of features in each view 
#' @param n_samples number of samples in each group
#' @param n_groups number of groups
#' @param n_factors number of factors
#' @param likelihood likelihood for each view, one of "gaussian" (default), "bernoulli", "poisson",
#'  or a character vector of length n_views
#' @param lscales vector of lengthscales, needs to be of length n_factors (default is 0 - no smooth factors)
#' @param sample_cov (only for use with MEFISTO) matrix of sample covariates for one group with covariates in rows and samples in columns 
#' or "equidistant" for sequential ordering, default is NULL (no smooth factors)
#' @param as.data.frame return data and covariates as long dataframe 
#' @return Returns a list containing the simulated data and simulation parameters.
#' @importFrom stats rnorm rbinom rpois
#' @importFrom dplyr left_join
#' @importFrom stats dist
#' @export
#' @examples
#' # Generate a simulated data set
#' MOFAexample <- make_example_data()

#' @param alpha_z[["main_MO"]] precision for multi-omic latent variable
#' @param alpha_z[["conf_MO"]] precision for multi-omic confouding latent variable
#' @param alpha_z_ds precision for data specific latent variables

compute_signal_to_noise_ratio <- function(MOdata){
	lapply(1:length(MOdata$data), function(iview){
		apply(MOdata$data[[iview]],1,function(line) {
			df = data.frame(y=as.numeric(line), MOdata$Z_covar[,c(1:2,iview+2)])
			resAOV = aov(y ~ ., data = df)
			return(omega_sq(resAOV))
		})
	})
}

compute_random_forest <- function(data, status){
	rf <- foreach(ntree=rep(50, 10), .combine=randomForest::combine,
			  .multicombine=TRUE, .packages='randomForest') %dopar% {
		randomForest(as.matrix(data), status, ntree=ntree)
	}
	roc(status,rf$votes[,2])
}

make_example_data <- function(n_views, scenario, likelihood = "gaussian") {
	
	n_factors = length(scenario$alpha_z)
	n_groups = length(scenario$alpha_z[["main_MO"]])
	likelihood = rep(likelihood,n_views)
	
    # simulate multi-omics factors
    Z = do.call("cbind",lapply(scenario$alpha_z, function(alphaz){
		if(is.list(alphaz)){
			do.call("cbind",lapply(alphaz, function(alz) as.vector(unlist(sapply(seq_len(n_groups), function(gp) abs(rnorm(scenario$n_samples[gp], 0, sqrt(1/alz[gp]))))))))
		}else{
			as.vector(unlist(sapply(seq_len(n_groups), function(gp) abs(rnorm(scenario$n_samples[gp], 0, sqrt(1/alphaz[gp]))))))
		}
	}))
	# Generate permutation matrix so that main, confounding and data specific effects are not correlated
	covIdx = rep(c(1,0),time=scenario$n_samples)
	Z_covar = NULL
	for(fct in 1:ncol(Z)){
		iPerm = sample(sum(scenario$n_samples))
		Z_covar = cbind(Z_covar,covIdx[iPerm])
		Z[,fct] = Z[iPerm,fct]
	}
	colnames(Z_covar) = colnames(Z)
    
	# simulate weights
	iShared = lapply(seq_len(n_views), function(vw) sample(1:scenario$n_features[vw], round(scenario$n_features[vw]*scenario$theta_w[vw]*scenario$percentageSharedFeat)))
	S_w = lapply(seq_len(n_views), function(vw){
		Sw = matrix(0,nrow=scenario$n_features[vw], ncol=n_factors)
		Sw[iShared[[vw]],] = 1
		for(fct in 1:ncol(Sw)){
			Sw[sample((1:scenario$n_features[vw])[-iShared[[vw]]], round(scenario$n_features[vw]*scenario$theta_w[vw]*(1-scenario$percentageSharedFeat))),fct] = 1
			iShared[[vw]] = which(rowSums(Sw)>0)
		}
		return(Sw)
	})
	for(vw in which(scenario$DisableMultiomicsFromView)){
		S_w[[vw]][,1] = S_w[[vw]][,1]/8
	}

	W <- lapply(seq_len(n_views), function(vw){
		sapply(1:n_factors, function(fct){rnorm(scenario$n_features[vw], 0, sqrt(1/scenario$alpha_w[vw]))})
	})

	# pre-compute linear term and rbind groups
	mu <- lapply(seq_len(n_views), function(vw){(Z[,c(1:2,vw+2)]) %*% t(S_w[[vw]]*W[[vw]])})


	# simulate data according to the likelihood
	data <- lapply(seq_len(n_views), function(vw){
		lk <- likelihood[vw]
		if (lk == "gaussian"){
		  dd <- t(mu[[vw]] + rnorm(length(mu[[vw]]),0,sqrt(1/scenario$tau)))
		}
		else if (lk == "poisson"){
		  term <- log(1+exp(mu[[vw]]))
		  dd <- t(apply(term, 2, function(tt) rpois(length(tt),tt)))
		}
		else if (lk == "bernoulli") {
		  term <- 1/(1+exp(-mu[[vw]]))
		  dd <- t(apply(term, 2, function(tt) rbinom(length(tt),1,tt)))
		}
		colnames(dd) <- paste0("sample_", seq_len(ncol(dd)))
		rownames(dd) <- paste0("feature_", seq_len(nrow(dd)),"_view", vw)
		dd
	})

	names(data) <- paste0("view_", seq_len(n_views))
	for(i in seq_len(n_views)){
		rownames(data[[i]]) = paste0("Feature_",1:nrow(data[[i]]))
		colnames(data[[i]]) = paste0("Sample_",1:ncol(data[[i]]))
	}

	return(list(data = data, Z = Z, S_w=S_w, W=W, Z_covar=Z_covar))
}
