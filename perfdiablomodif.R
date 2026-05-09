perf.sgccdamodif<-
function (object, dist = c("all", "max.dist", "centroids.dist", 
    "mahalanobis.dist"), validation = c("Mfold", "loo"), folds = 10, 
    nrepeat = 1, auc = FALSE, progressBar = FALSE, signif.threshold = 0.01, 
    cpus = 1, ...) 
{
    X = object$X
    level.Y = object$names$colnames$Y
    J = length(X)
    Y = object$Y
    n = nrow(X[[1]])
    indY = object$indY
    max.iter = object$max.iter
    ncomp <- max(object$ncomp[-indY])
    near.zero.var = !is.null(object$nzv)
    if (any(dist == "all")) {
        dist.select = c("max.dist", "centroids.dist", "mahalanobis.dist")
    }else {
        dist.select = dist
    }
    dist = match.arg(dist.select, choices = c("all", "max.dist", 
        "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
    signif.threshold <- .check_alpha(signif.threshold)
    if (length(validation) > 1) 
        validation = validation[1]
    if (!(validation %in% c("Mfold", "loo"))) 
        stop("Choose 'validation' among the two following possibilities: 'Mfold' or 'loo'")
    keepX = object$keepX
    error.mat = error.mat.class = Y.all = predict.all = Y.predict = list.features = final.features = weights = crit = list()
    if (length(X) > 1) 
        Y.mean = Y.mean.res = Y.weighted.vote = Y.weighted.vote.res = Y.vote = Y.vote.res = Y.WeightedPredict = Y.WeightedPredict.res = list()
    if (progressBar == TRUE) {
        cat(sprintf("\nPerforming repeated cross-validation with nrepeat = %s...\n", 
            nrepeat))
        pb = txtProgressBar(style = 3)
    }
    repeat_cv_perf.diablo <- function(nrep) {
        if (progressBar == TRUE) {
            setTxtProgressBar(pb = pb, value = nrep/nrepeat)
        }
        if (validation == "Mfold") {
            if (is.null(folds) || !is.numeric(folds) || folds < 
                2 || folds > n) {
                stop("Invalid number of folds.")
            }else {
                M = round(folds)
                temp = stratified.subsampling(Y, folds = M)
                folds = temp$SAMPLE
                if (temp$stop > 0) 
                  warning("At least one class is not represented in one fold, which may unbalance the error rate.\n  Consider a number of folds lower than the minimum in table(Y): ", 
                    min(table(Y)))
            }
        }else if (validation == "loo") {
            folds = split(1:n, rep(1:n, length = n))
            M = n
        }else {
            stop("validation can be only 'Mfold' or 'loo'")
        }
        M = length(folds)
        X.training = lapply(folds, function(x) {
            out = lapply(1:J, function(y) {
                X[[y]][-x, ]
            })
            names(out) = names(X)
            out
        })
        Y.training = lapply(folds, function(x) {
            Y[-x]
        })
        X.test = lapply(folds, function(x) {
            out = lapply(1:J, function(y) {
                X[[y]][x, , drop = FALSE]
            })
            names(out) = names(X)
            out
        })
        Y.test = lapply(folds, function(x) {
            Y[x]
        })
        block.splsda.args <- function(ind) {
            return(list(X = X.training[[ind]], Y = Y.training[[ind]], 
                ncomp = ncomp, keepX = keepX, design = object$design, 
                max.iter = max.iter, tol = object$tol, init = object$init, 
                scheme = object$scheme, near.zero.var = near.zero.var))
        }
        model = lapply(1:M, function(x) {
            suppressWarnings(do.call("block.splsda", block.splsda.args(x)))
        })
        if (auc) {
            auc.rep.block.comp.fold <- lapply(1:M, function(x) {
                auroc.sgccda(object = model[[x]], plot = FALSE, 
                  print = FALSE)
            })
            blocks <- .name_list(names(auc.rep.block.comp.fold[[1]]))
            comps <- .name_list(names(auc.rep.block.comp.fold[[1]][[1]]))
            auc.rep.block.comp <- lapply(blocks, function(block) {
                lapply(comps, function(comp) {
                  Reduce("+", lapply(auc.rep.block.comp.fold, 
                    function(fold) {
                      fold[[block]][[comp]]
                    }))/M
                })
            })
            auc.rep.comp <- lapply(comps, function(comp) {
                Reduce("+", lapply(blocks, function(block) {
                  auc.rep.block.comp[[block]][[comp]]
                }))/length(blocks)
            })
        }
        crit = lapply(1:M, function(x) {
            model[[x]]$crit
        })
        weights = lapply(1:M, function(fold_i) {
            fold_weights <- model[[fold_i]]$weights
            fold_weights$fold <- fold_i
            fold_weights$rep <- nrep
            fold_weights$block <- rownames(fold_weights)
            rownames(fold_weights) <- NULL
            fold_weights
        })
        weights <- Reduce(rbind, weights)
        features = lapply(1:J, function(x) {
            lapply(1:object$ncomp[x], function(y) {
                unlist(lapply(1:M, function(z) {
                  if (is.null(colnames(X[[x]]))) {
                    paste0("X", which(model[[z]]$loadings[[x]][, 
                      y] != 0))
                  }else {
                    if (!is.null(colnames(X[[x]]))) 
                      colnames(X[[x]])[model[[z]]$loadings[[x]][, 
                        y] != 0]
                  }
                }))
            })
        })
        list.features = lapply(1:J, function(x) {
            lapply(features[[x]], function(y) {
                (sort(table(factor(y))/M, decreasing = TRUE))
            })
        })
        final.features = lapply(1:J, function(x) {
            lapply(1:object$ncomp[x], function(y) {
                temp = as.data.frame(object$loadings[[x]][which(object$loadings[[x]][, 
                  y] != 0), y, drop = FALSE])
                if (is.null(colnames(X[[x]]))) {
                  row.names(temp) = paste0("X", which(object$loadings[[x]][, 
                    y] != 0))
                }else {
                  if (!is.null(colnames(X[[x]]))) 
                    row.names(temp) = colnames(X[[x]])[which(object$loadings[[x]][, 
                      y] != 0)]
                }
                names(temp) = "value.var"
                return(temp[sort(abs(temp[, 1]), index.return = TRUE, 
                  decreasing = TRUE)$ix, 1, drop = FALSE])
            })
        })
        list.features = lapply(1:J, function(x) {
            names(list.features[[x]]) = paste0("comp", 1:object$ncomp[x])
            return(list.features[[x]])
        })
        final.features = lapply(1:J, function(x) {
            names(final.features[[x]]) = paste0("comp", 1:object$ncomp[x])
            return(final.features[[x]])
        })
        predict.all = lapply(1:M, function(x) {
            predict.block.spls(model[[x]], X.test[[x]], dist = "all")
        })
        Y.predict = lapply(1:M, function(x) {
            predict.all[[x]]$class
        })
        Y.all = lapply(1:M, function(x) {
            predict.all[[x]]$predict
        })
        Y.all = lapply(1:J, function(x) {
            lapply(1:object$ncomp[x], function(y) {
                lapply(1:M, function(z) {
                  Y.all[[z]][[x]][, , y]
                })
            })
        })
        Y.all = lapply(1:J, function(x) {
            lapply(1:object$ncomp[x], function(y) {
                do.call(rbind, Y.all[[x]][[y]])
            })
        })
        Y.all = lapply(1:J, function(x) {
            lapply(1:object$ncomp[x], function(y) {
                row.names(Y.all[[x]][[y]]) = row.names(X[[x]])[unlist(folds)]
                return(Y.all[[x]][[y]])
            })
        })
        Y.all = lapply(1:J, function(x) {
            lapply(1:object$ncomp[x], function(y) {
                colnames(Y.all[[x]][[y]]) = levels(Y)
                return(Y.all[[x]][[y]])
            })
            names(Y.all[[x]]) = paste0("comp", 1:object$ncomp[x])
            return(Y.all[[x]])
        })
        Y.predict = lapply(1:J, function(x) {
            lapply(1:M, function(y) {
                lapply(which(c("max.dist", "centroids.dist", 
                  "mahalanobis.dist") %in% dist.select), function(z) {
                  Y.predict[[y]][[z]][[x]]
                })
            })
        })
        Y.predict = lapply(1:J, function(x) {
            lapply(1:M, function(y) {
                lapply(1:length(dist.select), function(z) {
                  if (!is.null(row.names(X[[x]])[folds[[y]]])) {
                    row.names(Y.predict[[x]][[y]][[z]]) = row.names(X[[x]])[folds[[y]]]
                  }
                  else {
                    row.names(Y.predict[[x]][[y]][[z]]) = paste0("Ind", 
                      unlist(folds))
                  }
                  return(Y.predict[[x]][[y]][[z]])
                })
            })
        })
        Y.predict = lapply(1:J, function(x) {
            lapply(1:length(dist.select), function(y) {
                lapply(1:M, function(z) {
                  Y.predict[[x]][[z]][[y]]
                })
            })
        })
        Y.predict = lapply(1:J, function(x) {
            lapply(1:length(dist.select), function(y) {
                do.call(rbind, Y.predict[[x]][[y]])
            })
        })
        Y.predict = lapply(1:J, function(x) {
            names(Y.predict[[x]]) = dist.select
            return(Y.predict[[x]])
        })
        error.mat = lapply(1:J, function(x) {
            lapply(dist.select, function(y) {
                apply(Y.predict[[x]][[y]], 2, function(z) {
                  1 - sum(diag(table(factor(z, levels = levels(Y)), 
                    Y[unlist(folds)])))/length(Y)
                })
            })
        })
        error.mat = lapply(1:J, function(x) {
            do.call(cbind, error.mat[[x]])
        })
        error.mat = lapply(1:J, function(x) {
            colnames(error.mat[[x]]) = dist.select
            return(error.mat[[x]])
        })
        error.mat.class = lapply(1:J, function(x) {
            lapply(dist.select, function(y) {
                apply(Y.predict[[x]][[y]], 2, function(z) {
                  temp = diag(table(factor(z, levels = levels(Y)), 
                    Y[unlist(folds)]))
                  1 - c(temp/summary(Y))
                })
            })
        })
        error.mat.class = lapply(1:J, function(x) {
            names(error.mat.class[[x]]) = dist.select
            return(error.mat.class[[x]])
        })
        if (!is.null(names(X))) {
            names(error.mat) = names(X)
            names(error.mat.class) = names(X)
            names(list.features) = names(X)
            names(final.features) = names(X)
            names(Y.all) = names(X)
            names(Y.predict) = names(X)
        }
        if (length(X) > 1) {
            Y.mean = lapply(1:M, function(y) {
                predict.all[[y]][["AveragedPredict.class"]][[1]]
            })
            Y.mean = do.call(rbind, Y.mean)
            Y.mean = Y.mean[sort(unlist(folds), index.return = TRUE)$ix, 
                , drop = FALSE]
            Y.mean.res = sapply(1:ncomp, function(x) {
                mat = table(factor(Y.mean[, x], levels = levels(Y)), 
                  Y)
                mat2 <- mat
                diag(mat2) <- 0
                err = c(c(colSums(mat2)/summary(Y), sum(mat2)/length(Y)), 
                  mean(colSums(mat2)/colSums(mat)))
            })
            colnames(Y.mean.res) = paste0("comp", 1:ncomp)
            row.names(Y.mean.res) = c(levels(Y), "Overall.ER", 
                "Overall.BER")
            Y.WeightedPredict = lapply(1:M, function(y) {
                predict.all[[y]][["WeightedPredict.class"]][[1]]
            })
            Y.WeightedPredict = do.call(rbind, Y.WeightedPredict)
            Y.WeightedPredict = Y.WeightedPredict[sort(unlist(folds), 
                index.return = TRUE)$ix, , drop = FALSE]
            Y.WeightedPredict.res = sapply(1:ncomp, function(x) {
                mat = table(factor(Y.WeightedPredict[, x], levels = levels(Y)), 
                  Y)
                mat2 <- mat
                diag(mat2) <- 0
                err = c(c(colSums(mat2)/summary(Y), sum(mat2)/length(Y)), 
                  mean(colSums(mat2)/colSums(mat)))
            })
            colnames(Y.WeightedPredict.res) = paste0("comp", 
                1:ncomp)
            row.names(Y.WeightedPredict.res) = c(levels(Y), "Overall.ER", 
                "Overall.BER")
            Y.weighted.vote = lapply(dist.select, function(x) {
                lapply(1:M, function(y) {
                  predict.all[[y]][["WeightedVote"]][[x]]
                })
            })
            Y.weighted.vote = lapply(1:length(dist.select), function(x) {
                do.call(rbind, Y.weighted.vote[[x]])
            })
            Y.weighted.vote = lapply(1:length(dist.select), function(x) {
                Y.weighted.vote[[x]][sort(unlist(folds), index.return = TRUE)$ix, 
                  , drop = FALSE]
            })
            names(Y.weighted.vote) = dist.select
            Y.weighted.vote.res = lapply(1:length(dist.select), 
                function(x) {
                  apply(Y.weighted.vote[[x]], 2, function(y) {
                    y[is.na(y)] <- nlevels(Y) + 5
                    temp = table(factor(y, levels = c(levels(Y), 
                      nlevels(Y) + 5)), Y)
                    diag(temp) <- 0
                    err = c(colSums(temp)/summary(Y), sum(temp)/length(Y), 
                      mean(colSums(temp)/summary(Y)))
                    return(err = err)
                  })
                })
            Y.weighted.vote.res = lapply(1:length(dist.select), 
                function(x) {
                  colnames(Y.weighted.vote.res[[x]]) = paste0("comp", 
                    1:max(object$ncomp[-(J + 1)]))
                  row.names(Y.weighted.vote.res[[x]]) = c(levels(Y), 
                    "Overall.ER", "Overall.BER")
                  return((Y.weighted.vote.res[[x]]))
                })
            names(Y.weighted.vote) = dist.select
            names(Y.weighted.vote.res) = dist.select
            Y.vote = lapply(dist.select, function(x) {
                lapply(1:M, function(y) {
                  predict.all[[y]][["MajorityVote"]][[x]]
                })
            })
            Y.vote = lapply(1:length(dist.select), function(x) {
                do.call(rbind, Y.vote[[x]])
            })
            Y.vote = lapply(1:length(dist.select), function(x) {
                Y.vote[[x]][sort(unlist(folds), index.return = TRUE)$ix, 
                  , drop = FALSE]
            })
            names(Y.vote) = dist.select
            Y.vote.res = lapply(1:length(dist.select), function(x) {
                apply(Y.vote[[x]], 2, function(y) {
                  y[is.na(y)] <- nlevels(Y) + 5
                  temp = table(factor(y, levels = c(levels(Y), 
                    nlevels(Y) + 5)), Y)
                  diag(temp) <- 0
                  err = c(colSums(temp)/summary(Y), sum(temp)/length(Y), 
                    mean(colSums(temp)/summary(Y)))
                  return(err = err)
                })
            })
            Y.vote.res = lapply(1:length(dist.select), function(x) {
                colnames(Y.vote.res[[x]]) = paste0("comp", 1:max(object$ncomp[-(J + 
                  1)]))
                row.names(Y.vote.res[[x]]) = c(levels(Y), "Overall.ER", 
                  "Overall.BER")
                return((Y.vote.res[[x]]))
            })
            names(Y.vote) = dist.select
            names(Y.vote.res) = dist.select
            repeat_cv_res <- list(error.mat = error.mat, error.mat.class = error.mat.class, 
                Y.mean = Y.mean, Y.mean.res = Y.mean.res, Y.WeightedPredict.res = Y.WeightedPredict.res, Y.WeightedPredict=Y.WeightedPredict,
                Y.weighted.vote.res = Y.weighted.vote.res, Y.vote.res = Y.vote.res, Y.vote = Y.vote, Y.weighted.vote=Y.weighted.vote,
                Y.all = Y.all, predict.all = predict.all, Y.predict = Y.predict, 
                list.features = list.features, final.features = final.features, 
                crit = crit, weights = weights)
            if (auc) {
                repeat_cv_res$auc.rep.comp <- auc.rep.comp
            }
            return(repeat_cv_res)
        }
    }
    nrep_list <- as.list(seq_len(nrepeat))
    names(nrep_list) <- paste0("nrep", nrep_list)
    cpus <- .check_cpus(cpus)
    parallel <- cpus > 1
    if (parallel) {
        if (progressBar == TRUE) {
            message("progressBar unavailable in parallel computing mode for perf.block.splsda ...")
        }
        if (.onUnix()) {
            cl <- makeForkCluster(cpus)
        }else {
            cl <- makePSOCKcluster(cpus)
        }
        on.exit(stopCluster(cl))
        clusterEvalQ(cl, library(mixOmics))
        clusterExport(cl, ls(), environment())
        if (!is.null(list(...)$seed)) {
            RNGversion(.mixo_rng())
            clusterSetRNGStream(cl = cl, iseed = list(...)$seed)
        }
        repeat_cv_perf.diablo_res <- parLapply(cl, nrep_list, 
            function(nrep) repeat_cv_perf.diablo(nrep))
    }else {
        repeat_cv_perf.diablo_res <- lapply(nrep_list, function(nrep) {
            repeat_cv_perf.diablo(nrep)
        })
    }
    repeat_cv_perf.diablo_res <- .unlist_repeat_cv_output(repeat_cv_perf.diablo_res)
    auc.rep.comp <- NULL
    list2env(repeat_cv_perf.diablo_res, envir = environment())
    #rm(repeat_cv_perf.diablo_res)
    if (nrepeat > 1) {
        error.rate.all = error.mat
        error.rate = list()
        error.rate.sd = list()
        error.rate.per.class.all = error.mat.class
        error.rate.per.class = relist(0, skeleton = error.mat.class[[1]])
        error.rate.per.class.sd = relist(0, skeleton = error.mat.class[[1]])
        temp.error.rate = array(0, c(dim(error.rate.all[[1]][[1]]), 
            nrepeat))
        temp.error.rate.per.class = array(0, c(dim(error.rate.per.class.all[[1]][[1]][[1]]), 
            nrepeat))
        for (i in 1:J) {
            for (nrep in 1:nrepeat) temp.error.rate[, , nrep] = error.rate.all[[nrep]][[i]]
            temp.error.rate.mean = apply(temp.error.rate, c(1, 
                2), mean)
            temp.error.rate.sd = apply(temp.error.rate, c(1, 
                2), sd)
            dimnames(temp.error.rate.mean) = dimnames(temp.error.rate.sd) = dimnames(error.rate.all[[nrep]][[i]])
            error.rate[[i]] = temp.error.rate.mean
            error.rate.sd[[i]] = temp.error.rate.sd
            for (j in 1:length(dist.select)) {
                for (nrep in 1:nrepeat) temp.error.rate.per.class[, 
                  , nrep] = error.rate.per.class.all[[nrep]][[i]][[j]]
                temp.error.rate.per.class.mean = apply(temp.error.rate.per.class, 
                  c(1, 2), mean)
                temp.error.rate.per.class.sd = apply(temp.error.rate.per.class, 
                  c(1, 2), sd)
                dimnames(temp.error.rate.per.class.mean) = dimnames(temp.error.rate.per.class.sd) = dimnames(error.rate.per.class.all[[nrep]][[i]][[j]])
                error.rate.per.class[[i]][[j]] = temp.error.rate.per.class.mean
                error.rate.per.class.sd[[i]][[j]] = temp.error.rate.per.class.sd
            }
        }
        names(error.rate) = names(error.rate.sd) = names(X)
        AveragedPredict.error.rate.all = Y.mean.res
        temp.AveragedPredict.error.rate = array(unlist(AveragedPredict.error.rate.all), 
            c(dim(AveragedPredict.error.rate.all[[1]]), nrepeat))
        AveragedPredict.error.rate = apply(temp.AveragedPredict.error.rate, 
            c(1, 2), mean)
        AveragedPredict.error.rate.sd = apply(temp.AveragedPredict.error.rate, 
            c(1, 2), sd)
        dimnames(AveragedPredict.error.rate) = dimnames(AveragedPredict.error.rate.sd) = dimnames(AveragedPredict.error.rate.all[[1]])
        WeightedPredict.error.rate.all = Y.WeightedPredict.res
        temp.WeightedPredict.error.rate = array(unlist(WeightedPredict.error.rate.all), 
            c(dim(WeightedPredict.error.rate.all[[1]]), nrepeat))
        WeightedPredict.error.rate = apply(temp.WeightedPredict.error.rate, 
            c(1, 2), mean)
        WeightedPredict.error.rate.sd = apply(temp.WeightedPredict.error.rate, 
            c(1, 2), sd)
        dimnames(WeightedPredict.error.rate) = dimnames(WeightedPredict.error.rate.sd) = dimnames(WeightedPredict.error.rate.all[[1]])
        MajorityVote.error.rate.all = Y.vote.res
        MajorityVote.error.rate = relist(0, skeleton = MajorityVote.error.rate.all[[1]])
        MajorityVote.error.rate.sd = relist(0, skeleton = MajorityVote.error.rate.all[[1]])
        WeightedVote.error.rate.all = Y.weighted.vote.res
        WeightedVote.error.rate = relist(0, skeleton = WeightedVote.error.rate.all[[1]])
        WeightedVote.error.rate.sd = relist(0, skeleton = WeightedVote.error.rate.all[[1]])
        temp.MajorityVote.error.rate = array(0, c(dim(MajorityVote.error.rate.all[[1]][[1]]), 
            nrepeat))
        temp.WeightedVote.error.rate = array(0, c(dim(WeightedVote.error.rate.all[[1]][[1]]), 
            nrepeat))
        for (j in 1:length(dist.select)) {
            for (nrep in 1:nrepeat) {
                temp.MajorityVote.error.rate[, , nrep] = MajorityVote.error.rate.all[[nrep]][[j]]
                temp.WeightedVote.error.rate[, , nrep] = WeightedVote.error.rate.all[[nrep]][[j]]
            }
            temp.MajorityVote.error.rate.mean = apply(temp.MajorityVote.error.rate, 
                c(1, 2), mean)
            temp.MajorityVote.error.rate.sd = apply(temp.MajorityVote.error.rate, 
                c(1, 2), sd)
            dimnames(temp.MajorityVote.error.rate.mean) = dimnames(temp.MajorityVote.error.rate.sd) = dimnames(MajorityVote.error.rate.all[[nrep]][[j]])
            MajorityVote.error.rate[[j]] = temp.MajorityVote.error.rate.mean
            MajorityVote.error.rate.sd[[j]] = temp.MajorityVote.error.rate.sd
            temp.WeightedVote.error.rate.mean = apply(temp.WeightedVote.error.rate, 
                c(1, 2), mean)
            temp.WeightedVote.error.rate.sd = apply(temp.WeightedVote.error.rate, 
                c(1, 2), sd)
            dimnames(temp.WeightedVote.error.rate.mean) = dimnames(temp.WeightedVote.error.rate.sd) = dimnames(WeightedVote.error.rate.all[[nrep]][[j]])
            WeightedVote.error.rate[[j]] = temp.WeightedVote.error.rate.mean
            WeightedVote.error.rate.sd[[j]] = temp.WeightedVote.error.rate.sd
        }
    }
    else {
        error.rate = error.mat[[1]]
        error.rate.per.class = error.mat.class[[1]]
        AveragedPredict.error.rate = Y.mean.res[[1]]
        WeightedPredict.error.rate = Y.WeightedPredict.res[[1]]
        MajorityVote.error.rate = Y.vote.res[[1]]
        WeightedVote.error.rate = Y.weighted.vote.res[[1]]
    }
    result = list()
    result$error.rate = error.rate
    if (nrepeat > 1) {
        result$error.rate.sd = error.rate.sd
        result$error.rate.all = error.rate.all
    }
    result$error.rate.per.class = error.rate.per.class
    if (nrepeat > 1) {
        result$error.rate.per.class.sd = error.rate.per.class.sd
        result$error.rate.per.class.all = error.rate.per.class.all
    }
    result$predict = Y.all
    result$class = Y.predict
    result$features$stable = list.features
    if (length(X) > 1) {
        result$AveragedPredict.class = Y.mean
        result$AveragedPredict.error.rate = AveragedPredict.error.rate
        if (nrepeat > 1) {
            result$AveragedPredict.error.rate.sd = AveragedPredict.error.rate.sd
            result$AveragedPredict.error.rate.all = AveragedPredict.error.rate.all
        }
        result$WeightedPredict.class = Y.WeightedPredict
        result$WeightedPredict.error.rate = WeightedPredict.error.rate
        if (nrepeat > 1) {
            result$WeightedPredict.error.rate.sd = WeightedPredict.error.rate.sd
            result$WeightedPredict.error.rate.all = WeightedPredict.error.rate.all
        }
        result$MajorityVote = Y.vote
        result$MajorityVote.error.rate = MajorityVote.error.rate
        if (nrepeat > 1) {
            result$MajorityVote.error.rate.sd = MajorityVote.error.rate.sd
            result$MajorityVote.error.rate.all = MajorityVote.error.rate.all
        }
        result$WeightedVote = Y.weighted.vote
        result$WeightedVote.error.rate = WeightedVote.error.rate
        if (nrepeat > 1) {
            result$WeightedVote.error.rate.sd = WeightedVote.error.rate.sd
            result$WeightedVote.error.rate.all = WeightedVote.error.rate.all
        }
        weights <- Reduce(rbind, weights)
        weights <- weights[order(weights$block), ]
        result$weights <- weights
    }
    measure = c("Overall.ER", "Overall.BER")
    ncomp_opt = vector("list", length = 4)
    names(ncomp_opt) = c("AveragedPredict", "WeightedPredict", 
        "MajorityVote", "WeightedVote")
    if (nrepeat > 2 & min(object$ncomp) > 1) {
        for (prediction_framework in names(ncomp_opt)) {
            if (prediction_framework %in% c("AveragedPredict", 
                "WeightedPredict")) {
                ncomp_opt[[prediction_framework]] = matrix(NA, 
                  nrow = length(measure), ncol = 1, dimnames = list(measure))
                for (measure_i in measure) {
                  mat.error.rate = sapply(get(paste0(prediction_framework, 
                    ".error.rate.all")), function(x) {
                    x[measure_i, ]
                  })
                  ncomp_opt[[prediction_framework]][measure_i, 
                    ] = t.test.process(t(mat.error.rate), alpha = signif.threshold)
                }
            }
            else {
                ncomp_opt[[prediction_framework]] = matrix(NA, 
                  nrow = length(measure), ncol = length(dist.select), 
                  dimnames = list(measure, dist.select))
                for (measure_i in measure) {
                  for (ijk in dist.select) {
                    mat.error.rate = sapply(get(paste0(prediction_framework, 
                      ".error.rate.all")), function(x) {
                      x[[ijk]][measure_i, ]
                    })
                    ncomp_opt[[prediction_framework]][measure_i, 
                      ijk] = t.test.process(t(mat.error.rate), 
                      alpha = signif.threshold)
                  }
                }
            }
        }
    }
    result$meth = "sgccda.mthd"
    class(result) = "perf.sgccda.mthd"
    result$call = match.call()
    result$crit = crit
    result$choice.ncomp = ncomp_opt
    result$predictall<-predict.all
    result$repcvres<-repeat_cv_perf.diablo_res
    if (auc) {
        auc.comp <- lapply(.name_list(names(auc.rep.comp[[1]])), 
            function(comp) {
                Reduce("+", lapply(seq_len(nrepeat), function(rep) {
                  auc.rep.comp[[rep]][[comp]]
                }))/nrepeat
            })
        result$auc <- auc.comp
    }
    return(invisible(result))
}

stratified.subsampling <- function(Y, folds = 10)
{
    
    stop <- 0
    for (i in seq_len(nlevels(Y)))
    {
        ai <- sample(which(Y == levels(Y)[i]), replace = FALSE) # random sampling of the samples from level i
        aai <- suppressWarnings(split(ai, factor(seq_len(
            min(folds, length(ai))
        ))))                       # split of the samples in k-folds
        if (length(ai) < folds)
            # if one level doesn't have at least k samples, the list is completed with "integer(0)"
        {
            for (j in (length(ai) + 1):folds)
                aai[[j]] <- integer(0)
            stop <- stop + 1
        }
        assign(paste("aa", i, sep = "_"), sample(aai, replace = FALSE))         # the `sample(aai)' is to avoid the first group to have a lot more data than the rest
    }
      # combination of the different split aa_i into SAMPLE
    SAMPLE <- list()
    for (j in seq_len(folds))
    {
        SAMPLE[[j]] <- integer(0)
        for (i in seq_len(nlevels(Y)))
        {
            SAMPLE[[j]] <- c(SAMPLE[[j]], get(paste("aa", i, sep = "_"))[[j]])
        }
    }# SAMPLE is a list of k splits
    
    ind0 <- sapply(SAMPLE, length)
    if (any(ind0 == 0))
    {
        SAMPLE <- SAMPLE [-which(ind0 == 0)]
        message(
            "Because of a too high number of 'folds' required, ",
            length(which(ind0 == 0)),
            " folds were randomly assigned no data: the number of 'folds' is reduced to ",
            length(SAMPLE)
        )
    }
    
    return(list(SAMPLE = SAMPLE, stop = stop))
}
  
.onUnix <- function() {
    return(ifelse(.Platform$OS.type == "unix", TRUE, FALSE))
}

.check_alpha <- function(alpha=NULL) {
    if (is.null(alpha))
        alpha <- 0.01
    
    if (!is.numeric(alpha) || alpha < 0 || alpha > 1)
        stop("invalid 'signif.threshold'. Use 0.01 or 0.05.", call. = FALSE)
    alpha
}

.check_cpus <- function(cpus) {
    if (!is.numeric(cpus) || cpus <= 0)
        stop("Number of CPUs should be a positive integer.", call. = FALSE)
    
    if (cpus < 0 || is.null(cpus))
        cpus <- 1
    
    if (cpus > parallel::detectCores())
        message(sprintf("\nOnly %s CPUs available for parallel processing.\n", parallel::detectCores()))
    return(cpus)
    
}

predict.block.spls<-
function (object, newdata, study.test, dist = c("all", "max.dist", 
    "centroids.dist", "mahalanobis.dist"), multilevel = NULL, 
    ...) 
{
    time = FALSE
    newdata.scale = misdata.all = is.na.X = is.na.newdata = noAveragePredict = NULL
    if (is(object, c("rgcca", "sparse.rgcca"))) 
        stop("no prediction for RGCCA methods")
    if (!any(dist %in% c("all", "max.dist", "centroids.dist", 
        "mahalanobis.dist")) & is(object, "DA")) 
        stop("ERROR : choose one of the four following modes: 'all', 'max.dist', 'centroids.dist' or 'mahalanobis.dist'")
    ncomp = object$ncomp
    if (is(object, "DA")) {
        Y.factor = object$Y
        Y = object$ind.mat
    }
    else {
        Y = object$Y
        if (is.null(Y)) {
            Y = object$X[[object$indY]]
        }
    }
    q = ncol(Y)
    if (time) 
        time1 = proc.time()
    if (hasArg(newdata.scale)) {
        newdata = list(...)$newdata.scale
        object$logratio = NULL
        multilevel = NULL
    }
    mint.object = c("mint.pls", "mint.spls", "mint.plsda", "mint.splsda")
    block.object = c("block.pls", "block.spls", "block.plsda", 
        "block.splsda")
    if (!is(object, block.object)) {
        p = ncol(object$X)
        if (is.list(object$X)) 
            stop("Something is wrong, object$X should be a matrix and it appears to be a list")
        if (is.list(newdata) & !is.data.frame(newdata)) 
            stop("'newdata' must be a numeric matrix")
        if (length(object$nzv$Position) > 0) 
            newdata = newdata[, -object$nzv$Position, drop = FALSE]
        if (all.equal(colnames(newdata), colnames(object$X)) != 
            TRUE) 
            stop("'newdata' must include all the variables of 'object$X'")
        if (length(dim(newdata)) == 2) {
            if (ncol(newdata) != p) 
                stop("'newdata' must be a numeric matrix with ncol = ", 
                  p, " or a vector of length = ", p, ".")
        }
        if (length(dim(newdata)) == 0) {
            if (length(newdata) != p) 
                stop("'newdata' must be a numeric matrix with ncol = ", 
                  p, " or a vector of length = ", p, ".")
            dim(newdata) = c(1, p)
        }
        check = Check.entry.single(newdata, ncomp, q = 1)
        newdata = check$X
        if (length(rownames(newdata)) == 0) 
            rownames(newdata) = 1:nrow(newdata)
        if (max(table(rownames(newdata))) > 1) 
            stop("samples should have a unique identifier/rowname")
        X = list(X = object$X)
        object$X = X
        newdata = list(newdata = newdata)
        object$indY = 2
        ind.match = 1
    }
    else {
        if (!is.list(newdata)) 
            stop("'newdata' should be a list")
        if (!is.null(object$indY) && !is(object, "DA")) {
            X = object$X[-object$indY]
        }
        else {
            X = object$X
        }
        object$X = X
        p = lapply(X, ncol)
        if (length(unique(names(newdata))) != length(newdata) | 
            sum(is.na(match(names(newdata), names(X)))) > 0) 
            stop("Each entry of 'newdata' must have a unique name corresponding to a block of 'X'")
        if (length(unique(sapply(newdata, nrow))) != 1) 
            stop("All entries of 'newdata' must have the same number of rows")
        ind.match = match(names(X), names(newdata))
        if (any(is.na(ind.match))) 
            warning("Some blocks are missing in 'newdata'; the prediction is based on the following blocks only: ", 
                paste(names(X)[!is.na(ind.match)], collapse = ", "))
        newdataA = list()
        for (q in 1:length(X)) {
            if (!is.na(ind.match[q])) {
                newdataA[[q]] = newdata[[ind.match[q]]]
            }
            else {
                newdataA[[q]] = matrix(0, nrow = unique(sapply(newdata, 
                  nrow)), ncol = ncol(X[[q]]), dimnames = list(rownames(newdata[[1]]), 
                  colnames(X[[q]])))
            }
        }
        names(newdataA) = names(X)
        newdata = newdataA
        if (!is.null(object$nzv)) {
            newdata = lapply(1:(length(object$nzv) - 1), function(x) {
                if (length(object$nzv[[x]]$Position > 0)) {
                  newdata[[x]][, -object$nzv[[x]]$Position, drop = FALSE]
                }
                else {
                  newdata[[x]]
                }
            })
        }
        if (length(newdata) != length(object$X)) 
            stop("'newdata' must have as many blocks as 'object$X'")
        names(newdata) = names(X)
        if (any(lapply(newdata, function(x) {
            length(dim(x))
        }) != 2)) {
            if (any(unlist(lapply(newdata, ncol)) != unlist(p))) 
                stop("'newdata' must be a list with ", length(p), 
                  " numeric matrix and ncol respectively equal to ", 
                  paste(p, collapse = ", "), ".")
        }
        if (any(lapply(newdata, function(x) {
            length(dim(x))
        }) == 0)) {
            if (any(unlist(lapply(newdata, ncol)) != unlist(p))) 
                stop("'newdata' must be a list with ", length(p), 
                  " numeric matrix and ncol respectively equal to ", 
                  paste(p, collapse = ", "), ".")
            dim(newdata) = c(1, p)
        }
        for (q in 1:length(newdata)) {
            check = Check.entry.single(newdata[[q]], ncomp[q], 
                q = q)
            newdata[[q]] = check$X
        }
        names(newdata) = names(X)
        if (all.equal(lapply(newdata, colnames), lapply(X, colnames)) != 
            TRUE) 
            stop("Each 'newdata[[i]]' must include all the variables of 'object$X[[i]]'")
        if (!is.null(object$indY)) {
            indY = object$indY
            object$variates = c(object$variates[-indY], object$variates[indY])
            object$loadings = c(object$loadings[-indY], object$loadings[indY])
        }
    }
    if (!is.null(object$logratio)) 
        newdata = lapply(newdata, logratio.transfo, logratio = object$logratio)
    if (!is.null(multilevel)) 
        newdata = lapply(newdata, withinVariation, design = data.frame(multilevel))
    p = lapply(X, ncol)
    q = ncol(Y)
    J = length(X)
    variatesX = object$variates[-(J + 1)]
    loadingsX = object$loadings[-(J + 1)]
    scale = object$scale
    if (!hasArg(newdata.scale)) {
        if (!is(object, mint.object)) {
            if (!is.null(attr(X[[1]], "scaled:center"))) 
                newdata[which(!is.na(ind.match))] = lapply(which(!is.na(ind.match)), 
                  function(x) {
                    sweep(newdata[[x]], 2, STATS = attr(X[[x]], 
                      "scaled:center"))
                  })
            if (scale) 
                newdata[which(!is.na(ind.match))] = lapply(which(!is.na(ind.match)), 
                  function(x) {
                    sweep(newdata[[x]], 2, FUN = "/", STATS = attr(X[[x]], 
                      "scaled:scale"))
                  })
            means.Y = matrix(attr(Y, "scaled:center"), nrow = nrow(newdata[[1]]), 
                ncol = q, byrow = TRUE)
            if (scale) {
                sigma.Y = matrix(attr(Y, "scaled:scale"), nrow = nrow(newdata[[1]]), 
                  ncol = q, byrow = TRUE)
            }
            else {
                sigma.Y = matrix(1, nrow = nrow(newdata[[1]]), 
                  ncol = q)
            }
            concat.newdata = newdata
            names(concat.newdata) = names(X)
        }
        else {
            if (missing(study.test)) {
                if (nlevels(object$study) == 1) {
                  study.test = factor(rep(1, nrow(newdata[[1]])))
                }
                else {
                  stop("'study.test' is missing")
                }
            }
            else {
                study.test = as.factor(study.test)
            }
            if (any(unlist(lapply(newdata, nrow)) != length(study.test))) 
                stop(paste0("'study' must be a factor of length ", 
                  nrow(newdata[[1]]), "."))
            M = nlevels(study.test)
            study.learn = factor(object$study)
            names.study.learn = levels(study.learn)
            names.study.test = levels(study.test)
            match.study = match(names.study.test, names.study.learn)
            match.study.indice = which(!is.na(match.study))
            newdata.list.study = lapply(newdata, study_split, 
                study.test)
            newdata.list.study.scale.temp = NULL
            concat.newdata = vector("list", length = J)
            for (j in 1:J) {
                for (m in 1:M) {
                  if (m %in% match.study.indice) {
                    if (scale == TRUE) {
                      if (nlevels(object$study) > 1) {
                        newdata.list.study.scale.temp = scale(newdata.list.study[[j]][[m]], 
                          center = attr(X[[j]], paste0("means:", 
                            levels(study.test)[m])), scale = attr(X[[j]], 
                            paste0("sigma:", levels(study.test)[m])))
                      }
                      else {
                        newdata.list.study.scale.temp = scale(newdata.list.study[[j]][[m]], 
                          center = attr(X[[j]], "scaled:center"), 
                          scale = attr(X[[j]], "scaled:scale"))
                      }
                    }
                    if (scale == FALSE) {
                      if (nlevels(object$study) > 1) {
                        newdata.list.study.scale.temp = scale(newdata.list.study[[j]][[m]], 
                          center = attr(X[[j]], paste0("means:", 
                            levels(study.test)[m])), scale = FALSE)
                      }
                      else {
                        newdata.list.study.scale.temp = scale(newdata.list.study[[j]][[m]], 
                          center = attr(X[[j]], "scaled:center"), 
                          scale = FALSE)
                      }
                    }
                  }
                  else {
                    if (scale == TRUE & dim(newdata.list.study[[j]][[m]])[1] == 
                      1) {
                      warning("Train data are scaled but the test data include a single sample from a new study (which cannot be scaled). Mkaing prediction without scaling this test sample and thus prediction for this sample should be given with care. Consider scale=FALSE for model, or using more samples for prediction, or using a sample from model studies.\n")
                      newdata.list.study.scale.temp = newdata.list.study[[j]][[m]]
                    }
                    else {
                      newdata.list.study.scale.temp = scale(newdata.list.study[[j]][[m]], 
                        center = TRUE, scale = scale)
                    }
                  }
                  concat.newdata[[j]] = rbind(concat.newdata[[j]], 
                    unlist(newdata.list.study.scale.temp))
                }
            }
            for (j in 1:J) {
                indice.match = match(rownames(newdata[[j]]), 
                  rownames(concat.newdata[[j]]))
                concat.newdata[[j]] = concat.newdata[[j]][indice.match, 
                  , drop = FALSE]
                concat.newdata[[j]][which(is.na(concat.newdata[[j]]))] = 0
            }
            names(concat.newdata) = names(X)
            means.Y = matrix(0, nrow = nrow(concat.newdata[[1]]), 
                ncol = q)
            sigma.Y = matrix(1, nrow = nrow(concat.newdata[[1]]), 
                ncol = q)
            for (m in 1:M) {
                if (m %in% match.study.indice) {
                  if (nlevels(object$study) > 1) {
                    means.Y[which(study.test %in% levels(study.learn)[match.study[m]]), 
                      ] = matrix(attr(Y, paste0("means:", levels(study.test)[m])), 
                      nrow = length(which(study.test %in% levels(study.learn)[match.study[m]])), 
                      ncol = q, byrow = TRUE)
                  }
                  else {
                    means.Y[which(study.test %in% levels(study.learn)[match.study[m]]), 
                      ] = matrix(attr(Y, "scaled:center"), nrow = length(which(study.test %in% 
                      levels(study.learn)[match.study[m]])), 
                      ncol = q, byrow = TRUE)
                  }
                  if (scale == TRUE) {
                    if (nlevels(object$study) > 1) {
                      sigma.Y[which(study.test %in% levels(study.learn)[match.study[m]]), 
                        ] = matrix(attr(Y, paste0("sigma:", levels(study.test)[m])), 
                        nrow = length(which(study.test %in% levels(study.learn)[match.study[m]])), 
                        ncol = q, byrow = TRUE)
                    }
                    else {
                      sigma.Y[which(study.test %in% levels(study.learn)[match.study[m]]), 
                        ] = matrix(attr(Y, "scaled:scale"), nrow = length(which(study.test %in% 
                        levels(study.learn)[match.study[m]])), 
                        ncol = q, byrow = TRUE)
                    }
                  }
                }
            }
        }
    }
    else {
        means.Y = matrix(attr(Y, "scaled:center"), nrow = nrow(newdata[[1]]), 
            ncol = q, byrow = TRUE)
        if (scale) {
            sigma.Y = matrix(attr(Y, "scaled:scale"), nrow = nrow(newdata[[1]]), 
                ncol = q, byrow = TRUE)
        }
        else {
            sigma.Y = matrix(1, nrow = nrow(newdata[[1]]), ncol = q)
        }
        concat.newdata = newdata
        names(concat.newdata) = names(X)
    }
    if (time) 
        time2 = proc.time()
    if (time) 
        print("scaling")
    if (time) 
        print(time2 - time1)
    if ((hasArg(misdata.all) & hasArg(is.na.X)) && any(list(...)$misdata.all)) {
        for (j in c(1:J)[list(...)$misdata.all]) X[[j]][list(...)$is.na.X[[j]]] = 0
    }
    else {
        X = lapply(X, function(x) {
            if (anyNA(x)) {
                ind = is.na(x)
                x[ind] = 0
            }
            x
        })
    }
    if ((hasArg(misdata.all) & hasArg(is.na.newdata)) && any(list(...)$misdata.all)) {
        concat.newdata = lapply(1:J, function(q) {
            replace(concat.newdata[[q]], list(...)$is.na.newdata[[q]], 
                0)
        })
    }
    else {
        concat.newdata = lapply(concat.newdata, function(x) {
            if (anyNA(x)) {
                ind = is.na(x)
                x[ind] = 0
            }
            x
        })
    }
    if (any(sapply(concat.newdata, anyNA))) 
        stop("Some missing values are present in the test data")
    Y[is.na(Y)] = 0
    if (time) 
        time3 = proc.time()
    if (time) 
        print("NA")
    if (time) 
        print(time3 - time2)
    B.hat = t.pred = Y.hat = list()
    for (i in 1:J) {
        Pmat = Cmat = Wmat = NULL
        Pmat = crossprod(X[[i]], variatesX[[i]])
        Cmat = crossprod(Y, variatesX[[i]])
        Wmat = loadingsX[[i]]
        Ypred = lapply(1:ncomp[i], function(x) {
            concat.newdata[[i]] %*% Wmat[, 1:x] %*% solve(t(Pmat[, 
                1:x]) %*% Wmat[, 1:x]) %*% t(Cmat)[1:x, ]
        })
        Ypred = sapply(Ypred, function(x) {
            x * sigma.Y + means.Y
        }, simplify = "array")
        Y.hat[[i]] = array(Ypred, c(nrow(newdata[[i]]), ncol(Y), 
            ncomp[i]))
        t.pred[[i]] = concat.newdata[[i]] %*% Wmat %*% solve(t(Pmat) %*% 
            Wmat)
        t.pred[[i]] = matrix(data = sapply(1:ncol(t.pred[[i]]), 
            function(x) {
                t.pred[[i]][, x] * apply(variatesX[[i]], 2, function(y) {
                  (norm(y, type = "2"))^2
                })[x]
            }), nrow = nrow(concat.newdata[[i]]), ncol = ncol(t.pred[[i]]))
        B.hat[[i]] = sapply(1:ncomp[i], function(x) {
            Wmat[, 1:x] %*% solve(t(Pmat[, 1:x]) %*% Wmat[, 1:x]) %*% 
                t(Cmat)[1:x, ]
        }, simplify = "array")
        rownames(t.pred[[i]]) = rownames(newdata[[i]])
        colnames(t.pred[[i]]) = paste0("dim", c(1:ncomp[i]))
        rownames(Y.hat[[i]]) = rownames(newdata[[i]])
        colnames(Y.hat[[i]]) = colnames(Y)
        dimnames(Y.hat[[i]])[[3]] = paste0("dim", c(1:ncomp[i]))
        rownames(B.hat[[i]]) = colnames(newdata[[i]])
        colnames(B.hat[[i]]) = colnames(Y)
        dimnames(B.hat[[i]])[[3]] = paste0("dim", c(1:ncomp[i]))
    }
    if (time) 
        time4 = proc.time()
    if (time) 
        print("Prediction")
    if (time) 
        print(time4 - time3)
    names(Y.hat) = names(t.pred) = names(B.hat) = names(object$X)
    if (time) 
        time4 = proc.time()
    if (is(object, block.object) & length(object$X) > 1) {
        out = list(predict = Y.hat[which(!is.na(ind.match))], 
            variates = t.pred[which(!is.na(ind.match))], B.hat = B.hat[which(!is.na(ind.match))])
        temp.all = list()
        for (comp in 1:min(ncomp[-object$indY])) {
            temp = array(0, c(nrow(Y.hat[[1]]), ncol(Y.hat[[1]]), 
                J), dimnames = list(rownames(newdata[[1]]), colnames(Y), 
                names(object$X)))
            for (i in 1:J) temp[, , i] = Y.hat[[i]][, , comp]
            temp.all[[comp]] = temp
        }
        names(temp.all) = paste0("dim", c(1:min(ncomp[-object$indY])))
        if (!hasArg(noAveragePredict)) {
            out$AveragedPredict = array(unlist(lapply(temp.all, 
                function(x) {
                  apply(x, c(1, 2), mean)
                })), dim(Y.hat[[1]]), dimnames = list(rownames(newdata[[1]]), 
                colnames(Y), paste0("dim", c(1:min(ncomp[-object$indY])))))
            out$WeightedPredict = array(unlist(lapply(temp.all, 
                function(x) {
                  apply(x, c(1, 2), function(z) {
                    temp = aggregate(rowMeans(object$weights), 
                      list(z), sum)
                    ind = which(temp[, 2] == max(temp[, 2]))
                    if (length(ind) == 1) {
                      res = temp[ind, 1]
                    }
                    else {
                      res = NA
                    }
                    res
                  })
                })), dim(Y.hat[[1]]), dimnames = list(rownames(newdata[[1]]), 
                colnames(Y), paste0("dim", c(1:min(ncomp[-object$indY])))))
        }
    }
    else if (is(object, block.object)) {
        out = list(predict = Y.hat, variates = t.pred, B.hat = B.hat)
    }
    else {
        out = list(predict = Y.hat[[1]], variates = t.pred[[1]], 
            B.hat = B.hat[[1]])
    }
    if (time) 
        time5 = proc.time()
    if (time) 
        print("Y.hat")
    if (time) 
        print(time5 - time4)
    if (is(object, "DA")) {
        if (is(object, block.object) & length(object$X) > 1) {
            if (!hasArg(noAveragePredict)) {
                out$AveragedPredict.class$max.dist = matrix(sapply(1:ncomp[1], 
                  function(y) {
                    apply(out$AveragedPredict[, , y, drop = FALSE], 
                      1, function(z) {
                        paste(levels(Y.factor)[which(z == max(z))], 
                          collapse = "/")
                      })
                  }), nrow = nrow(newdata[[1]]), ncol = ncomp[1])
                out$WeightedPredict.class$max.dist = matrix(sapply(1:ncomp[1], 
                  function(y) {
                    apply(out$WeightedPredict[, , y, drop = FALSE], 
                      1, function(z) {
                        paste(levels(Y.factor)[which(z == max(z))], 
                          collapse = "/")
                      })
                  }), nrow = nrow(newdata[[1]]), ncol = ncomp[1])
                rownames(out$AveragedPredict.class$max.dist) = rownames(out$WeightedPredict.class$max.dist) = rownames(newdata[[1]])
                colnames(out$AveragedPredict.class$max.dist) = colnames(out$WeightedPredict.class$max.dist) = paste0("dim", 
                  c(1:min(ncomp[-object$indY])))
            }
        }
        out.temp = list(predict = Y.hat[which(!is.na(ind.match))], 
            variates = t.pred[which(!is.na(ind.match))], B.hat = B.hat[which(!is.na(ind.match))])
        out.temp$newdata = concat.newdata[which(!is.na(ind.match))]
        object.temp = object
        object.temp$X = object.temp$X[which(!is.na(ind.match))]
        object.temp$variates = object.temp$variates[c(which(!is.na(ind.match)), 
            J + 1)]
        if (!is.null(object$weights)) 
            weights <- rowMeans(object$weights)[which(!is.na(ind.match))]
        else weights <- NULL
        classif.DA = internal_predict.DA(object = object.temp, 
            q = q, out = out.temp, dist = dist, weights = weights)
        out = c(out, classif.DA)
    }
    if (time) 
        time6 = proc.time()
    if (time) 
        print("DA")
    if (time) 
        print(time6 - time5)
    out$call = match.call()
    class(out) = paste("predict")
    out
}

Check.entry.single = function(X,  ncomp, q)
{

    #-- validation des arguments --#
    if (length(dim(X)) != 2)
    stop(paste0("'X[[", q, "]]' must be a numeric matrix."))
    
    if(! any(class(X) %in% "matrix"))
    X = as.matrix(X)
    
    if (!is.numeric(X))
    stop(paste0("'X[[", q, "]]'  must be a numeric matrix."))
    
    N = nrow(X)
    P = ncol(X)
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
    stop(paste0("invalid number of variates 'ncomp' for matrix 'X[[", q, "]]'."))
    
    ncomp = round(ncomp)
    
    # add colnames and rownames if missing
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:P, sep = "")
        dimnames(X)[[2]] = X.names
    }
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = 1:N
        rownames(X)  = ind.names
    }
    
    if (length(unique(rownames(X))) != nrow(X))
        stop("samples should have a unique identifier/rowname")
    if (length(unique(X.names)) != P)
        stop("Unique indentifier is needed for the columns of X")
    
    return(list(X=X, ncomp=ncomp, X.names=X.names, ind.names=ind.names))
}

internal_predict.DA<-
function (object, out, q, dist, weights) 
{
    if (!is(object, "DA")) 
        stop("'Object' is not from a Discriminant Analysis", 
            call. = FALSE)
    out.DA = list()
    J = length(object$X)
    p = lapply(object$X, ncol)
    t.pred = out$variates
    Y.hat = out$predict
    newdata = out$newdata
    variatesX = object$variates[-(J + 1)]
    ncomp = object$ncomp
    Y = object$Y
    Y.prim = unmap(object$Y)
    G = cls = list()
    for (i in 1:J) {
        G[[i]] = sapply(1:q, function(x) {
            apply(as.matrix(variatesX[[i]][Y.prim[, x] == 1, 
                , drop = FALSE]), 2, mean)
        })
        if (ncomp[i] == 1) 
            G[[i]] = t(t(G[[i]]))
        else G[[i]] = t(G[[i]])
        colnames(G[[i]]) = paste0("dim", c(1:ncomp[i]))
        rownames(G[[i]]) = colnames(object$ind.mat)
    }
    names(G) = names(object$X)
    if (any(dist == "all") || any(dist == "max.dist")) {
        cls$max.dist = lapply(1:J, function(x) {
            matrix(sapply(1:ncomp[x], function(y) {
                apply(Y.hat[[x]][, , y, drop = FALSE], 1, function(z) {
                  paste(levels(Y)[which(z == max(z))], collapse = "/")
                })
            }), nrow = nrow(newdata[[x]]), ncol = ncomp[x])
        })
        cls$max.dist = lapply(1:J, function(x) {
            colnames(cls$max.dist[[x]]) = paste0(rep("comp", 
                ncomp[x]), 1:ncomp[[x]])
            rownames(cls$max.dist[[x]]) = rownames(newdata[[x]])
            return(cls$max.dist[[x]])
        })
        names(cls$max.dist) = names(object$X)
    }
    if (any(dist == "all") || any(dist == "centroids.dist")) {
        cl = list()
        centroids.fun = function(x, G, h, i) {
            q = nrow(G[[i]])
            x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
            if (h > 1) {
                d = apply((x - G[[i]][, 1:h])^2, 1, sum)
            }
            else {
                d = (x - G[[i]][, 1])^2
            }
            cl.id = paste(levels(Y)[which(d == min(d))], collapse = "/")
        }
        for (i in 1:J) {
            cl[[i]] = matrix(nrow = nrow(newdata[[1]]), ncol = ncomp[i])
            for (h in 1:ncomp[[i]]) {
                cl.id = apply(matrix(t.pred[[i]][, 1:h], ncol = h), 
                  1, function(x) {
                    centroids.fun(x = x, G = G, h = h, i = i)
                  })
                cl[[i]][, h] = cl.id
            }
        }
        cls$centroids.dist = lapply(1:J, function(x) {
            colnames(cl[[x]]) = paste0(rep("comp", ncomp[x]), 
                1:ncomp[[x]])
            rownames(cl[[x]]) = rownames(newdata[[x]])
            return(cl[[x]])
        })
        names(cls$centroids.dist) = names(object$X)
    }
    if (any(dist == "all") || any(dist == "mahalanobis.dist")) {
        cl = list()
        Sr.fun = function(x, G, Yprim, h, i) {
            q = nrow(G[[i]])
            Xe = Yprim %*% G[[i]][, 1:h]
            Xr = variatesX[[i]][, 1:h] - Xe
            Sr = t(Xr) %*% Xr/nrow(Yprim)
            Sr.inv = solve(Sr)
            x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
            if (h > 1) {
                mat = (x - G[[i]][, 1:h]) %*% Sr.inv %*% t(x - 
                  G[[i]][, 1:h])
                d = apply(mat^2, 1, sum)
            }
            else {
                d = drop(Sr.inv) * (x - G[[i]][, 1])^2
            }
            cl.id = paste(levels(Y)[which(d == min(d))], collapse = "/")
        }
        for (i in 1:J) {
            cl[[i]] = matrix(nrow = nrow(newdata[[1]]), ncol = ncomp[i])
            for (h in 1:ncomp[[i]]) {
                cl.id = apply(matrix(t.pred[[i]][, 1:h], ncol = h), 
                  1, Sr.fun, G = G, Yprim = Y.prim, h = h, i = i)
                cl[[i]][, h] = cl.id
            }
        }
        cls$mahalanobis.dist = lapply(1:J, function(x) {
            colnames(cl[[x]]) = paste0(rep("comp", ncomp[x]), 
                1:ncomp[[x]])
            rownames(cl[[x]]) = rownames(newdata[[x]])
            return(cl[[x]])
        })
        names(cls$mahalanobis.dist) = names(object$X)
    }
    out.DA$class = cls
    if (length(object$X) > 1) {
        for (ijk in 1:length(out.DA$class)) {
            temp = array(, c(nrow(newdata[[1]]), min(ncomp), 
                J))
            for (i in 1:J) {
                temp[, , i] = out.DA$class[[ijk]][[i]][, 1:min(ncomp)]
            }
            table.temp = apply(temp, c(1, 2), function(x) {
                a = table(x)
                if (length(which(a == max(a))) == 1) {
                  b = names(which.max(a))
                }
                else {
                  b = NA
                }
                b
            })
            colnames(table.temp) = colnames(out.DA$class[[ijk]][[i]])[1:min(ncomp)]
            rownames(table.temp) = rownames(out.DA$class[[ijk]][[i]])
            out.DA$MajorityVote[[ijk]] = table.temp
        }
        names(out.DA$MajorityVote) = names(out.DA$class)
        if (!is.null(weights)) {
            out.WeightedVote = vector("list", length = length(out.DA$class))
            Group.2 = n = x = NULL
            for (i in 1:length(out.DA$class)) {
                out = matrix(NA_real_, nrow = nrow(newdata[[1]]), 
                  ncol = min(ncomp))
                rownames(out) = rownames(newdata[[1]])
                colnames(out) = paste0("comp", 1:min(ncomp))
                for (comp in 1:min(ncomp)) {
                  data.temp = NULL
                  for (j in 1:J) {
                    data.temp = rbind(data.temp, out.DA$class[[i]][[j]][, 
                      comp, drop = FALSE])
                  }
                  colnames(data.temp) = "pred"
                  temp = data.frame(data.temp, indiv = rownames(data.temp), 
                    weights = rep(weights, each = nrow(out.DA$class[[1]][[1]])))
                  ag = aggregate(temp$weights, by = list(temp$pred, 
                    temp$indiv), FUN = sum)
                  data_max <- group_by(.data = ag, Group.2)
                  data_max <- filter(.data = data_max, row_number(x) == 
                    n())
                  out.comp = as.matrix(data_max[, 1])
                  rownames(out.comp) = as.matrix(data_max[, 2])
                  colnames(out.comp) = paste0("comp", comp)
                  out[, comp] = out.comp[match(rownames(out), 
                    rownames(out.comp)), ]
                }
                out.WeightedVote[[i]] = out
            }
            names(out.WeightedVote) = names(out.DA$class)
            out.DA$WeightedVote = out.WeightedVote
            if (FALSE) {
                out.DA$WeightedVote = lapply(out.DA$class, function(x) {
                  class.per.comp = lapply(1:min(ncomp), function(y) {
                    matrix(sapply(x, function(z) z[, y, drop = FALSE]), 
                      ncol = J)
                  })
                  names(class.per.comp) = paste0("comp", 1:min(ncomp))
                  class.weighted.per.comp = sapply(class.per.comp, 
                    function(y) {
                      apply(y, 1, function(z) {
                        temp = aggregate(weights, list(as.character(z)), 
                          sum)
                        ind = which(temp[, 2] == max(temp[, 2]))
                        if (length(ind) == 1) {
                          res = as.character(temp[ind, 1])
                        }
                        else {
                          res = NA
                        }
                        res
                      })
                    })
                  class.weighted.per.comp = matrix(class.weighted.per.comp, 
                    nrow = nrow(class.per.comp[[1]]))
                  colnames(class.weighted.per.comp) = names(class.per.comp)
                  rownames(class.weighted.per.comp) = rownames(out.DA$MajorityVote[[1]])
                  class.weighted.per.comp
                })
            }
            out.DA$weights = weights
        }
    }
    else {
        out.DA$MajorityVote = lapply(out.DA$class, function(x) {
            x[[1]]
        })
    }
    block.object = c("block.pls", "block.spls", "block.plsda", 
        "block.spsda")
    if (is(object, block.object) & J > 1) {
        out.DA$centroids = G
    }
    else {
        out.DA$centroids = G[[1]]
        out.DA$class = out.DA$MajorityVote
    }
    if (any(dist == "all")) 
        dist = "all"
    out.DA$dist = dist
    out.DA
}

.unlist_repeat_cv_output<-
function (list_rep = NULL) 
{
    names_list <- as.list(names(list_rep[[1]]))
    names(names_list) <- names(list_rep[[1]])
    list_nrep <- lapply(names_list, function(x) {
        lapply(list_rep, function(y) {
            y[[x]]
        })
    })
    return(list_nrep)
}

t.test.process<-
function (mat.error.rate, alpha = 0.01, alternative = "greater") 
{
    if (nrow(mat.error.rate) < 3) {
        return(NA)
    }
    .calc_pval <- function(df, col1, col2) {
        x <- df[, col1]
        y <- df[, col2]
        pval <- tryCatch({
            t.test(x, y, alternative = alternative)$p.value
        }, error = function(e) e)
        if (!is.numeric(pval)) {
            if (is(pval, "error") && grepl("data are essentially constant", 
                x = pval$message)) {
                if (mean(y)/(mean(x) + 0.01) < 0.95) {
                  pval <- 0
                }
                else {
                  pval <- 1
                }
            }
            else {
                .unexpected_err(trying_to = "choose the optimum number of components")
            }
        }
        else if (is.nan(pval)) {
            pval <- 1
        }
        return(pval)
    }
    max_comp <- ncol(mat.error.rate)
    pval <- 1
    opt_comp <- 1
    next_comp <- 2
    while (opt_comp < max_comp & next_comp <= max_comp) {
        pval <- .calc_pval(mat.error.rate, opt_comp, next_comp)
        if (pval < alpha) {
            opt_comp <- next_comp
        }
        next_comp <- next_comp + 1
    }
    return(opt_comp)
}

