avgNormDiff <-
function (x) 
{
    normDiff <- function(x) {
        nm <- colnames(x)
        x <- as.numeric(x)
        n <- length(x)
        rngX <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
        out <- matrix(NA, nrow = n, ncol = n)
        for (j in seq_len(n)) out[, j] <- 1 - (abs((x - x[j])/rngX))
        rownames(out) <- nm
        colnames(out) <- nm
        out
    }
    sim <- matrix(0, nrow = ncol(x), ncol = ncol(x))
    for (k in seq_len(nrow(x))) {
        tmp <- normDiff(x[k, , drop = FALSE])
        sim <- sim + tmp
        rownames(sim) <- rownames(tmp)
        colnames(sim) <- colnames(tmp)
    }
    sim <- sim/nrow(x)
    sim
}
buildPredictor <-
function (dataList, groupList, outDir = tempdir(), makeNetFunc = NULL, 
    sims = NULL, featScoreMax = 10L, trainProp = 0.80000000000000004, 
    numSplits = 10L, numCores, JavaMemory = 4L, featSelCutoff = 9L, 
    keepAllData = FALSE, startAt = 1L, preFilter = FALSE, impute = FALSE, 
    preFilterGroups = NULL, imputeGroups = NULL, logging = "default", 
    debugMode = FALSE) 
{
    verbose_default <- TRUE
    verbose_runQuery <- FALSE
    verbose_compileNets <- FALSE
    verbose_runFS <- TRUE
    verbose_predict <- FALSE
    verbose_compileFS <- FALSE
    verbose_makeFeatures <- FALSE
    if (logging == "all") {
        verbose_runQuery <- TRUE
        verbose_compileNets <- TRUE
        verbose_compileFS <- TRUE
        verbose_makeFeatures <- TRUE
    }
    else if (logging == "none") {
        verbose_runFS <- FALSE
        verbose_default <- FALSE
        verbose_predict <- FALSE
    }
    if (missing(dataList)) 
        stop("dataList must be supplied.\n")
    if (missing(groupList)) 
        stop("groupList must be supplied.\n")
    if (length(groupList) < 1) 
        stop("groupList must be of length 1+\n")
    tmp <- unlist(lapply(groupList, class))
    not_list <- sum(tmp == "list") < length(tmp)
    nm1 <- setdiff(names(groupList), "clinical")
    if (!is(dataList, "MultiAssayExperiment")) 
        stop("dataList must be a MultiAssayExperiment")
    names_nomatch <- any(!nm1 %in% names(dataList))
    if (!is(groupList, "list") || not_list || names_nomatch) {
        msg <- c("groupList must be a list of lists.", " Names must match those in dataList, and each entry should be a list", 
            " of networks for this group.")
        stop(paste(msg, sep = ""))
    }
    x <- checkMakeNetFuncSims(makeNetFunc = makeNetFunc, sims = sims, 
        groupList = groupList)
    if (!is(dataList, "MultiAssayExperiment")) 
        stop("dataList must be a MultiAssayExperiment")
    if (trainProp <= 0 | trainProp >= 1) 
        stop("trainProp must be greater than 0 and less than 1")
    if (startAt > numSplits) 
        stop("startAt should be between 1 and numSplits")
    megaDir <- outDir
    if (file.exists(megaDir)) {
        stop(paste("outDir seems to already exist!", "Please provide a new directory, as its contents will be overwritten", 
            sprintf("You provided: %s", outDir), sep = "\n"))
    }
    else {
        dir.create(megaDir, recursive = TRUE)
    }
    pheno_all <- colData(dataList)
    pheno_all <- as.data.frame(pheno_all)
    message("Predictor started at:")
    message(Sys.time())
    subtypes <- unique(pheno_all$STATUS)
    exprs <- experiments(dataList)
    datList2 <- list()
    for (k in seq_len(length(exprs))) {
        tmp <- exprs[[k]]
        df <- sampleMap(dataList)[which(sampleMap(dataList)$assay == 
            names(exprs)[k]), ]
        colnames(tmp) <- df$primary[match(df$colname, colnames(tmp))]
        if ("matrix" %in% class(tmp)) {
            datList2[[names(exprs)[k]]] <- tmp
        }
        else {
            tmp <- as.matrix(assays(tmp)[[1]])
            datList2[[names(exprs)[k]]] <- tmp
        }
    }
    if ("clinical" %in% names(groupList)) {
        tmp <- colData(dataList)
        vars <- unique(unlist(groupList[["clinical"]]))
        datList2[["clinical"]] <- t(as.matrix(tmp[, vars, drop = FALSE]))
    }
    dataList <- datList2
    rm(datList2)
    if (verbose_default) {
        message(sprintf("-------------------------------"))
        message(sprintf("# patients = %i", nrow(pheno_all)))
        message(sprintf("# classes = %i { %s }", length(subtypes), 
            paste(subtypes, collapse = ",")))
        message("Sample breakdown by class")
        message(table(pheno_all$STATUS))
        message(sprintf("%i train/test splits", numSplits))
        message(sprintf("Feature selection cutoff = %i of %i", 
            featSelCutoff, featScoreMax))
        message(sprintf("Datapoints:"))
        for (nm in names(dataList)) {
            message(sprintf("\t%s: %i units", nm, nrow(dataList[[nm]])))
        }
    }
    outList <- list()
    tmp <- list()
    for (nm in names(groupList)) {
        curNames <- names(groupList[[nm]])
        tmp[[nm]] <- cbind(rep(nm, length(curNames)), curNames)
    }
    tmp <- do.call("rbind", tmp)
    if (length(nm) < 2) 
        tmp <- as.matrix(tmp)
    colnames(tmp) <- c("NetType", "NetName")
    outList[["inputNets"]] <- tmp
    if (verbose_default) {
        if (!is.null(makeNetFunc)) {
            message("\n\nCustom function to generate input nets:")
            print(makeNetFunc)
        }
        else {
            message("Similarity metrics provided:")
            print(sims)
        }
        message(sprintf("-------------------------------\n"))
    }
    for (rngNum in startAt:numSplits) {
        curList <- list()
        if (verbose_default) {
            message(sprintf("-------------------------------"))
            message(sprintf("Train/test split # %i", rngNum))
            message(sprintf("-------------------------------"))
        }
        outDir <- paste(megaDir, sprintf("rng%i", rngNum), sep = "/")
        dir.create(outDir)
        pheno_all$TT_STATUS <- splitTestTrain(pheno_all, pctT = trainProp, 
            verbose = verbose_default)
        pheno <- pheno_all[which(pheno_all$TT_STATUS %in% "TRAIN"), 
            ]
        dats_train <- lapply(dataList, function(x) x[, which(colnames(x) %in% 
            pheno$ID), drop = FALSE])
        if (impute) {
            if (verbose_default) 
                message("**** IMPUTING ****")
            if (is.null(imputeGroups)) 
                imputeGroups <- names(dats_train)
            if (!any(imputeGroups %in% names(dats_train))) 
                stop("imputeGroups must match names in dataList")
            nmset <- names(dats_train)
            dats_train <- lapply(names(dats_train), function(nm) {
                x <- dats_train[[nm]]
                print(class(x))
                if (nm %in% imputeGroups) {
                  missidx <- which(rowSums(is.na(x)) > 0)
                  for (i in missidx) {
                    na_idx <- which(is.na(x[i, ]))
                    x[i, na_idx] <- median(x[i, ], na.rm = TRUE)
                  }
                }
                x
            })
            names(dats_train) <- nmset
        }
        if (preFilter) {
            if (is.null(preFilterGroups)) 
                preFilterGroups <- names(dats_train)
            if (!any(preFilterGroups %in% names(dats_train))) {
                stop("preFilterGroups must match names in dataList")
            }
            message("Prefiltering enabled")
            for (nm in preFilterGroups) {
                message(sprintf("%s: %i variables", nm, nrow(dats_train[[nm]])))
                if (nrow(dats_train[[nm]]) < 2) 
                  vars <- rownames(dats_train[[nm]])
                else {
                  newx <- na.omit(dats_train[[nm]])
                  tmp <- pheno[which(pheno$ID %in% colnames(newx)), 
                    ]
                  tryCatch({
                    fit <- cv.glmnet(x = t(newx), y = factor(tmp$STATUS), 
                      family = "binomial", alpha = 1)
                  }, error = function(ex) {
                    print(ex)
                    message("*** You may need to set impute=TRUE for prefiltering ***")
                  }, finally = {
                  })
                  wt <- abs(coef(fit, s = "lambda.min")[, 1])
                  vars <- setdiff(names(wt)[which(wt > .Machine$double.eps)], 
                    "(Intercept)")
                }
                if (length(vars) > 0) {
                  tmp <- dats_train[[nm]]
                  tmp <- tmp[which(rownames(tmp) %in% vars), 
                    , drop = FALSE]
                  dats_train[[nm]] <- tmp
                }
                else {
                }
                message(sprintf("rngNum %i: %s: %s pruned", rngNum, 
                  nm, length(vars)))
            }
        }
        if (verbose_default) {
            message("# values per feature (training)")
            for (nm in names(dats_train)) {
                message(sprintf("\tGroup %s: %i values", nm, 
                  nrow(dats_train[[nm]])))
            }
        }
        netDir <- paste(outDir, "tmp", sep = "/")
        dir.create(netDir)
        pheno_id <- setupFeatureDB(pheno, netDir)
        if (verbose_default) 
            message("** Creating features")
        createPSN_MultiData(dataList = dats_train, groupList = groupList, 
            pheno = pheno_id, netDir = netDir, makeNetFunc = makeNetFunc, 
            sims = sims, numCores = numCores, verbose = verbose_makeFeatures)
        if (verbose_default) 
            message("** Compiling features")
        dbDir <- compileFeatures(netDir, outDir, numCores = numCores, 
            verbose = verbose_compileFS, debugMode = debugMode)
        if (verbose_default) 
            message("\n** Running feature selection")
        curList[["featureScores"]] <- list()
        for (g in subtypes) {
            pDir <- paste(outDir, g, sep = "/")
            if (file.exists(pDir)) 
                unlink(pDir, recursive = TRUE)
            dir.create(pDir)
            if (verbose_default) 
                message(sprintf("\tClass: %s", g))
            pheno_subtype <- pheno
            pheno_subtype$STATUS[which(!pheno_subtype$STATUS %in% 
                g)] <- "nonpred"
            trainPred <- pheno_subtype$ID[which(pheno_subtype$STATUS %in% 
                g)]
            if (verbose_default) {
                print(table(pheno_subtype$STATUS, useNA = "always"))
            }
            resDir <- paste(pDir, "GM_results", sep = "/")
            message(sprintf("\tScoring features"))
            runFeatureSelection(trainPred, outDir = resDir, dbPath = dbDir$dbDir, 
                nrow(pheno_subtype), verbose = verbose_runFS, 
                numCores = numCores, verbose_runQuery = TRUE, 
                featScoreMax = featScoreMax, JavaMemory = JavaMemory, 
                debugMode = debugMode)
            tmp <- dir(path = resDir, pattern = "RANK$")[1]
            tmp <- sprintf("%s/%s", resDir, tmp)
            if (sum(grepl(pattern = ",", readLines(tmp, n = 6)) > 
                0)) {
                replacePattern(path = resDir, fileType = "RANK$")
            }
            nrank <- dir(path = resDir, pattern = "NRANK$")
            if (verbose_default) 
                message("\tCompiling feature scores")
            pTally <- compileFeatureScores(paste(resDir, nrank, 
                sep = "/"), verbose = verbose_compileFS)
            tallyFile <- paste(resDir, sprintf("%s_pathway_CV_score.txt", 
                g), sep = "/")
            write.table(pTally, file = tallyFile, sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE)
            curList[["featureScores"]][[g]] <- pTally
            if (verbose_default) 
                message("")
        }
        if (verbose_default) 
            message("\n** Predicting labels for test")
        pheno <- pheno_all
        predRes <- list()
        curList[["featureSelected"]] <- list()
        for (g in subtypes) {
            if (verbose_default) 
                message(sprintf("%s", g))
            pDir <- paste(outDir, g, sep = "/")
            pTally <- read.delim(paste(pDir, "GM_results", sprintf("%s_pathway_CV_score.txt", 
                g), sep = "/"), sep = "\t", header = TRUE, 
                as.is = TRUE)
            idx <- which(pTally[, 2] >= featSelCutoff)
            pTally <- pTally[idx, 1]
            pTally <- sub(".profile", "", pTally)
            pTally <- sub("_cont.txt", "", pTally)
            curList[["featureSelected"]][[g]] <- pTally
            if (verbose_default) 
                message(sprintf("\t%i feature(s) selected", length(pTally)))
            netDir <- paste(pDir, "networks", sep = "/")
            dats_tmp <- list()
            for (nm in names(dataList)) {
                passed <- rownames(dats_train[[nm]])
                tmp <- dataList[[nm]]
                dats_tmp[[nm]] <- tmp[which(rownames(tmp) %in% 
                  passed), , drop = FALSE]
            }
            if (impute) {
                train_samp <- pheno_all$ID[which(pheno_all$TT_STATUS %in% 
                  "TRAIN")]
                test_samp <- pheno_all$ID[which(pheno_all$TT_STATUS %in% 
                  "TEST")]
                nmSet <- names(dats_tmp)
                dats_tmp <- lapply(names(dats_tmp), function(nm) {
                  x <- dats_tmp[[nm]]
                  if (nm %in% imputeGroups) {
                    missidx <- which(rowSums(is.na(x)) > 0)
                    train_idx <- which(colnames(x) %in% train_samp)
                    test_idx <- which(colnames(x) %in% test_samp)
                    for (i in missidx) {
                      na_idx <- intersect(which(is.na(x[i, ])), 
                        train_idx)
                      na_idx1 <- na_idx
                      x[i, na_idx] <- median(x[i, train_idx], 
                        na.rm = TRUE)
                      na_idx <- intersect(which(is.na(x[i, ])), 
                        test_idx)
                      na_idx2 <- na_idx
                      x[i, na_idx] <- median(x[i, test_idx], 
                        na.rm = TRUE)
                    }
                  }
                  x
                })
                names(dats_tmp) <- nmSet
            }
            if (verbose_default) 
                message(sprintf("\t%s: Create & compile features", 
                  g))
            if (length(pTally) >= 1) {
                netDir <- paste(pDir, "tmp", sep = "/")
                dir.create(netDir)
                pheno_id <- setupFeatureDB(pheno, netDir)
                createPSN_MultiData(dataList = dats_tmp, groupList = groupList, 
                  pheno = pheno_id, netDir = netDir, makeNetFunc = makeNetFunc, 
                  sims = sims, numCores = numCores, filterSet = pTally, 
                  verbose = verbose_default)
                dbDir <- compileFeatures(netDir, outDir = pDir, 
                  numCores = numCores, verbose = verbose_compileNets, 
                  debugMode = debugMode)
                qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
                  pheno$TT_STATUS %in% "TRAIN")]
                qFile <- paste(pDir, sprintf("%s_query", g), 
                  sep = "/")
                writeQueryFile(qSamps, "all", nrow(pheno), qFile)
                if (verbose_default) 
                  message(sprintf("\t** %s: Compute similarity", 
                    g))
                resFile <- runQuery(dbDir$dbDir, qFile, resDir = pDir, 
                  JavaMemory = JavaMemory, numCores = numCores, 
                  verbose = verbose_runQuery, debugMode = debugMode)
                predRes[[g]] <- getPatientRankings(sprintf("%s.PRANK", 
                  resFile), pheno, g)
            }
            else {
                predRes[[g]] <- NA
            }
        }
        if (verbose_default) 
            message("")
        if (sum(is.na(predRes)) > 0 & verbose_default) {
            str <- sprintf("RNG %i : One or more classes have no selected features.", 
                rngNum)
            str <- sprintf("%s Not classifying.", str)
            message(str)
        }
        else {
            if (verbose_default) 
                message("** Predict labels")
            predClass <- predictPatientLabels(predRes, verbose = verbose_predict)
            out <- merge(x = pheno_all, y = predClass, by = "ID")
            outFile <- paste(outDir, "predictionResults.txt", 
                sep = "/")
            acc <- sum(out$STATUS == out$PRED_CLASS)/nrow(out)
            if (verbose_default) 
                message(sprintf("Split %i: ACCURACY (N=%i test) = %2.1f%%", 
                  rngNum, nrow(out), acc * 100))
            curList[["predictions"]] <- out
            curList[["accuracy"]] <- acc
        }
        if (!keepAllData) {
            unlink(outDir, recursive = TRUE)
        }
        if (verbose_default) {
            message("\n----------------------------------------")
        }
        outList[[sprintf("Split%i", rngNum)]] <- curList
    }
    message("Predictor completed at:")
    message(Sys.time())
    return(outList)
}
buildPredictor_sparseGenetic <-
function (phenoDF, cnv_GR, predClass, group_GRList, outDir = tempdir(), 
    numSplits = 3L, featScoreMax = 10L, filter_WtSum = 100L, 
    enrichLabels = TRUE, enrichPthresh = 0.070000000000000007, 
    numPermsEnrich = 2500L, minEnr = -1, numCores = 1L, FS_numCores = NULL, 
    ...) 
{
    if (file.exists(outDir)) {
        stop("outDir exists. Please specify another output directory.")
    }
    dir.create(outDir, recursive = TRUE)
    netDir <- paste(outDir, "networks_orig", sep = "/")
    message("making rangesets")
    netList <- makePSN_RangeSets(cnv_GR, group_GRList, netDir, 
        verbose = FALSE)
    message("counting patients in net")
    p <- countPatientsInNet(netDir, netList, phenoDF$ID)
    message("updating nets")
    tmp <- updateNets(p, phenoDF, writeNewNets = FALSE, verbose = FALSE)
    netmat <- tmp[[1]]
    phenoDF <- tmp[[2]]
    if (is.null(FS_numCores)) 
        FS_numCores <- max(1, numCores - 1)
    message("* Resampling train/test samples")
    TT_STATUS <- splitTestTrain_resampling(phenoDF, nFold = numSplits, 
        predClass = predClass, verbose = TRUE)
    p_full <- netmat
    pheno_full <- phenoDF
    if (any(colnames(pheno_full) %in% "TT_STATUS")) {
        message(paste("** Warning, found TT_STATUS column. ", 
            "netDx adds its own column so this one will be removed **", 
            sep = ""))
        pheno_full <- pheno_full[, -which(colnames(pheno_full) %in% 
            "TT_STATUS")]
    }
    pScore <- list()
    enrichedNets <- list()
    for (k in seq_len(length(TT_STATUS))) {
        p <- p_full
        pheno <- pheno_full
        pheno <- cbind(pheno, TT_STATUS = TT_STATUS[[k]])
        message("----------------------------------------")
        message(sprintf("Resampling round %i", k))
        message("----------------------------------------")
        print(table(pheno[, c("STATUS", "TT_STATUS")]))
        newOut <- paste(outDir, sprintf("part%i", k), sep = "/")
        dir.create(newOut)
        outF <- paste(newOut, "TT_STATUS.txt", sep = "/")
        write.table(pheno, file = outF, sep = "\t", col.names = TRUE, 
            row.names = FALSE, quote = FALSE)
        message("# patients: train only")
        pheno_train <- subset(pheno, TT_STATUS %in% "TRAIN")
        p_train <- p[which(rownames(p) %in% pheno_train$ID), 
            ]
        print(nrow(pheno_train))
        message("Training only:")
        trainNetDir <- paste(newOut, "networks", sep = "/")
        tmp <- updateNets(p_train, pheno_train, oldNetDir = netDir, 
            newNetDir = trainNetDir, verbose = FALSE)
        p_train <- tmp[[1]]
        pheno_train <- tmp[[2]]
        if (enrichLabels) {
            message("Running label enrichment")
            tmpDir <- paste(outDir, "tmp", sep = "/")
            if (!file.exists(tmpDir)) 
                dir.create(tmpDir)
            netInfo <- enrichLabelNets(trainNetDir, pheno_train, 
                newOut, predClass = predClass, numReps = numPermsEnrich, 
                numCores = numCores, tmpDir = tmpDir, verbose = FALSE)
            pvals <- as.numeric(netInfo[, "pctl"])
            netInfo <- netInfo[which(pvals < enrichPthresh), 
                ]
            print(nrow(netInfo))
            p_train <- p_train[, which(colnames(p_train) %in% 
                rownames(netInfo))]
            trainNetDir <- paste(newOut, "networksEnriched", 
                sep = "/")
            tmp <- updateNets(p_train, pheno_train, oldNetDir = netDir, 
                newNetDir = trainNetDir, verbose = FALSE)
            p_train <- tmp[[1]]
            pheno_train <- tmp[[2]]
            enrichedNets[[k]] <- rownames(netInfo)
        }
        pheno_train <- setupFeatureDB(pheno_train, trainNetDir)
        moveInteractionNets(netDir = trainNetDir, outDir = paste(trainNetDir, 
            "INTERACTIONS", sep = "/"), pheno = pheno_train)
        tmp <- dir(path = sprintf("%s/INTERACTIONS", trainNetDir), 
            pattern = "txt$")[1]
        tmp <- sprintf("%s/INTERACTIONS/%s", trainNetDir, tmp)
        if (sum(grepl(pattern = ",", readLines(tmp, n = 1)) > 
            0)) {
            replacePattern(path = sprintf("%s/INTERACTIONS", 
                trainNetDir))
        }
        x <- compileFeatures(trainNetDir, newOut, verbose = FALSE)
        trainPred <- pheno_train$ID[which(pheno_train$STATUS %in% 
            predClass)]
        resDir <- paste(newOut, "GM_results", sep = "/")
        dbPath <- paste(newOut, "dataset", sep = "/")
        t0 <- Sys.time()
        runFeatureSelection(trainID_pred = trainPred, outDir = resDir, 
            dbPath = dbPath, numTrainSamps = nrow(p_train), verbose = FALSE, 
            numCores = FS_numCores, featScoreMax = featScoreMax, 
            ...)
        t1 <- Sys.time()
        message("Score features for this train/test split")
        print(t1 - t0)
        tmp <- dir(path = resDir, pattern = "RANK$")[1]
        tmp <- sprintf("%s/%s", resDir, tmp)
        if (sum(grepl(pattern = ",", readLines(tmp, n = 6)) > 
            0)) {
            replacePattern(path = resDir, fileType = "RANK$")
        }
        nrankFiles <- paste(resDir, dir(path = resDir, pattern = "NRANK$"), 
            sep = "/")
        pathwayRank <- compileFeatureScores(nrankFiles, filter_WtSum = filter_WtSum, 
            verbose = FALSE)
        write.table(pathwayRank, file = paste(resDir, "pathwayScore.txt", 
            sep = "/"), col.names = TRUE, row.names = FALSE, 
            quote = FALSE)
        pScore[[k]] <- pathwayRank
    }
    phenoDF$TT_STATUS <- NULL
    out2 <- RR_featureTally(netmat, phenoDF, TT_STATUS, predClass, 
        pScore, outDir, enrichLabels, enrichedNets, verbose = FALSE, 
        maxScore = numSplits * featScoreMax)
    colnames(netmat) <- sub("_cont.txt", "", colnames(netmat))
    out <- list(netmat = netmat, pheno = phenoDF, TT_STATUS = TT_STATUS, 
        pathwayScores = pScore, enrichedNets = enrichedNets)
    for (k in names(out2)) {
        out[[k]] <- out2[[k]]
    }
    return(out)
}
callFeatSel <-
function (netScores, fsCutoff, fsPctPass) 
{
    fs_nets <- c()
    for (index in seq_len(nrow(netScores))) {
        cur_pathway <- netScores[index, ]
        pass_thresh <- length(which(cur_pathway >= fsCutoff))
        percent_pass <- pass_thresh/length(cur_pathway)
        if (percent_pass >= fsPctPass) {
            fs_nets <- c(fs_nets, netScores[, 1][index])
        }
    }
    return(fs_nets)
}
callOverallSelectedFeatures <-
function (featScores, featureSelCutoff, featureSelPct, cleanNames = TRUE) 
{
    featScores2 <- lapply(featScores, getNetConsensus)
    if (cleanNames) {
        featScores2 <- lapply(featScores2, function(x) {
            x$PATHWAY_NAME <- sub(".profile", "", x$PATHWAY_NAME)
            x$PATHWAY_NAME <- sub("_cont.txt", "", x$PATHWAY_NAME)
            colnames(x)[1] <- "Feature"
            x
        })
    }
    featSelNet <- lapply(featScores2, function(x) {
        x <- callFeatSel(x, fsCutoff = featureSelCutoff, fsPctPass = featureSelPct)
    })
    return(list(featScores = featScores2, selectedFeatures = featSelNet))
}
cleanPathwayName <-
function (curP) 
{
    pforfile <- gsub(" ", "_", curP)
    pforfile <- gsub("<", "_", pforfile)
    pforfile <- gsub(">", "_", pforfile)
    pforfile <- gsub("\\(", "_", pforfile)
    pforfile <- gsub("\\)", "_", pforfile)
    pforfile <- gsub("&", "_", pforfile)
    pforfile <- gsub(";", "_", pforfile)
    pforfile <- gsub(":", "_", pforfile)
    pforfile <- gsub("\\/", "_", pforfile)
    pforfile <- gsub("\\\xec", "X", pforfile)
    pforfile <- gsub("\\\xc2\\\xa0", "_", pforfile)
    pforfile <- gsub("\\\xa0", "X", pforfile)
    pforfile <- gsub("\\\xca", "_", pforfile)
    pforfile <- gsub("\\+", "plus", pforfile)
    pforfile <- gsub(",", ".", pforfile)
    return(pforfile)
}
compareShortestPath <-
function (net, pheno, plotDist = FALSE, verbose = TRUE) 
{
    colnames(net) <- c("source", "target", "weight")
    if (verbose) {
        message("Weight distribution:")
        print(summary(net[, 3]))
    }
    .getAvgD <- function(mat) {
        tmp <- mat[upper.tri(mat, diag = FALSE)]
        idx <- which(is.infinite(tmp))
        if (any(idx)) 
            tmp <- tmp[-idx]
        c(mean(tmp, na.rm = TRUE), sd(tmp, na.rm = TRUE), length(tmp))
    }
    .getAllD <- function(mat) {
        tmp <- mat[upper.tri(mat, diag = FALSE)]
        idx <- which(is.infinite(tmp))
        if (any(idx)) 
            tmp <- tmp[-idx]
        tmp
    }
    g <- igraph::graph_from_data_frame(net, vertices = pheno$ID)
    d_overall <- igraph::shortest.paths(g, algorithm = "dijkstra")
    if (verbose) 
        message(sprintf("Overall: %i nodes", length(pheno$ID)))
    tmp <- .getAvgD(d_overall)
    if (verbose) {
        message(sprintf("All-all shortest path = %2.3f (SD=%2.3f)", 
            tmp[1], tmp[2]))
        message(sprintf("(N=%i distances)", tmp[3]))
    }
    cnames <- unique(pheno$GROUP)
    dset <- list()
    dall <- list()
    for (curr_cl in cnames) {
        cl <- pheno$ID[which(pheno$GROUP %in% curr_cl)]
        if (verbose) {
            message(sprintf("\n%s: N=%i nodes", curr_cl, length(cl)))
        }
        tmp <- net[which(net[, 1] %in% cl & net[, 2] %in% cl), 
            ]
        g2 <- igraph::graph_from_data_frame(d = tmp, vertices = cl)
        tmp <- igraph::shortest.paths(g2, algorithm = "dijkstra")
        if (verbose) 
            message(sprintf("%s", curr_cl))
        dset[[curr_cl]] <- .getAvgD(tmp)
        dall[[curr_cl]] <- .getAllD(tmp)
        tmp <- dset[[curr_cl]]
        if (verbose) {
            message(sprintf("\t%s-%s: Mean shortest = %2.3f (SD= %2.3f)", 
                curr_cl, curr_cl, tmp[1], tmp[2]))
            message("(N=%i dist)", tmp[3])
        }
    }
    cpairs <- as.matrix(combinat::combn(cnames, 2))
    if (verbose) 
        message("Pairwise classes:")
    for (k in seq_len(ncol(cpairs))) {
        type1 <- pheno$ID[which(pheno$GROUP %in% cpairs[1, k])]
        type2 <- pheno$ID[which(pheno$GROUP %in% cpairs[2, k])]
        idx <- which(net[, 1] %in% type1 & net[, 2] %in% type2)
        idx2 <- which(net[, 1] %in% type2 & net[, 2] %in% type1)
        idx <- c(idx, idx2)
        g <- igraph::graph_from_data_frame(d = net[idx, ])
        tmp <- igraph::shortest.paths(g, algorithm = "dijkstra")
        cur <- sprintf("%s-%s", cpairs[1, k], cpairs[2, k])
        dset[[cur]] <- .getAvgD(tmp)
        dall[[cur]] <- .getAllD(tmp)
        tmp <- dset[[curr_cl]]
        if (verbose) {
            message(sprintf("\t%s-%s: Mean shortest = %2.3f (SD= %2.3f)", 
                cpairs[1, k], cpairs[2, k], tmp[1], tmp[2]))
            message(sprintf("(N=%i dist)", tmp[3]))
        }
    }
    dset[["overall"]] <- .getAvgD(d_overall)
    dall[["overall"]] <- .getAllD(d_overall)
    out <- list(avg = dset, all = dall)
    if (plotDist) {
        par(las = 1, bty = "n")
        dl <- data.frame(intType = rep(names(dall), lapply(dall, 
            length)), dijk = unlist(dall))
        plotList <- list()
        p <- ggplot(dl, aes(dl[, "intType"], dl[, "dijk"]))
        p <- p + ylab("Pairwise Dijkstra distance\n(smaller is better)")
        p <- p + xlab("Pair groups")
        p2 <- p + geom_violin(scale = "width") + geom_boxplot(width = 0.02)
        print(p2)
        out[["plot"]] <- p2
    }
    return(out)
}
compileFeatures <-
function (netDir, outDir = tempdir(), simMetric = "pearson", 
    netSfx = "txt$", verbose = TRUE, numCores = 1L, P2N_threshType = "off", 
    P2N_maxMissing = 100, JavaMemory = 4L, altBaseDir = NULL, 
    debugMode = FALSE, ...) 
{
    dataDir <- paste(outDir, "dataset", sep = "/")
    GM_jar <- getGMjar_path()
    if (P2N_maxMissing < 5) 
        PSN_maxMissing <- 5
    if (P2N_maxMissing > 100) 
        PSN_maxMissing <- 100
    if (!P2N_threshType %in% c("off", "auto")) 
        P2N_threshType <- "off"
    if (!file.exists(dataDir)) 
        dir.create(dataDir)
    curwd <- getwd()
    setwd(netDir)
    netList1 <- dir(path = paste(netDir, "profiles", sep = "/"), 
        pattern = "profile$")
    netList2 <- dir(path = paste(netDir, "INTERACTIONS", sep = "/"), 
        pattern = netSfx)
    netList <- c(netList1, netList2)
    if (verbose) 
        message(sprintf("Got %i networks", length(netList)))
    idFile <- paste(outDir, "ids.txt", sep = "/")
    writeQueryBatchFile(netDir, netList, netDir, idFile, ...)
    if (length(netList1) > 0) {
        if (verbose) 
            message("\t* Converting profiles to interaction networks")
        cl <- makeCluster(numCores, outfile = paste(netDir, "P2N_log.txt", 
            sep = "/"))
        registerDoParallel(cl)
        if (simMetric == "pearson") {
            corType <- "PEARSON"
        }
        else if (simMetric == "MI") {
            corType <- "MUTUAL_INFORMATION"
        }
        args <- c(sprintf("-Xmx%iG", JavaMemory), "-cp", GM_jar)
        args <- c(args, paste("org.genemania.engine.core.", "evaluation.ProfileToNetworkDriver", 
            sep = ""))
        args <- c(args, c("-proftype", "continuous", "-cor", 
            corType))
        args <- c(args, c("-threshold", P2N_threshType, "-maxmissing", 
            sprintf("%1.1f", P2N_maxMissing)))
        profDir <- paste(netDir, "profiles", sep = "/")
        netOutDir <- paste(netDir, "INTERACTIONS", sep = "/")
        tmpsfx <- sub("\\$", "", netSfx)
        curProf <- ""
        `%myinfix%` <- ifelse(debugMode, `%do%`, `%dopar%`)
        foreach(curProf = dir(path = profDir, pattern = "profile$")) %myinfix% 
            {
                args2 <- c("-in", paste(profDir, curProf, sep = "/"))
                args2 <- c(args2, "-out", paste(netOutDir, sub(".profile", 
                  ".txt", curProf), sep = "/"))
                args2 <- c(args2, "-syn", paste(netDir, "1.synonyms", 
                  sep = "/"), "-keepAllTies", "-limitTies")
                if (debugMode) {
                  message("Making Java call")
                  tmp <- paste(c(args, args2), collapse = " ")
                  message(sprintf("java %s", tmp))
                  system2("java", args = c(args, args2), wait = TRUE)
                }
                else {
                  system2("java", args = c(args, args2), wait = TRUE, 
                    stdout = NULL)
                }
            }
        stopCluster(cl)
        netSfx = ".txt"
        netList2 <- dir(path = netOutDir, pattern = netSfx)
        msg2 <- paste("This problem usually occurs because of a failed", 
            "Java call. Try upgrading to Java 11. If that doesn't", 
            "work, contact shraddha.pai@utoronto.ca with a copy", 
            "of the complete log file after running buildPredict()", 
            "with debugMode=TRUE", sep = "\n")
        if (length(netList2) < length(netList)) {
            browser()
            warnings(paste("", "---------------------------------", 
                "One or more profiles did not successfully convert to PSNs!", 
                "This usually happens because of the underlying call to a", 
                "Java library failed. Upgrading to Java 11 usually fixes", 
                "this problem. If not, please send a copy of the detailed", 
                "error message from the call below to", "shraddha.pai@utoronto.ca", 
                "", sep = "\n"))
            curProf <- dir(profDir, "profile$")[1]
            args2 <- c("-in", paste(profDir, curProf, sep = "/"))
            args2 <- c(args2, "-out", paste(netOutDir, sub(".profile", 
                ".txt", curProf), sep = "/"))
            args2 <- c(args2, "-syn", paste(netDir, "1.synonyms", 
                sep = "/"), "-keepAllTies", "-limitTies")
            tmp <- paste(c(args, args2), collapse = " ")
            print(sprintf("java %s", tmp))
            system2("java", args = c(args, args2), wait = TRUE)
            stop("Stopping netDx now. See error message above.")
        }
        if (verbose) 
            message(sprintf("Got %i networks from %i profiles", 
                length(netList2), length(netList)))
        netList <- netList2
        rm(netOutDir, netList2)
    }
    if (verbose) 
        message("\t* Build GeneMANIA index")
    setwd(dataDir)
    args <- c("-Xmx10G", "-cp", GM_jar)
    args <- c(args, paste("org.genemania.mediator.lucene.", "exporter.Generic2LuceneExporter", 
        sep = ""))
    args <- c(args, paste(netDir, "db.cfg", sep = "/"), 
        netDir, paste(netDir, "colours.txt", sep = "/"))
    if (debugMode) {
        tmp <- paste(args, collapse = " ")
        message(sprintf("java %s", tmp))
        system2("java", args, wait = TRUE)
    }
    else {
        system2("java", args, wait = TRUE, stdout = NULL)
    }
    olddir <- paste(dataDir, "lucene_index", sep = "/")
    flist <- list.files(olddir, recursive = TRUE)
    dirs <- list.dirs(olddir, recursive = TRUE, full.names = FALSE)
    dirs <- setdiff(dirs, "")
    for (d in dirs) dir.create(paste(dataDir, d, sep = "/"))
    file.copy(from = paste(olddir, flist, sep = "/"), 
        to = paste(dataDir, flist, sep = "/"))
    unlink(olddir)
    tmp <- dir(path = sprintf("%s/INTERACTIONS", netDir), pattern = "txt$")[1]
    tmp <- sprintf("%s/INTERACTIONS/%s", netDir, tmp)
    if (sum(grepl(pattern = ",", readLines(tmp, n = 1)) > 0)) {
        replacePattern(path = sprintf("%s/INTERACTIONS", netDir))
    }
    if (verbose) 
        message("\t* Build GeneMANIA cache")
    args <- c("-Xmx10G", "-cp", GM_jar, "org.genemania.engine.apps.CacheBuilder")
    args <- c(args, "-cachedir", "cache", "-indexDir", ".", "-networkDir", 
        paste(netDir, "INTERACTIONS", sep = "/"), "-log", 
        paste(netDir, "test.log", sep = "/"))
    if (debugMode) {
        tmp <- paste(args, collapse = " ")
        message(sprintf("java %s", tmp))
        system2("java", args = args)
    }
    else {
        system2("java", args = args, stdout = NULL)
    }
    if (verbose) 
        message("\t * Cleanup")
    GM_xml <- system.file("extdata", "genemania.xml", package = "netDx")
    file.copy(from = GM_xml, to = paste(dataDir, ".", sep = "/"))
    setwd(curwd)
    return(list(dbDir = dataDir, netDir = netDir))
}
compileFeatureScores <-
function (fList, filter_WtSum = 100, verbose = FALSE) 
{
    if (filter_WtSum < 5) {
        message("filter_WtSum cannot be < 5 ; setting to 5")
        filter_WtSum <- 5
    }
    pathwayTally <- list()
    ctr <- 1
    for (fName in fList) {
        tmp <- basename(fName)
        try(dat <- read.delim(fName, sep = "\t", header = TRUE, 
            as.is = TRUE, skip = 1), silent = TRUE)
        ctr <- ctr + 1
        if (!inherits(dat, "try-error")) {
            if (verbose) {
                message("Net weight distribution:")
                print(summary(dat$Weight))
            }
            dat <- dat[order(dat$Weight, decreasing = TRUE), 
                ]
            cs <- cumsum(dat$Weight)
            keep_max <- which.min(abs(cs - filter_WtSum))
            dat <- dat[seq_len(keep_max), ]
            if (verbose) {
                message(sprintf(paste("filter_WtSum = %1.1f; ", 
                  "%i of %i networks left", sep = ""), filter_WtSum, 
                  nrow(dat), length(cs)))
            }
            for (k in dat$Network) {
                if (!k %in% names(pathwayTally)) 
                  pathwayTally[[k]] <- 0
                pathwayTally[[k]] <- pathwayTally[[k]] + 1
            }
        }
    }
    out <- unlist(pathwayTally)
    out <- sort(out, decreasing = TRUE)
    out <- data.frame(name = names(out), score = as.integer(out), 
        stringsAsFactors = FALSE)
    out[, 2] <- as.integer(as.character(out[, 2]))
    out
}
confusionMatrix <-
function (model) 
{
    nmList <- names(model)[grep("Split", names(model))]
    cl <- sort(unique(model$Split1$STATUS))
    conf <- list()
    mega <- NULL
    for (nm in nmList) {
        pred <- model[[nm]][["predictions"]][, c("ID", "STATUS", 
            "TT_STATUS", "PRED_CLASS")]
        m <- as.matrix(table(pred[, c("STATUS", "PRED_CLASS")]))
        conf[[nm]] <- m/colSums(m)
        if (is.null(mega)) 
            mega <- conf[[nm]]
        else mega <- mega + conf[[nm]]
    }
    mega <- mega/length(conf)
    mega <- round(mega * 100, 2)
    mega <- t(mega)
    metric <- "%% Accuracy"
    tbl <- table(model$Split1$predictions$STATUS)
    nm <- names(tbl)
    val <- as.integer(tbl)
    ttl <- sprintf("%s\n(N=%i)", rownames(mega), val[match(rownames(mega), 
        nm)])
    par(mar = c(4, 8, 2, 2))
    color2D.matplot(mega, show.values = TRUE, border = "white", 
        extremes = c(1, 2), axes = FALSE, xlab = "Predicted class", 
        ylab = "")
    axis(1, at = seq_len(ncol(mega)) - 0.5, labels = colnames(mega))
    axis(2, at = seq_len(ncol(mega)) - 0.5, labels = rev(ttl), 
        las = 2)
    title(sprintf("Confusion matrix: Accuracy (avg of %i splits)", 
        length(conf)))
    return(list(splitWiseConfMatrix = conf, average = mega))
}
convertProfileToNetworks <-
function (netDir, outDir = tempdir(), simMetric = "pearson", 
    numCores = 1L, JavaMemory = 4L, GM_jar = NULL, P2N_threshType = "off", 
    P2N_maxMissing = 100, netSfx = "txt$", debugMode = FALSE) 
{
    if (is.null(GM_jar)) 
        GM_jar <- getGMjar_path()
    cl <- makeCluster(numCores, outfile = paste(netDir, "P2N_log.txt", 
        sep = "/"))
    registerDoParallel(cl)
    if (simMetric == "pearson") {
        corType <- "PEARSON"
    }
    else if (simMetric == "MI") {
        corType <- "MUTUAL_INFORMATION"
    }
    args <- c(sprintf("-Xmx%iG", JavaMemory), "-cp", GM_jar)
    args <- c(args, paste("org.genemania.engine.core.", "evaluation.ProfileToNetworkDriver", 
        sep = ""))
    args <- c(args, c("-proftype", "continuous", "-cor", corType))
    args <- c(args, c("-threshold", P2N_threshType, "-maxmissing", 
        sprintf("%1.1f", P2N_maxMissing)))
    profDir <- netDir
    tmpsfx <- sub("\\$", "", netSfx)
    curProf <- ""
    `%myinfix%` <- ifelse(debugMode, `%do%`, `%dopar%`)
    foreach(curProf = dir(path = profDir, pattern = "profile$")) %myinfix% 
        {
            if (debugMode) 
                print(curProf)
            args2 <- c("-in", paste(profDir, curProf, sep = "/"))
            args2 <- c(args2, "-out", paste(outDir, sub(".profile", 
                ".txt", curProf), sep = "/"))
            args2 <- c(args2, "-syn", paste(netDir, "..", "1.synonyms", 
                sep = "/"), "-keepAllTies", "-limitTies")
            if (debugMode) 
                stdout <- ""
            else stdout <- NULL
            system2("java", args = c(args, args2), wait = TRUE, 
                stdout = stdout)
        }
    tmp <- dir(path = outDir, pattern = "txt$")[1]
    tmp <- sprintf("%s/%s", outDir, tmp)
    if (sum(grepl(pattern = ",", readLines(tmp, n = 1)) > 0)) {
        replacePattern(path = outDir, fileType = "txt$")
    }
    stopCluster(cl)
}
convertToMAE <-
function (dataList) 
{
    if (class(dataList) != "list") {
        stop("dataList must be a list. \n")
    }
    if (is.null(dataList$pheno)) {
        stop("dataList must have key-value pair labelled pheno.\n")
    }
    if (length(dataList) == 1) {
        stop("dataList must have assay data to incorporate into a \n         MultiAssayExperiment object")
    }
    patientPheno <- dataList$pheno
    tmp <- NULL
    track <- c()
    datType <- names(dataList)
    for (k in 1:length(dataList)) {
        if (names(dataList[k]) != "pheno") {
            if (sum(duplicated(colnames(dataList[[k]]))) != 0) {
                dataList[[k]] <- dataList[[k]][, !duplicated(colnames(dataList[[k]]))]
            }
            track <- c(track, k)
            tmp <- c(tmp, list(dataList[[k]]))
        }
    }
    names(tmp) <- datType[track]
    MAE <- MultiAssayExperiment(experiments = tmp, colData = patientPheno)
    return(MAE)
}
countIntType <-
function (inFile, plusID, minusID) 
{
    dat <- read.delim(inFile, sep = "\t", header = FALSE, as.is = TRUE)
    pp <- sum(dat[, 1] %in% plusID & dat[, 2] %in% plusID)
    return(c(pp, nrow(dat) - pp))
}
countIntType_batch <-
function (inFiles, plusID, minusID, tmpDir = tempdir(), enrType = "binary", 
    numCores = 1L) 
{
    randString <- function(n = 1) {
        a <- do.call(paste0, replicate(5, sample(LETTERS, n, 
            TRUE), FALSE))
        paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, 
            n, TRUE))
    }
    curRand <- randString()
    bkFile <- sprintf("%s.bk", curRand)
    descFile <- sprintf("%s.desc", curRand)
    bkFull <- paste(tmpDir, bkFile, sep = "/")
    descFull <- paste(tmpDir, descFile, sep = "/")
    if (file.exists(bkFull)) 
        file.remove(bkFull)
    if (file.exists(descFull)) 
        file.remove(descFull)
    out <- big.matrix(NA, nrow = length(inFiles), ncol = 2, type = "double", 
        backingfile = bkFile, backingpath = tmpDir, descriptorfile = descFile)
    cl <- makeCluster(numCores, outfile = paste(tmpDir, "shuffled_log.txt", 
        sep = "/"))
    registerDoParallel(cl)
    k <- 0
    foreach(k = seq_len(length(inFiles))) %dopar% {
        m <- attach.big.matrix(paste(tmpDir, descFile, sep = "/"))
        if (enrType == "binary") 
            m[k, ] <- countIntType(inFiles[k], plusID, minusID)
        else if (enrType == "corr") 
            m[k, ] <- getCorrType(inFiles[k], plusID, minusID)
    }
    stopCluster(cl)
    out <- as.matrix(out)
    unlink(bkFull)
    unlink(descFull)
    return(out)
}
countPatientsInNet <-
function (netDir, fList, ids) 
{
    outmat <- matrix(0, nrow = length(ids), ncol = length(fList))
    colnames(outmat) <- fList
    rownames(outmat) <- ids
    ctr <- 1
    for (f in fList) {
        dat <- read.delim(paste(netDir, f, sep = "/"), 
            sep = "\t", header = FALSE, as.is = TRUE)
        memb <- c(dat[, 1], dat[, 2])
        outmat[which(ids %in% memb), ctr] <- 1
        ctr <- ctr + 1
    }
    return(outmat)
}
createNetFuncFromSimList <-
function (dataList, groupList, netDir, sims, verbose = TRUE, 
    ...) 
{
    if (length(groupList) != length(sims)) {
        stop("groupList and sims need to be of same length.")
    }
    if (all.equal(sort(names(groupList)), sort(names(sims))) != 
        TRUE) {
        stop("names(groupList) needs to match names(sims).")
    }
    settings <- list(dataList = dataList, groupList = groupList, 
        netDir = netDir, sims = sims)
    if (verbose) 
        message("Making nets from sims")
    netList <- c()
    for (nm in names(sims)) {
        csim <- sims[[nm]]
        netList_cur <- NULL
        if (verbose) 
            message(sprintf("\t%s", nm))
        cur_set <- settings
        cur_set[["name"]] <- nm
        cur_set[["similarity"]] <- csim
        if (!is.null(groupList[[nm]])) {
            if (class(csim) == "function") {
                netList_cur <- psn__custom(cur_set, csim, verbose, 
                  ...)
            }
            else if (csim == "pearsonCorr") {
                netList_cur <- psn__corr(cur_set, verbose, ...)
            }
            else {
                netList_cur <- psn__builtIn(cur_set, verbose, 
                  ...)
            }
            netList <- c(netList, netList_cur)
        }
    }
    if (verbose) {
        message("Net construction complete!")
    }
    unlist(netList)
}
createPSN_MultiData <-
function (dataList, groupList, pheno, netDir = tempdir(), filterSet = NULL, 
    verbose = TRUE, makeNetFunc = NULL, sims = NULL, ...) 
{
    if (missing(dataList)) 
        stop("dataList must be supplied.\n")
    if (missing(groupList)) 
        stop("groupList must be supplied.\n")
    dataList <- lapply(dataList, function(x) {
        midx <- match(colnames(x), pheno$ID)
        colnames(x) <- pheno$INTERNAL_ID[midx]
        x
    })
    if (!is.null(filterSet)) {
        if (length(filterSet) < 1) {
            s1 <- "filterSet is empty."
            s2 <- "It needs to have at least one net to proceed."
            stop(paste(s1, s2, sep = " "))
        }
    }
    if (!is.null(filterSet)) {
        if (verbose) 
            message("\tFilter set provided")
        groupList2 <- list()
        for (nm in names(groupList)) {
            idx <- which(names(groupList[[nm]]) %in% filterSet)
            if (verbose) {
                message(sprintf("\t\t%s: %i of %i nets pass", 
                  nm, length(idx), length(groupList[[nm]])))
            }
            if (length(idx) > 0) {
                groupList2[[nm]] <- groupList[[nm]][idx]
            }
        }
        groupList <- groupList2
        sims <- sims[which(names(sims) %in% names(groupList))]
        rm(groupList2)
    }
    if (!is.null(makeNetFunc)) {
        netList <- makeNetFunc(dataList = dataList, groupList = groupList, 
            netDir = netDir, ...)
    }
    else {
        netList <- createNetFuncFromSimList(dataList = dataList, 
            groupList = groupList, netDir = netDir, sims = sims, 
            ...)
    }
    if (length(netList) < 1) 
        stop("\n\nNo features created! Filters may be too stringent.\n")
    netID <- data.frame(ID = seq_len(length(netList)), name = netList, 
        ID = seq_len(length(netList)), name2 = netList, 0, 1, 
        stringsAsFactors = TRUE)
    fsep = "/"
    prof <- grep(".profile$", netList)
    if (length(prof) > 0) {
        prof <- netList[prof]
        dir.create(paste(netDir, "profiles", sep = fsep))
        for (p in prof) {
            file.rename(from = paste(netDir, p, sep = fsep), 
                to = paste(netDir, "profiles", sprintf("1.%i.profile", 
                  netID$ID[which(netID$name == p)]), sep = fsep))
        }
    }
    dir.create(paste(netDir, "INTERACTIONS", sep = fsep))
    cont <- grep("_cont.txt$", netList)
    if (length(cont) > 0) {
        cont <- netList[cont]
        for (p in cont) {
            file.rename(from = paste(netDir, p, sep = fsep), 
                to = paste(netDir, "INTERACTIONS", sprintf("1.%i.txt", 
                  netID$ID[which(netID$name == p)]), sep = fsep))
        }
    }
    write.table(netID, file = paste(netDir, "NETWORKS.txt", sep = fsep), 
        sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    con <- file(paste(netDir, "NETWORK_GROUPS.txt", sep = fsep), 
        "w")
    write(paste(1, "dummy_group", "geneset_1", "dummy_group", 
        1, sep = "\t"), file = con)
    close(con)
    con <- file(paste(netDir, "NETWORK_METADATA.txt", sep = fsep), 
        "w")
    tmp <- paste(netID$ID, "", "", "", "", "", "", "", "", "", 
        0, "", "", 0, "", "", "", "", "", sep = "\t")
    write.table(tmp, file = con, sep = "\t", col.names = FALSE, 
        row.names = FALSE, quote = FALSE)
    close(con)
    return(netList)
}
dataList2List <-
function (dat, groupList) 
{
    exprs <- experiments(dat)
    datList2 <- list()
    for (k in seq_len(length(exprs))) {
        tmp <- exprs[[k]]
        df <- sampleMap(dat)[which(sampleMap(dat)$assay == names(exprs)[k]), 
            ]
        colnames(tmp) <- df$primary[match(df$colname, colnames(tmp))]
        if ("SimpleList" %in% class(tmp)) {
            tmp <- as.matrix(assays(tmp)[[1]])
        }
        else if ("SummarizedExperiment" %in% class(tmp)) {
            tmp <- as.matrix(assays(tmp)[[1]])
        }
        datList2[[names(exprs)[k]]] <- tmp
    }
    if ("clinical" %in% names(groupList)) {
        tmp <- colData(dat)
        vars <- unique(unlist(groupList[["clinical"]]))
        datList2[["clinical"]] <- t(as.matrix(tmp[, vars, drop = FALSE]))
    }
    pheno_all <- colData(dat)
    pheno_all <- as.data.frame(pheno_all)
    out <- list(assays = datList2, pheno = pheno_all)
}
enrichLabelNets <-
function (netDir, pheno_DF, outDir, numReps = 50L, minEnr = -1, 
    outPref = "enrichLabelNets", verbose = TRUE, setSeed = 42L, 
    enrType = "binary", numCores = 1L, predClass, tmpDir = tempdir(), 
    netGrep = "_cont.txt$", getShufResults = FALSE, ...) 
{
    today <- format(Sys.Date(), "%y%m%d")
    runtime <- format(Sys.time(), "%H%M")
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    orig_enr <- getEnr(netDir, pheno_DF, predClass = predClass, 
        tmpDir = tmpDir, enrType = enrType, netGrep = netGrep, 
        ...)
    plusID <- orig_enr[["plusID"]]
    minusID <- orig_enr[["minusID"]]
    orig <- orig_enr[["orig"]]
    orig_rat <- orig_enr[["orig_rat"]]
    fList <- orig_enr[["fList"]]
    print(summary(orig_rat))
    idx <- which(orig_rat >= minEnr)
    message(sprintf("\n\t%i of %i networks have ENR >= %1.1f -> filter", 
        length(idx), length(fList), minEnr))
    fList <- fList[idx]
    orig_rat <- orig_rat[idx]
    orig <- orig[idx, ]
    message("* Computing shuffled")
    both <- c(plusID, minusID)
    n1 <- length(plusID)
    n <- length(both)
    drop_out <- integer()
    N <- length(fList)
    to_run <- seq_len(N)
    currRep <- 1
    shuf_rat <- matrix(NA, nrow = length(fList), ncol = numReps)
    x0 <- system.time(while ((length(to_run) > 0) & (currRep <= 
        numReps)) {
        shuf <- sample(both, replace = FALSE)
        tmp <- countIntType_batch(fList[to_run], shuf[seq_len(n1)], 
            shuf[(n1 + 1):n], tmpDir = tmpDir, enrType = enrType, 
            numCores = numCores, ...)
        if (length(to_run) < 2) 
            tmp <- matrix(tmp, ncol = 2)
        if (enrType == "binary") {
            shuf_rat[to_run, currRep] <- (tmp[, 1] - tmp[, 2])/(tmp[, 
                1] + tmp[, 2])
        }
        else if (enrType == "corr") {
            shuf_rat[to_run, currRep] <- (tmp[, 1] - tmp[, 2])/2
        }
        else {
            shuf_rat[to_run, currRep] <- NA
        }
        if (currRep%%10 == 0) {
            orig_pct <- numeric(length(to_run))
            ctr <- 1
            cols <- seq_len(currRep)
            for (k in to_run) {
                num <- sum(shuf_rat[k, cols] >= orig_rat[k])
                orig_pct[ctr] <- num/currRep
                ctr <- ctr + 1
            }
            y <- exp(-(currRep - 1)/100)
            y2 <- 3/currRep
            if (y - y2 > 0) 
                chk_thresh <- y
            else chk_thresh <- y2
            idx <- which(orig_pct >= chk_thresh)
            if (length(idx) > 0) 
                drop_out <- c(drop_out, to_run[idx])
            to_run <- setdiff(seq_len(N), drop_out)
            if (verbose) {
                message(sprintf("%i (%1.5f) %i drop out; %i left", 
                  currRep, chk_thresh, length(idx), length(to_run)))
            }
            else {
                message(".")
            }
        }
        if (currRep%%100 == 0 & !verbose) 
            message("\t")
        currRep <- currRep + 1
    })
    stopCluster(cl)
    mu <- rowMeans(shuf_rat, na.rm = TRUE)
    sigma <- apply(shuf_rat, 1, sd, na.rm = TRUE)
    orig_z <- (orig_rat - mu)/sigma
    orig_pct <- numeric(length(orig_rat))
    for (i in seq_len(length(orig_rat))) {
        tmp <- na.omit(shuf_rat[i, ])
        orig_pct[i] <- sum(tmp >= orig_rat[i])/length(tmp)
    }
    maxShufs <- numeric(nrow(shuf_rat))
    for (i in seq_len(nrow(shuf_rat))) {
        maxShufs[i] <- max(which(!is.na(shuf_rat[i, ])))
    }
    qval <- p.adjust(orig_pct, method = "BH")
    out <- data.frame(NETWORK = basename(fList), orig_pp = orig[, 
        1], orig_rest = orig[, 2], ENR = orig_rat, TOTAL_INT = log10(orig[, 
        1] + orig[, 2]), numPerm = maxShufs, shuf_mu = mu, shuf_sigma = sigma, 
        Z = orig_z, pctl = orig_pct, Q = qval)
    write.table(out, file = paste(outDir, sprintf("%s.stats.txt", 
        outPref), sep = "/"), sep = "\t", col.names = TRUE, 
        row.names = FALSE, quote = FALSE)
    if (getShufResults) {
        out <- list(shufres = shuf_rat, res = out)
    }
    else {
        return(out)
    }
    message("about to leave enrich nets")
}
fetchPathwayDefinitions <-
function (month = NULL, year = NULL, day = 1, verbose = FALSE) 
{
    if (is.null(month) || is.null(year)) {
        stop("Please provide a month and year.")
    }
    if (year < 2020) {
        stop("Currently, year must equal 2020 or greater.")
    }
    if (class(month) %in% c("numeric", "integer")) {
        month <- month.name[month]
    }
    pdate <- sprintf("%s_%02d_%i", month, day, year)
    pathwayURL <- paste("https://downloads.res.oicr.on.ca/pailab/public/EM_Genesets/", 
        sprintf("%s/Human/symbol/", pdate), sprintf("Human_AllPathways_%s_symbol.gmt", 
            pdate), sep = "")
    message(sprintf("Fetching %s", pathwayURL))
    bfc <- .get_cache()
    chk <- httr::HEAD(pathwayURL)
    if (chk$status_code == 404) {
        stop(paste(sprintf("The pathway file for %02d %s %i doesn't exist.", 
            day, month, year), "Select a different date. ", "See https://downloads.res.oicr.on.ca/pailab/public/EM_Genesets/Human/symbol for options.", 
            sep = " "))
    }
    bfcrpath(bfc, pathwayURL)
}
getEMapInput <-
function (featScores, namedSets, netInfo, pctPass = 0.69999999999999996, 
    minScore = 1, maxScore = 10, trimFromName = c(".profile", 
        "_cont"), verbose = FALSE) 
{
    netNames <- featScores[, 1]
    featScores <- as.matrix(featScores[, -1])
    maxNetS <- matrix(NA, nrow = length(netNames), ncol = 1)
    for (sc in minScore:maxScore) {
        if (ncol(featScores) >= 2) {
            tmp <- rowSums(featScores >= sc)
            idx <- which(tmp >= floor(pctPass * ncol(featScores)))
        }
        else {
            idx <- which(featScores >= sc)
        }
        if (verbose) 
            message(sprintf("\t%i : %i pass", sc, length(idx)))
        maxNetS[idx, 1] <- sc
    }
    idx <- which(!is.na(maxNetS))
    maxNetS <- maxNetS[idx, , drop = FALSE]
    netNames <- netNames[idx]
    for (tr in trimFromName) netNames <- sub(tr, "", netNames)
    df1 <- data.frame(netName = netNames, maxScore = maxNetS)
    colnames(netInfo) <- c("netType", "netName")
    df2 <- merge(x = df1, y = netInfo, by = "netName")
    featSet <- list()
    for (cur in df2$netName) {
        k2 <- simpleCap(cur)
        if (is.null(namedSets[[cur]])) 
            namedSets[[cur]] <- k2
        featSet[[k2]] <- namedSets[[cur]]
    }
    out <- list(nodeAttrs = df2, featureSets = featSet)
    return(out)
}
getEMapInput_many <-
function (featScores, namedSets_valid, netTypes, outDir, ...) 
{
    out <- list()
    for (gp in names(featScores)) {
        cur_out_files <- getEMapInput(featScores[[gp]], namedSets_valid, 
            netTypes, ...)
        out[[gp]] <- cur_out_files
    }
    return(out)
}
getEnr <-
function (netDir, pheno_DF, predClass, netGrep = "_cont.txt$", 
    enrType = "binary", ...) 
{
    if (missing(predClass)) 
        stop("predClass must be supplied.\n")
    fList <- dir(path = netDir, pattern = netGrep)
    fList <- paste(netDir, fList, sep = "/")
    message(sprintf("Got %i networks", length(fList)))
    pheno <- pheno_DF
    pheno <- pheno[!duplicated(pheno$ID), ]
    plus_idx <- which(pheno$STATUS %in% predClass)
    plusID <- pheno$ID[plus_idx]
    minusID <- pheno$ID[setdiff(seq_len(nrow(pheno)), plus_idx)]
    message(sprintf("Total %i subjects ; %i of class %s, %i other", 
        nrow(pheno), length(plusID), predClass, length(minusID)))
    message("* Computing real (+,+) (+,-)")
    t0 <- system.time(orig <- countIntType_batch(fList, plusID, 
        minusID, enrType = enrType, ...))
    print(t0)
    if (enrType == "binary") {
        orig_rat <- (orig[, 1] - orig[, 2])/(orig[, 1] + orig[, 
            2])
    }
    else if (enrType == "corr") {
        orig_rat <- (orig[, 1] - orig[, 2])/2
    }
    else {
        orig_rat <- NA
    }
    names(orig_rat) <- basename(fList)
    return(list(plusID = plusID, minusID = minusID, orig_rat = orig_rat, 
        fList = fList, orig = orig))
}
getFeatureScores <-
function (inDir, predClasses, getFullCons = TRUE) 
{
    if (missing(inDir)) 
        stop("inDir not provided")
    if (missing(predClasses)) 
        stop("predClasses missing; please specify classes")
    out <- list()
    for (gp in predClasses) {
        message(sprintf("%s\n", gp))
        if (is(inDir, "character")) {
            message("\tSingle directory provided, retrieving CV score files\n")
            rngDirs <- dir(path = inDir, pattern = "^rng")
            fList <- paste(inDir, rngDirs, gp, "GM_results", 
                sprintf("%s_pathway_CV_score.txt", gp), sep = "/")
        }
        else {
            message("\tList of filenames provided\n")
            fList <- inDir[[gp]]
        }
        message(sprintf("Got %i iterations", length(fList)))
        netColl <- list()
        for (scoreFile in fList) {
            tmp <- read.delim(scoreFile, sep = "\t", header = TRUE, 
                as.is = TRUE)
            colnames(tmp)[1] <- "PATHWAY_NAME"
            netColl[[scoreFile]] <- tmp
        }
        spos <- gregexpr("\\/", fList)
        fNames <- lapply(seq_len(length(spos)), function(x) {
            n <- length(spos[[x]])
            y <- substr(fList[x], spos[[x]][n - 3] + 1, spos[[x]][n - 
                2] - 1)
            y
        })
        fNames <- unlist(fNames)
        names(netColl) <- fNames
        message("* Computing consensus\n")
        cons <- getNetConsensus(netColl)
        x1 <- nrow(cons)
        na_sum <- rowSums(is.na(cons))
        full_cons <- cons
        cons <- cons[which(na_sum < 1), ]
        if (getFullCons) 
            out[[gp]] <- full_cons
        else out[[gp]] <- cons
    }
    return(out)
}
getFileSep <-
function () 
{
    if (.Platform$OS.type == "windows") 
        return("\\")
    else return(.Platform$file.sep)
}
getGMjar_path <-
function (verbose = FALSE) 
{
    java_ver <- suppressWarnings(system2("java", args = "--version", 
        stdout = TRUE, stderr = NULL))
    if (any(grep(" 11", java_ver)) || any(grep(" 12", java_ver)) || 
        any(grep(" 13", java_ver)) || any(grep(" 14", java_ver)) || 
        any(grep(" 16", java_ver))) {
        if (verbose) 
            message("Java 11+ detected")
        fileURL <- paste("https://downloads.res.oicr.on.ca/pailab/netDx/java11/", 
            "genemania-netdx.jar", sep = "")
    }
    else {
        if (verbose) 
            message("Java 8 detected")
        fileURL <- paste("https://downloads.res.oicr.on.ca/pailab/netDx/java8/", 
            "genemania-netdx.jar", sep = "")
    }
    bfc <- .get_cache()
    bfcrpath(bfc, fileURL)
}
getNetConsensus <-
function (scorelist) 
{
    out <- scorelist[[1]]
    colnames(out)[2] <- names(scorelist)[1]
    for (k in 2:length(scorelist)) {
        x <- merge(x = out, y = scorelist[[k]], by = "PATHWAY_NAME", 
            all.x = TRUE, all.y = TRUE)
        colnames(x)[k + 1] <- names(scorelist)[k]
        out <- x
    }
    out
}
getOR <-
function (pNetworks, pheno_DF, predClass, netFile, verbose = TRUE) 
{
    predSamps <- pheno_DF$ID[pheno_DF$STATUS %in% predClass]
    otherSamps <- pheno_DF$ID[!pheno_DF$STATUS %in% predClass]
    idx <- which(colnames(pNetworks) %in% netFile)
    if (length(idx) < 1) 
        return(out = list(stats = matrix(NA, nrow = 2, ncol = 3), 
            relEnr = NA, OLsamps = NA))
    pNetworks <- pNetworks[, idx, drop = FALSE]
    OLsamps <- rownames(pNetworks)[which(rowSums(pNetworks) >= 
        1)]
    OLpred <- sum(OLsamps %in% predSamps)
    OLother <- sum(OLsamps %in% otherSamps)
    pctPred <- OLpred/length(predSamps)
    pctOther <- OLother/length(otherSamps)
    if (pctPred < .Machine$double.eps) 
        pctPred <- .Machine$double.eps
    if (pctOther < .Machine$double.eps) 
        pctOther <- .Machine$double.eps
    relEnr <- pctPred/pctOther
    outmat <- matrix(nrow = 2, ncol = 3)
    colnames(outmat) <- c("total", "num OL", "pct OL")
    rownames(outmat) <- c(predClass, "(other)")
    outmat[1, ] <- c(length(predSamps), OLpred, round(pctPred * 
        100, digits = 1))
    outmat[2, ] <- c(length(otherSamps), OLother, round(pctOther * 
        100, digits = 1))
    if (verbose) 
        print(outmat)
    if (verbose) 
        message(sprintf("Relative enrichment of %s: %1.3f", predClass, 
            relEnr))
    out <- list(stats = outmat, relEnr = relEnr, OLsamps = OLsamps)
    out
}
getPatientPredictions <-
function (predFiles, pheno, plotAccuracy = FALSE) 
{
    if (length(predFiles) == 1) {
        message("predFiles is of length 1. Assuming directory\n")
        all_rngs <- list.dirs(predFiles, recursive = FALSE)
        all_rngs <- all_rngs[grep("rng", all_rngs)]
        predFiles <- unlist(lapply(all_rngs, function(x) {
            paste(x, "predictionResults.txt", sep = "/")
        }))
    }
    else {
        message("predFiles is of length > 1. Assuming filenames provided\n")
    }
    output_mat <- matrix(NA, nrow = nrow(pheno), ncol = length(predFiles) + 
        2)
    patient_list <- list()
    for (cur_pat in pheno$ID) patient_list[[cur_pat]] <- c()
    uq_mat <- matrix(NA, nrow = nrow(pheno), ncol = length(predFiles))
    rownames(uq_mat) <- pheno$ID
    for (ctr in seq_len(length(predFiles))) {
        curFile <- predFiles[ctr]
        dat <- read.delim(curFile, sep = "\t", header = TRUE, 
            as.is = TRUE)
        for (k in seq_len(nrow(dat))) {
            tmp <- which(rownames(uq_mat) == dat$ID[k])
            uq_mat[tmp, ctr] <- dat$PRED_CLASS[k]
        }
    }
    uq_mat <- as.data.frame(uq_mat)
    pctCorr <- c()
    for (k in seq_len(nrow(uq_mat))) {
        testCt <- sum(!is.na(uq_mat[k, ]))
        cur <- sum(uq_mat[k, ] == pheno$STATUS[k], na.rm = TRUE)/testCt
        pctCorr <- c(pctCorr, cur * 100)
    }
    uq_mat <- cbind(uq_mat, pheno$STATUS, pctCorr)
    spos <- gregexpr("\\/", predFiles)
    fNames <- lapply(seq_len(length(spos)), function(x) {
        n <- length(spos[[x]])
        y <- substr(predFiles[x], spos[[x]][n - 1] + 1, spos[[x]][n] - 
            1)
        y
    })
    fNames <- unlist(fNames)
    output_mat <- uq_mat
    rownames(output_mat) <- pheno$ID
    colnames(output_mat) <- c(fNames, "STATUS", "pctCorrect")
    if (plotAccuracy) {
        p <- ggplot(output_mat, aes(x = output_mat$pctCorrect)) + 
            geom_dotplot()
        msg <- sprintf("Patient-level classification accuracy (N=%i)", 
            length(predFiles))
        p <- p + ggtitle(msg)
        p <- p + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 13))
        print(p)
        return(list(predictions = output_mat, plot = p))
    }
    else return(list(predictions = output_mat))
}
getPatientRankings <-
function (pFile, pheno_DF, predClass, plotIt = FALSE, verbose = FALSE) 
{
    dat <- read.table(pFile, sep = "\t", header = TRUE, as.is = TRUE, 
        skip = 1)
    pheno_DF$ID <- as.character(pheno_DF$ID)
    pheno_DF$STATUS <- as.integer(pheno_DF$STATUS == predClass)
    if (verbose) 
        message(sprintf("%i total ; ", nrow(dat)))
    dat <- dat[which(!is.na(dat[, 2])), ]
    if (verbose) 
        message(sprintf("%i non-query entries in PRANK file\n", 
            nrow(dat)))
    midx <- match(dat[, 1], pheno_DF$ID)
    if (all.equal(pheno_DF$ID[midx], dat[, 1]) != TRUE) {
        stop("\tgetPatientRankings:IDs in GM results don't match pheno\n")
    }
    curlbl <- pheno_DF[midx, , drop = FALSE]
    curlbl <- cbind(curlbl, similarityScore = order(dat[, 2])/nrow(dat), 
        IsPredClass = curlbl$STATUS)
    curlbl <- curlbl[, -which(colnames(curlbl) %in% c("STATUS", 
        "TT_STATUS"))]
    pred <- NA
    perf <- NA
    auc <- NA
    precall <- NA
    f <- NA
    if (nrow(curlbl) >= 2 && length(unique(curlbl$IsPredClass)) == 
        2) {
        pred <- prediction(curlbl$similarityScore, curlbl$IsPredClass)
        perf <- performance(pred, "tpr", "fpr")
        auc <- performance(pred, "auc")@y.values[[1]]
        precall <- performance(pred, "prec", "rec")
        f <- performance(pred, "f")
        if (plotIt) {
            plot(perf, main = sprintf("%i predictions; AUC= %1.2f", 
                nrow(curlbl), auc), bty = "n", las = 1, cex.axis = 1.3)
        }
    }
    out <- merge(x = curlbl, y = pheno_DF, by = "ID", all.y = TRUE)
    return(list(predLbl = curlbl$similarityScore, realLbl = curlbl$IsPredClass, 
        fullmat = out, pred = pred, roc = perf, auc = auc, precall = precall, 
        f = f))
}
getPerformance <-
function (res, predClasses) 
{
    prauc <- function(res) {
        x <- res@x.values[[1]]
        y <- res@y.values[[1]]
        idx <- which(is.nan(y))
        if (any(idx)) {
            x <- x[-idx]
            y <- y[-idx]
        }
        pracma::trapz(x, y)
    }
    pred_col1 <- sprintf("%s_SCORE", predClasses[1])
    pred_col2 <- sprintf("%s_SCORE", predClasses[2])
    idx1 <- which(colnames(res) == pred_col1)
    idx2 <- which(colnames(res) == pred_col2)
    pred <- ROCR::prediction(res[, idx1] - res[, idx2], res$STATUS == 
        predClasses[1])
    st <- res$STATUS
    c1 <- predClasses[1]
    tp <- sum(res$STATUS == res$PRED_CLASS & res$STATUS == c1)
    tn <- sum(res$STATUS == res$PRED_CLASS & res$STATUS != c1)
    fp <- sum(res$STATUS != res$PRED_CLASS & res$STATUS != c1)
    fn <- sum(res$STATUS != res$PRED_CLASS & res$STATUS == c1)
    curRoc <- ROCR::performance(pred, "tpr", "fpr")
    curPr <- ROCR::performance(pred, "prec", "rec")
    tmp <- data.frame(score = 0, tp = tp, tn = tn, fp = fp, fn = fn)
    auroc <- ROCR::performance(pred, "auc")@y.values[[1]]
    aupr <- prauc(curPr)
    corr <- sum(res$STATUS == res$PRED_CLASS)
    acc <- (corr/nrow(res)) * 100
    return(list(rocCurve = curRoc, prCurve = curPr, auroc = auroc, 
        aupr = aupr, accuracy = acc))
}
getPSN <-
function (dat, groupList, makeNetFunc = NULL, sims = NULL, selectedFeatures, 
    plotCytoscape = FALSE, aggFun = "MEAN", prune_pctX = 0.29999999999999999, 
    prune_useTop = TRUE, numCores = 1L, calcShortestPath = FALSE) 
{
    x <- checkMakeNetFuncSims(makeNetFunc = makeNetFunc, sims = sims, 
        groupList = groupList)
    topPath <- gsub(".profile", "", unique(unlist(selectedFeatures)))
    topPath <- gsub("_cont.txt", "", topPath)
    g2 <- list()
    s2 <- list()
    for (nm in names(groupList)) {
        cur <- groupList[[nm]]
        idx <- which(names(cur) %in% topPath)
        message(sprintf("%s: %i features", nm, length(idx)))
        if (length(idx) > 0) {
            g2[[nm]] <- cur[idx]
            s2[[nm]] <- sims[[nm]]
        }
    }
    message("* Making integrated PSN")
    psn <- plotIntegratedPatientNetwork(dataList = dat, groupList = g2, 
        makeNetFunc = makeNetFunc, sims = s2, aggFun = aggFun, 
        prune_pctX = prune_pctX, prune_useTop = prune_useTop, 
        numCores = numCores, calcShortestPath = calcShortestPath, 
        showStats = FALSE, verbose = TRUE, plotCytoscape = plotCytoscape)
    return(psn)
}
getRegionOL <-
function (gr, rngList) 
{
    rng <- GRanges()
    for (k in seq_len(length(rngList))) {
        cur <- rngList[[k]]
        seqlevels(rng) <- unique(c(seqlevels(rng), seqlevels(cur)))
        rng <- c(rng, cur)
    }
    tmp <- as.character(seqlevels(gr))
    rng <- rng[which(as.character(seqnames(rng)) %in% tmp)]
    seqlevels(rng) <- seqlevels(gr)
    ol <- findOverlaps(gr, rng)
    ol <- cbind(queryHits(ol), subjectHits(ol))
    ol_nm <- rng$name[ol[, 2]]
    LOCUS_NAMES <- rep("", length(gr))
    t0 <- Sys.time()
    for (k in unique(ol[, 1])) {
        idx <- which(ol[, 1] == k)
        LOCUS_NAMES[k] <- paste(unique(ol_nm[idx]), collapse = ",")
    }
    print(Sys.time() - t0)
    gr$LOCUS_NAMES = LOCUS_NAMES
    gr
}
getResults <-
function (res, status, featureSelCutoff = 1L, featureSelPct = 0) 
{
    numSplits <- length(grep("^Split", names(res)))
    st <- status
    message(sprintf("Detected %i splits and %i classes", numSplits, 
        length(st)))
    acc <- c()
    predList <- list()
    featScores <- list()
    for (cur in unique(st)) featScores[[cur]] <- list()
    for (k in 1:numSplits) {
        pred <- res[[sprintf("Split%i", k)]][["predictions"]]
        tmp <- pred[, c("ID", "STATUS", "TT_STATUS", "PRED_CLASS", 
            sprintf("%s_SCORE", st))]
        predList[[k]] <- tmp
        acc <- c(acc, sum(tmp$PRED == tmp$STATUS)/nrow(tmp))
        for (cur in unique(st)) {
            tmp <- res[[sprintf("Split%i", k)]][["featureScores"]][[cur]]
            colnames(tmp) <- c("PATHWAY_NAME", "SCORE")
            featScores[[cur]][[sprintf("Split%i", k)]] <- tmp
        }
    }
    auroc <- NULL
    aupr <- NULL
    #if (length(st) == 2) {
    #    message("* Plotting performance")
    #    predPerf <- plotPerf(predList, predClasses = st)
    #    auroc <- unlist(lapply(predPerf, function(x) x$auroc))
    #    aupr <- unlist(lapply(predPerf, function(x) x$aupr))
    #}
    message("* Compiling feature scores and calling selected features")
    feats <- callOverallSelectedFeatures(featScores, featureSelCutoff = featureSelCutoff, 
        featureSelPct = featureSelPct, cleanNames = TRUE)
    return(list(selectedFeatures = feats$selectedFeatures, featureScores = feats$featScores, 
        performance = list(meanAccuracy = mean(acc), splitAccuracy = acc))) 
    #        splitAUROC = auroc, splitAUPR = aupr)))
}
getSimilarity <-
function (x, type = "pearson", customFunc, ...) 
{
    switch(type, pearson = round(cor(na.omit(x), method = "pearson"), 
        digits = 3), custom = customFunc(x, ...))
}
makeInputForEnrichmentMap <-
function (model, results, pathwayList, EMapMinScore = 0L, EMapMaxScore = 1L, 
    EMapPctPass = 0.5, outDir) 
{
    featScores <- results$featureScores
    message("* Creating input files for EnrichmentMap")
    Emap_res <- getEMapInput_many(featScores, pathwayList, minScore = EMapMinScore, 
        maxScore = EMapMaxScore, pctPass = EMapPctPass, model$inputNets, 
        verbose = FALSE)
    gmtFiles <- list()
    nodeAttrFiles <- list()
    message("* Writing files for network visualization")
    for (g in names(Emap_res)) {
        outFile <- paste(outDir, sprintf("%s_nodeAttrs.txt", 
            g), sep = "/")
        write.table(Emap_res[[g]][["nodeAttrs"]], file = outFile, 
            sep = "\t", col.names = TRUE, row.names = FALSE, 
            quote = FALSE)
        nodeAttrFiles[[g]] <- outFile
        outFile <- paste(outDir, sprintf("%s.gmt", g), sep = "/")
        conn <- suppressWarnings(suppressMessages(base::file(outFile, 
            "w")))
        tmp <- Emap_res[[g]][["featureSets"]]
        gmtFiles[[g]] <- outFile
        for (cur in names(tmp)) {
            curr <- sprintf("%s\t%s\t%s", cur, cur, paste(tmp[[cur]], 
                collapse = "\t"))
            writeLines(curr, con = conn)
        }
        close(conn)
    }
    return(list(GMTfiles = gmtFiles, NodeStyles = nodeAttrFiles))
}
makePSN_NamedMatrix <-
function (xpr, nm, namedSets, outDir = tempdir(), simMetric = "pearson", 
    verbose = TRUE, numCores = 1L, writeProfiles = TRUE, sparsify = FALSE, 
    useSparsify2 = FALSE, cutoff = 0.29999999999999999, sparsify_edgeMax = Inf, 
    sparsify_maxInt = 50, minMembers = 1L, runSerially = FALSE, 
    ...) 
{
    if ((!simMetric %in% c("pearson", "MI")) & writeProfiles == 
        TRUE) {
        print(simMetric)
        stop(paste("writeProfiles must only be TRUE with simMetric", 
            " set to pearson or MI. For all other metrics, ", 
            "set writeProfiles=FALSE", sep = ""))
    }
    cl <- makeCluster(numCores, outfile = paste(outDir, "makePSN_log.txt", 
        sep = "/"))
    if (!runSerially) {
        registerDoParallel(cl)
    }
    else {
        message("running serially")
    }
    if (simMetric == "pearson") {
        message(paste("Pearson similarity chosen - ", "enforcing min. 5 patients per net.", 
            sep = ""))
        minMembers <- 5
    }
    `%myinfix%` <- ifelse(runSerially, `%do%`, `%dopar%`)
    outFiles <- foreach(curSet = names(namedSets)) %myinfix% 
        {
            if (verbose) 
                message(sprintf("%s: ", curSet))
            idx <- which(nm %in% namedSets[[curSet]])
            if (verbose) 
                message(sprintf("%i members", length(idx)))
            oFile <- NULL
            if (length(idx) >= minMembers) {
                if (writeProfiles) {
                  outFile <- paste(outDir, sprintf("%s.profile", 
                    curSet), sep = "/")
                  write.table(t(xpr[idx, , drop = FALSE]), file = outFile, 
                    sep = "\t", dec = ".", col.names = FALSE, 
                    row.names = TRUE, quote = FALSE)
                }
                else {
                  outFile <- paste(outDir, sprintf("%s_cont.txt", 
                    curSet), sep = "/")
                  message(sprintf("computing sim for %s", curSet))
                  sim <- getSimilarity(xpr[idx, , drop = FALSE], 
                    type = simMetric, ...)
                  if (is.null(sim)) {
                    stop(sprintf(paste("makePSN_NamedMatrix:%s: ", 
                      "similarity matrix is empty (NULL).\n", 
                      "Check that there isn't a mistake in the ", 
                      "input data or similarity method of choice.\n", 
                      sep = ""), curSet))
                  }
                  pat_pairs <- sim
                  if (sparsify) {
                    if (useSparsify2) {
                      tryCatch({
                        spmat <- sparsify2(pat_pairs, cutoff = cutoff, 
                          EDGE_MAX = sparsify_edgeMax, outFile = outFile, 
                          maxInt = sparsify_maxInt)
                      }, error = function(ex) {
                        stop("sparsify2 caught error\n")
                      })
                    }
                    else {
                      message("sparsify3")
                      tryCatch({
                        sp_t0 <- Sys.time()
                        spmat <- sparsify3(pat_pairs, cutoff = cutoff, 
                          EDGE_MAX = sparsify_edgeMax, outFile = outFile, 
                          maxInt = sparsify_maxInt, verbose = FALSE)
                        print(Sys.time() - sp_t0)
                      }, error = function(ex) {
                        stop("sparsify3 caught error\n")
                      })
                    }
                  }
                  else {
                    write.table(pat_pairs, file = outFile, sep = "\t", 
                      col.names = FALSE, row.names = FALSE, quote = FALSE)
                    print(basename(outFile))
                    message("done")
                  }
                }
                oFile <- basename(outFile)
            }
            oFile
        }
    stopCluster(cl)
    outFiles
}
makePSN_RangeSets <-
function (gr, rangeSet, netDir = tempdir(), simMetric = "coincide", 
    quorum = 2L, verbose = TRUE, numCores = 1L) 
{
    if (!file.exists(netDir)) 
        dir.create(netDir)
    TEST_MODE <- FALSE
    if (TEST_MODE) 
        verbose <- TRUE
    netCountFile <- paste(netDir, "patient_count.txt", sep = "/")
    incPatientFile <- paste(netDir, "inc_patients.txt", sep = "/")
    uq_loci <- unique(unlist(lapply(rangeSet, function(x) {
        x$name
    })))
    uq_patients <- unique(gr$ID)
    if (!simMetric %in% "coincide") 
        stop("Only value supported for simMetric is 'coincide'")
    if (!"LOCUS_NAMES" %in% names(elementMetadata(gr))) {
        message(paste("\tLOCUS_NAMES column not provided; computing ", 
            "overlap of patients\t\twith regions", sep = ""))
        gr <- getRegionOL(gr, rangeSet)
    }
    message("* Preparing patient-locus matrix\n")
    message(sprintf("\t%i unique patients, %i unique locus symbols\n", 
        length(uq_patients), length(uq_loci)))
    bkFile <- paste(tempdir(), "pgmat.bk", sep = "/")
    descFile <- paste(tempdir(), "pgmat.desc", sep = "/")
    if (file.exists(bkFile)) 
        unlink(bkFile)
    if (file.exists(descFile)) 
        unlink(descFile)
    pgMat <- big.matrix(nrow = length(uq_patients), ncol = length(uq_loci), 
        type = "integer", backingpath = tempdir(), backingfile = "pgmat.bk", 
        descriptorfile = "pgmat.desc")
    pgDesc <- describe(pgMat)
    cl <- makeCluster(numCores, outfile = "")
    registerDoParallel(cl)
    num <- length(uq_patients)
    ckSize <- 50
    t0 <- Sys.time()
    x <- foreach(spos = seq(1, num, ckSize)) %dopar% {
        inner_mat <- bigmemory::attach.big.matrix(pgDesc)
        epos <- spos + (ckSize - 1)
        for (k in spos:epos) {
            idx <- which(gr$ID %in% uq_patients[k])
            myloci <- unlist(strsplit(gr$LOCUS_NAMES[idx], ","))
            myloci <- setdiff(myloci, "")
            if (length(myloci) > 0) {
                inner_mat[k, which(uq_loci %in% myloci)] <- 1L
            }
            if (k%%100 == 0) 
                message(".")
        }
    }
    t1 <- Sys.time()
    print(t1 - t0)
    hit_p <- integer(length(rangeSet))
    inc_patients <- integer(length(uq_patients))
    names(inc_patients) <- uq_patients
    message("* Writing networks\n")
    `%myinfix%` <- ifelse(TEST_MODE, `%do%`, `%dopar%`)
    t0 <- Sys.time()
    outFiles <- foreach(idx = seq_len(length(rangeSet))) %myinfix% 
        {
            curP <- names(rangeSet)[idx]
            if (verbose) 
                message(sprintf("\t%s: ", curP))
            inner_mat <- bigmemory::attach.big.matrix(pgDesc)
            locus_idx <- which(uq_loci %in% rangeSet[[idx]]$name)
            if (length(locus_idx) >= 2) {
                hit_pathway <- rowSums(inner_mat[, locus_idx])
            }
            else {
                hit_pathway <- inner_mat[, locus_idx]
            }
            hit_p[idx] <- sum(hit_pathway > 0)
            if (verbose) 
                message(sprintf("%i patients with interactions", 
                  hit_p[idx]))
            pScore <- 1
            outFile <- ""
            if (hit_p[idx] >= quorum) {
                if (verbose) 
                  message(sprintf("\n\t\tlength=%i; score = %1.2f", 
                    length(rangeSet[[idx]]), pScore))
                x <- hit_pathway
                inc_patients[x > 0] <- inc_patients[x > 0] + 
                  1
                tmp <- uq_patients[hit_pathway > 0]
                pat_pairs <- t(combinat::combn(tmp, 2))
                pat_pairs <- cbind(pat_pairs, pScore)
                outFile <- paste(netDir, sprintf("%s_cont.txt", 
                  curP), sep = "/")
                write.table(pat_pairs, file = outFile, sep = "\t", 
                  col.names = FALSE, row.names = FALSE, quote = FALSE)
                outFile <- basename(outFile)
            }
            if (idx%%100 == 0) 
                message(".")
            if (verbose) 
                message("\n")
            if (!verbose) {
                if (idx%%100 == 0) 
                  message(".")
                if (idx%%1000 == 0) 
                  message("\n")
            }
            outFile
        }
    t1 <- Sys.time()
    print(t1 - t0)
    stopCluster(cl)
    if (file.exists(bkFile)) 
        unlink(bkFile)
    if (file.exists(descFile)) 
        unlink(descFile)
    outFiles <- unlist(outFiles)
    outFiles <- outFiles[which(outFiles != "")]
    outFiles
}
makeQueries <-
function (incPat, featScoreMax = 10L, verbose = TRUE) 
{
    incPat <- sample(incPat, replace = FALSE)
    num2samp <- floor(((featScoreMax - 1)/featScoreMax) * length(incPat))
    csize <- round((1/featScoreMax) * length(incPat))
    if (verbose) {
        message(sprintf("\t\t%i IDs; %i queries (%i sampled, %i test)", 
            length(incPat), featScoreMax, num2samp, csize))
    }
    out <- list()
    for (k in seq_len(featScoreMax)) {
        sidx <- ((k - 1) * csize) + 1
        eidx <- k * csize
        if (k == featScoreMax) 
            eidx <- length(incPat)
        p1 <- sprintf("\t\tQ%i: %i test; ", k, eidx - sidx + 
            1)
        out[[k]] <- setdiff(incPat, incPat[sidx:eidx])
        if (verbose) 
            message(sprintf("%s %i query", p1, length(out[[k]])))
    }
    out
}
makeSymmetric <-
function (x, verbose = FALSE) 
{
    samps <- unique(c(x[, 1], x[, 2]))
    newmat <- matrix(NA, nrow = length(samps), ncol = length(samps))
    rownames(newmat) <- samps
    colnames(newmat) <- samps
    i <- 1
    for (k in samps) {
        idx <- which(x[, 1] == k)
        if (verbose) 
            message(k)
        for (curr in idx) {
            j <- which(colnames(newmat) == x[curr, 2])
            newmat[i, j] <- x[curr, 3]
            newmat[j, i] <- x[curr, 3]
        }
        i <- i + 1
    }
    diag(newmat) <- 1
    return(newmat)
}
mapNamedRangesToSets <-
function (gr, rangeList, verbose = FALSE) 
{
    out <- list()
    for (nm in names(rangeList)) {
        my_gr <- gr[which(gr$name %in% rangeList[[nm]])]
        if (verbose) 
            message(sprintf("%s: %i ranges\n", nm, length(my_gr)))
        out[[nm]] <- my_gr
    }
    out
}
normDiff <-
function (x) 
{
    nm <- colnames(x)
    x <- as.numeric(x)
    n <- length(x)
    rngX <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    out <- matrix(NA, nrow = n, ncol = n)
    for (j in seq_len(n)) out[, j] <- 1 - (abs((x - x[j])/rngX))
    rownames(out) <- nm
    colnames(out) <- nm
    out
}
perfCalc <-
function (dat) 
{
    dat <- na.omit(dat)
    tp2 <- 2 * dat$tp
    f1 <- tp2/(tp2 + dat$fp + dat$fn)
    ppv <- dat$tp/(dat$tp + dat$fp)
    rec <- dat$tp/(dat$tp + dat$fn)
    prauc <- pracma::trapz(rev(rec), rev(ppv))
    x <- dat$fp/(dat$fp + dat$tn)
    y <- dat$tp/(dat$tp + dat$fn)
    x <- c(0, rev(x), 1)
    y <- c(0, rev(y), 1)
    auc <- pracma::trapz(x, y)
    out <- data.frame(score = dat$score, ppv = ppv, f1 = f1, 
        rec = rec)
    return(list(stats = out, auc = auc, prauc = prauc))
}
plotEmap <-
function (gmtFile, nodeAttrFile, netName = "generic", scoreCol = "maxScore", 
    minScore = 1, maxScore = 10, nodeFillStops = c(7, 9), colorScheme = "cont_heatmap", 
    imageFormat = "png", verbose = FALSE, createStyle = TRUE, 
    groupClusters = FALSE, hideNodeLabels = FALSE) 
{
    if (!requireNamespace("RCy3", quietly = TRUE)) {
        stop("Package \"RCy3\" needed for plotEmap() to work. Please install it and then make your call.", 
            call. = FALSE)
    }
    validColSchemes <- c("cont_heatmap", "netDx_ms")
    if (!colorScheme %in% validColSchemes) {
        stop(sprintf("colorScheme should be one of { %s }\n", 
            paste(validColSchemes, collapse = ",")))
    }
    if (netName %in% RCy3::getNetworkList()) {
        RCy3::deleteNetwork(netName)
    }
    em_command <- paste("enrichmentmap build analysisType=\"generic\"", 
        "gmtFile=", gmtFile, "pvalue=", 1, "qvalue=", 1, "similaritycutoff=", 
        0.050000000000000003, "coefficients=", "JACCARD")
    response <- RCy3::commandsGET(em_command)
    RCy3::renameNetwork(netName, RCy3::getNetworkSuid())
    aa_command <- paste("autoannotate annotate-clusterBoosted", 
        "clusterAlgorithm=MCL", "labelColumn=name", "maxWords=3", 
        "network=", netName)
    print(aa_command)
    response <- RCy3::commandsGET(aa_command)
    message("* Importing node attributes\n")
    table_command <- sprintf(paste("table import file file=%s ", 
        "keyColumnIndex=1 ", "firstRowAsColumnNames=true startLoadRow=1 TargetNetworkList=%s ", 
        "WhereImportTable=To%%20selected%%20networks%%20only", 
        sep = " "), nodeAttrFile, netName)
    response <- RCy3::commandsGET(table_command)
    message("* Creating or applying style\n")
    all_unique_scores_int <- sort(unique(read.delim(nodeAttrFile)[, 
        2]))
    all_unique_scores <- unlist(lapply(all_unique_scores_int, 
        toString))
    styleName <- "EMapStyle"
    scoreVals <- minScore:maxScore
    style_cols <- ""
    if (colorScheme == "cont_heatmap") {
        colfunc <- colorRampPalette(c("yellow", "red"))
        gradient_cols <- colfunc(length(scoreVals))
        style_cols <- colfunc(length(scoreVals))
    }
    else if (colorScheme == "netDx_ms") {
        style_cols <- rep("white", length(scoreVals))
        style_cols[which(scoreVals >= nodeFillStops[1])] <- "orange"
        style_cols[which(scoreVals >= nodeFillStops[2])] <- "red"
    }
    nodeLabels <- RCy3::mapVisualProperty("node label", "name", 
        "p")
    nodeFills <- RCy3::mapVisualProperty("node fill color", scoreCol, 
        "d", scoreVals, style_cols)
    defaults <- list(NODE_SHAPE = "ellipse", NODE_SIZE = 30, 
        EDGE_TRANSPARENCY = 200, NODE_TRANSPARENCY = 255, EDGE_STROKE_UNSELECTED_PAINT = "#999999")
    if (createStyle) {
        message("Making style\n")
        RCy3::createVisualStyle(styleName, defaults, list(nodeLabels, 
            nodeFills))
    }
    RCy3::setVisualStyle(styleName)
    if (groupClusters) {
        RCy3::layoutNetwork("attributes-layout NodeAttribute=__mclCLuster")
        redraw_command <- sprintf("autoannotate redraw network=%s", 
            RCy3::getNetworkSuid())
        response <- RCy3::commandsGET(redraw_command)
        RCy3::fitContent()
        redraw_command <- sprintf("autoannotate redraw network=%s", 
            RCy3::getNetworkSuid())
        response <- RCy3::commandsGET(redraw_command)
        RCy3::fitContent()
    }
    if (hideNodeLabels) {
        RCy3::setNodeFontSizeDefault(0, styleName)
    }
}
plotIntegratedPatientNetwork <-
function (dataList, groupList, makeNetFunc = NULL, sims = NULL, 
    setName = "predictor", prune_pctX = 0.050000000000000003, 
    prune_useTop = TRUE, aggFun = "MAX", calcShortestPath = FALSE, 
    showStats = FALSE, outDir = tempdir(), numCores = 1L, nodeSize = 50L, 
    edgeTransparency = 40L, nodeTransparency = 155L, plotCytoscape = FALSE, 
    verbose = FALSE) 
{
    checkMakeNetFuncSims(makeNetFunc = makeNetFunc, sims = sims, 
        groupList = groupList)
    if (missing(dataList)) 
        stop("dataList is missing.")
    dat <- dataList2List(dataList, groupList)
    pheno <- dat$pheno[, c("ID", "STATUS")]
    if (!file.exists(outDir)) 
        dir.create(outDir)
    profDir <- paste(outDir, "profiles", sep = "/")
    if (!file.exists(profDir)) 
        dir.create(profDir)
    pheno_id <- setupFeatureDB(pheno, outDir)
    createPSN_MultiData(dataList = dat$assays, groupList = groupList, 
        pheno = pheno_id, netDir = outDir, makeNetFunc = makeNetFunc, 
        sims = sims, numCores = numCores, verbose = FALSE)
    convertProfileToNetworks(netDir = profDir, outDir = paste(outDir, 
        "INTERACTIONS", sep = "/"))
    predClasses <- unique(pheno$STATUS)
    colnames(pheno)[which(colnames(pheno) == "STATUS")] <- "GROUP"
    pid <- read.delim(paste(outDir, "GENES.txt", sep = "/"), 
        sep = "\t", header = FALSE, as.is = TRUE)[, 1:2]
    colnames(pid)[1:2] <- c("GM_ID", "ID")
    netid <- read.delim(paste(outDir, "NETWORKS.txt", sep = "/"), 
        sep = "\t", header = FALSE, as.is = TRUE)
    colnames(netid)[1:2] <- c("NET_ID", "NETWORK")
    message("* Computing aggregate net")
    out <- writeWeightedNets(pid, netIDs = netid, netDir = paste(outDir, 
        "INTERACTIONS", sep = "/"), filterEdgeWt = 0, 
        limitToTop = Inf, plotEdgeDensity = FALSE, aggNetFunc = aggFun, 
        verbose = FALSE)
    aggNet <- out$aggNet
    aggNet <- aggNet[, 1:3]
    distNet <- aggNet
    distNet[, 3] <- 1 - distNet[, 3]
    if (calcShortestPath) {
        x <- compareShortestPath(distNet, pheno, verbose = showStats, 
            plotDist = TRUE)
        gp <- unique(pheno$GROUP)
        oppName <- paste(gp[1], gp[2], sep = "-")
        curDijk <- matrix(NA, nrow = 1, ncol = 6)
        colnames(curDijk) <- c(gp[1], gp[2], oppName, "overall", 
            sprintf("p%s-Opp", gp[1]), sprintf("p%s-Opp", gp[2]))
        if (showStats) {
            message("Shortest path averages &")
            message("p-values (one-sided WMW)")
            message("------------------------------------")
        }
        for (k in 1:length(x$all)) {
            cur <- names(x$all)[k]
            idx <- which(colnames(curDijk) %in% cur)
            curDijk[1, idx] <- median(x$all[[k]])
            if (showStats) 
                message(sprintf("\t%s: Median = %1.2f ", cur, 
                  curDijk[1, idx]))
            if (cur %in% gp) {
                tmp <- wilcox.test(x$all[[cur]], x$all[[oppName]], 
                  alternative = "less")$p.value
                curDijk[4 + k] <- tmp
                if (showStats) 
                  message(sprintf(" p(<Opp) = %1.2e\n", curDijk[4 + 
                    k]))
            }
        }
    }
    message("")
    message("* Prune network")
    colnames(aggNet) <- c("source", "target", "weight")
    aggNet_pruned <- pruneNet_pctX(aggNet, pheno$ID, pctX = prune_pctX, 
        useTop = prune_useTop)
    if (plotCytoscape) {
        if (!requireNamespace("RCy3", quietly = TRUE)) {
            stop("Package \"RCy3\" needed for this function to work. Please install it.", 
                call. = FALSE)
        }
        message("* Creating network in Cytoscape")
        colnames(pheno)[which(colnames(pheno) == "ID")] <- "id"
        colnames(aggNet_pruned) <- c("source", "target", "weight")
        if (prune_useTop) {
            topTxt <- "top"
        }
        else {
            topTxt <- "bottom"
        }
        RCy3::createNetworkFromDataFrames(nodes = pheno, edges = aggNet_pruned, 
            title = sprintf("%s_%s_%s%1.2f", setName, aggFun, 
                topTxt, prune_pctX), collName = setName)
        styleName <- "PSNstyle"
        pal <- suppressWarnings(brewer.pal(name = "Dark2", n = length(predClasses)))
        message("* Creating style")
        colLegend <- data.frame(STATUS = predClasses, colour = pal[1:length(predClasses)], 
            stringsAsFactors = FALSE)
        nodeFills <- RCy3::mapVisualProperty("node fill color", 
            "GROUP", mapping.type = "d", table.column.values = predClasses, 
            visual.prop.values = pal)
        defaults <- list(NODE_SHAPE = "ellipse", NODE_SIZE = nodeSize, 
            EDGE_TRANSPARENCY = edgeTransparency, EDGE_STROKE_UNSELECTED_PAINT = "#999999", 
            NODE_TRANSPARENCY = nodeTransparency)
        sty <- RCy3::createVisualStyle(styleName, defaults, list(nodeFills))
        RCy3::setVisualStyle(styleName)
        message("* Applying layout\n")
        RCy3::layoutNetwork("kamada-kawai column=weight")
        RCy3::setVisualStyle(styleName)
        out <- list(patientSimNetwork_unpruned = aggNet, patientSimNetwork_pruned = aggNet_pruned, 
            colLegend = colLegend, outDir = outDir)
    }
    else {
        message(paste("* plotCytoscape is set to FALSE.", "Set to TRUE to visualize patient network in Cytoscape", 
            sep = ""))
        out <- list(patientSimNetwork_unpruned = aggNet, patientDistNetwork_pruned = aggNet_pruned, 
            aggFun = aggFun, outDir = outDir)
    }
    return(out)
}
plotPerf <-
function (resList = NULL, inFiles, predClasses, plotSEM = FALSE) 
{
    if (is.null(resList)) {
        if (missing(inFiles)) 
            stop("inDir not provided")
    }
    if (missing(predClasses)) 
        stop("predClasses missing; please specify classes")
    prauc <- function(dat) {
        x <- dat@x.values[[1]]
        y <- dat@y.values[[1]]
        idx <- which(is.nan(y))
        if (any(idx)) {
            x <- x[-idx]
            y <- y[-idx]
        }
        pracma::trapz(x, y)
    }
    if (is.null(resList)) {
        resList <- list()
        ctr <- 1
        for (fName in inFiles) {
            resList[[ctr]] <- read.delim(fName, sep = "\t", header = TRUE, 
                as.is = TRUE)
            ctr <- ctr + 1
        }
    }
    mega <- list()
    for (ctr in seq_len(length(resList))) {
        dat <- resList[[ctr]]
        out <- list()
        overall_acc <- numeric()
        curRoc <- list()
        curPr <- list()
        pred_col1 <- sprintf("%s_SCORE", predClasses[1])
        pred_col2 <- sprintf("%s_SCORE", predClasses[2])
        idx1 <- which(colnames(dat) == pred_col1)
        idx2 <- which(colnames(dat) == pred_col2)
        pred <- ROCR::prediction(dat[, idx1] - dat[, idx2], dat$STATUS == 
            predClasses[1])
        c1 <- predClasses[1]
        tp <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS == 
            c1)
        tn <- sum(dat$STATUS == dat$PRED_CLASS & dat$STATUS != 
            c1)
        fp <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS != 
            c1)
        fn <- sum(dat$STATUS != dat$PRED_CLASS & dat$STATUS == 
            c1)
        curRoc <- ROCR::performance(pred, "tpr", "fpr")
        curPr <- ROCR::performance(pred, "prec", "rec")
        tmp <- data.frame(score = 0, tp = tp, tn = tn, fp = fp, 
            fn = fn)
        out <- perfCalc(tmp)
        auroc <- performance(pred, "auc")@y.values[[1]]
        aupr <- prauc(curPr)
        corr <- sum(dat$STATUS == dat$PRED_CLASS)
        overall_acc <- c(overall_acc, corr/nrow(dat) * 100)
        mega[[ctr]] <- list(stats = out$stats, roc_curve = curRoc, 
            pr_curve = curPr, auroc = auroc, aupr = aupr, accuracy = overall_acc)
    }
    .plotAvg <- function(res, name, plotSEM) {
        mu <- mean(res, na.rm = TRUE)
        if (plotSEM) {
            err <- sd(res, na.rm = TRUE)/sqrt(length(res))
            errnm <- "SEM"
        }
        else {
            err <- sd(res, na.rm = TRUE)
            errnm <- "SD"
        }
        plot(1, mu, type = "n", bty = "n", ylab = sprintf("%s (mean+/-%s)", 
            name, errnm), xaxt = "n", ylim = c(0.40000000000000002, 
            1), las = 1, xlim = c(0.80000000000000004, 1.2), 
            cex.axis = 1.3999999999999999, xlab = "")
        abline(h = c(0.69999999999999996, 0.80000000000000004), 
            col = "cadetblue3", lty = 3, lwd = 3)
        points(1, mu, type = "p", cex = 1.3999999999999999, pch = 16)
        segments(x0 = 1, y0 = mu - err, y1 = mu + err, lwd = 3)
        segments(x0 = 1 - 0.01, x1 = 1 + 0.01, y0 = mu - err, 
            y1 = mu - err)
        segments(x0 = 1 - 0.01, x1 = 1 + 0.01, y0 = mu + err, 
            y1 = mu + err)
        abline(h = 0.5, col = "red", lty = 1, lwd = 2)
        title(sprintf("%s: N=%i runs", name, length(res)))
    }
    par(mfrow = c(2, 2))
    x <- unlist(lapply(mega, function(x) x$auroc))
    .plotAvg(x, "AUROC", plotSEM)
    x <- unlist(lapply(mega, function(x) x$aupr))
    .plotAvg(x, "AUPR", plotSEM)
    rocCurves <- lapply(mega, function(x) x$roc_curve)
    plotPerf_multi(rocCurves, "ROC")
    prCurves <- lapply(mega, function(x) x$pr_curve)
    plotPerf_multi(prCurves, "PR", plotType = "PR")
    return(mega)
}
plotPerf_multi <-
function (inList, plotTitle = "performance", plotType = "ROC", 
    xlab = "TPR", ylab = "FPR", meanCol = "darkblue", xlim = c(0, 
        1), ylim = c(0, 1)) 
{
    if (class(inList) == "performance") {
        inList <- list(inList)
    }
    if (plotType == "ROC") {
        xlab <- "TPR"
        ylab <- "FPR"
    }
    else if (plotType == "PR") {
        xlab <- "Precision"
        ylab <- "Recall"
    }
    else {
        message("custom type plot\n")
    }
    plot(0, 0, type = "n", bty = "n", las = 1, xlim = xlim, ylim = ylim, 
        xlab = xlab, ylab = ylab, main = plotTitle, cex.axis = 1.8, 
        cex.lab = 1.6000000000000001)
    out <- list()
    is_empty <- 0
    for (k in seq_len(length(inList))) {
        if (length(slotNames(inList[[k]])) == 0) {
            is_empty <- is_empty + 1
            next
        }
        x <- inList[[k]]@x.values[[1]]
        y <- inList[[k]]@y.values[[1]]
        cur <- aggregate(y, by = list(xvals = x), FUN = mean, 
            na.rm = TRUE)
        colnames(cur) <- c("x", "y")
        out[[k]] <- cur
        points(x, y, type = "l", col = "gray90", lwd = 3)
    }
    x <- inList[[k]]@x.values[[1]]
    y <- inList[[k]]@y.values[[1]]
    cur_y <- aggregate(y, by = list(xvals = x), FUN = mean, na.rm = TRUE)
    cur <- cbind(cur_y, k)
    colnames(cur) <- c("x", "y", "k")
    out[[k]] <- cur
    points(x, y, type = "l", col = meanCol, lwd = 4)
    if (length(inList) > 1) {
        text(0.80000000000000004 * xlim[2], 0.10000000000000001 * 
            ylim[2], sprintf("%i splits", length(inList) - is_empty), 
            cex = 1.3)
    }
    if (plotType == "ROC") 
        abline(0, 1, col = "red", lwd = 3)
    else if (plotType == "PR") 
        abline(h = 0.5, col = "red", lwd = 3)
}
predict <-
function (trainMAE, testMAE, groupList, selectedFeatures, makeNetFunc = NULL, 
    sims = NULL, outDir, verbose = FALSE, numCores = 1L, JavaMemory = 4L, 
    debugMode = FALSE) 
{
    if (missing(trainMAE)) 
        stop("trainMAE must be supplied.\n")
    if (missing(testMAE)) 
        stop("testMAE must be supplied.\n")
    if (missing(groupList)) 
        stop("groupList must be supplied.\n")
    if (length(groupList) < 1) 
        stop("groupList must be of length 1+\n")
    if (class(selectedFeatures) != "list") 
        stop("selectedFeatures must be a list with patient labels as keys, and selected features as values")
    if (missing(outDir)) 
        stop("outDir must be supplied.\n")
    if (!is(trainMAE, "MultiAssayExperiment")) 
        stop("trainMAE must be a MultiAssayExperiment")
    if (!is(testMAE, "MultiAssayExperiment")) 
        stop("testMAE must be a MultiAssayExperiment")
    tmp <- unlist(lapply(groupList, class))
    not_list <- sum(tmp == "list") < length(tmp)
    nm1 <- setdiff(names(groupList), "clinical")
    names_nomatch <- any(!nm1 %in% names(trainMAE))
    if (!is(groupList, "list") || not_list || names_nomatch) {
        msg <- c("groupList must be a list of lists.", " Names must match those in trainMAE, and each entry should be a list", 
            " of networks for this group.")
        stop(paste(msg, sep = ""))
    }
    for (nm in names(selectedFeatures)) {
        selectedFeatures[[nm]] <- sub("_cont.txt", "", sub(".profile", 
            "", selectedFeatures[[nm]]))
    }
    fs <- unlist(selectedFeatures)
    names(fs) <- NULL
    gl <- c()
    for (k in names(groupList)) {
        m <- groupList[[k]]
        gl <- c(gl, names(m))
    }
    if (sum(!fs %in% gl) > 0) {
        stop("One or more entry in selectedFeatures was not found in groupList.")
    }
    trainList <- dataList2List(trainMAE, groupList)
    testList <- dataList2List(testMAE, groupList)
    ph <- trainList$pheno[, c("ID", "STATUS")]
    ph2 <- testList$pheno[, c("ID", "STATUS")]
    ph$TT_STATUS <- "TRAIN"
    ph2$TT_STATUS <- "TEST"
    message("* Merging metadata tables...")
    tryCatch({
        pheno <- rbind(ph, ph2)
    }, error = function(ex) {
        stop(paste("couldn't combine train and test pheno.", 
            "check that they have identical columns in same order", 
            sep = ""))
    })
    print(table(pheno[, c("STATUS", "TT_STATUS")]))
    message("* Merging assays ...")
    assays <- list()
    for (nm in names(trainList$assays)) {
        message(sprintf("\t%s", nm))
        tryCatch({
            assays[[nm]] <- cbind(trainList$assays[[nm]], testList$assays[[nm]])
        }, error = function(ex) {
            stop(sprintf(paste("Error while combining data type %s for train and test ", 
                "samples. Have you checked that measures are identical for both?", 
                sep = ""), nm))
        })
    }
    message("* Measuring similarity to each known class")
    subtypes <- unique(ph$STATUS)
    predRes <- list()
    for (g in subtypes) {
        if (verbose) 
            message(sprintf("\t%s", g))
        pDir <- paste(outDir, g, sep = "/")
        netDir <- paste(pDir, "networks", sep = "/")
        dir.create(pDir)
        dir.create(netDir)
        pheno_id <- setupFeatureDB(pheno, netDir)
        x <- checkMakeNetFuncSims(makeNetFunc = makeNetFunc, 
            sims = sims, groupList = groupList)
        if (verbose) 
            message("Creating PSN")
        createPSN_MultiData(dataList = assays, groupList = groupList, 
            pheno = pheno_id, netDir = netDir, makeNetFunc = makeNetFunc, 
            sims = sims, numCores = 1L, filterSet = selectedFeatures[[g]], 
            verbose = verbose)
        dbDir <- compileFeatures(netDir, outDir = pDir, numCores = numCores, 
            verbose = verbose, debugMode = debugMode)
        qSamps <- pheno$ID[which(pheno$STATUS %in% g & pheno$TT_STATUS %in% 
            "TRAIN")]
        qFile <- paste(pDir, sprintf("%s_query", g), sep = "/")
        message(sprintf("\t%s : %s training samples", g, prettyNum(length(qSamps), 
            big.mark = ",")))
        writeQueryFile(qSamps, "all", nrow(pheno), qFile)
        if (verbose) 
            message(sprintf("\t** %s: Compute similarity", g))
        resFile <- runQuery(dbDir$dbDir, qFile, resDir = pDir, 
            JavaMemory = JavaMemory, numCores = numCores, verbose = verbose, 
            debugMode = debugMode)
        predRes[[g]] <- getPatientRankings(sprintf("%s.PRANK", 
            resFile), pheno, g)
    }
    predClass <- predictPatientLabels(predRes, verbose = verbose)
    out <- merge(x = pheno, y = predClass, by = "ID")
    if (nrow(out) != nrow(colData(testMAE))) {
        warning(paste(rep("*", 25), "Not all patients provided in the test sample were classified.", 
            rep("*", 25), sep = "\n"))
    }
    acc <- sum(out$STATUS == out$PRED_CLASS)/nrow(out)
    message(sprintf("%s test patients", prettyNum(nrow(out), 
        big.mark = ",")))
    message(sprintf("ACCURACY (N=%i test) = %2.1f%%", nrow(out), 
        acc * 100))
    message("Confusion matrix")
    print(table(out[, c("STATUS", "PRED_CLASS")]))
    out <- out[, -which(colnames(out) == "TT_STATUS")]
    return(out)
}
predictPatientLabels <-
function (resSet, verbose = TRUE) 
{
    type_rank <- NULL
    for (k in seq_len(length(resSet))) {
        x <- resSet[[k]]$fullmat
        idx <- which(colnames(x) == "GM_score")
        if (any(idx)) 
            colnames(x)[idx] <- "similarityScore"
        if (is.null(type_rank)) 
            type_rank <- x[, c("ID", "similarityScore")]
        else {
            if (all.equal(x$ID, type_rank$ID) != TRUE) {
                stop("predictPatientLabels: ids don't match")
            }
            type_rank <- cbind(type_rank, x[, "similarityScore"])
        }
        rnkCol <- paste(names(resSet)[k], "SCORE", sep = "_")
        colnames(type_rank)[ncol(type_rank)] <- rnkCol
    }
    na_sum <- rowSums(is.na(type_rank[, -1]))
    if (verbose) {
        if (any(na_sum > 0)) 
            message(sprintf(paste("*** %i rows have an NA prediction ", 
                "(probably query samples that were not not ranked\n", 
                sep = ""), sum(na_sum > 0)))
    }
    type_rank <- na.omit(type_rank)
    maxScore <- rep(NA, nrow(type_rank))
    for (k in seq_len(nrow(type_rank))) {
        maxScore[k] <- colnames(type_rank)[which.max(type_rank[k, 
            -1]) + 1]
    }
    patClass <- sub("_SCORE", "", maxScore)
    type_rank <- cbind(type_rank, PRED_CLASS = patClass)
    type_rank$PRED_CLASS <- as.character(type_rank$PRED_CLASS)
    type_rank
}
pruneNet <-
function (net, vertices, pctX = 0.10000000000000001, useTop = TRUE) 
{
    g <- igraph::graph_from_data_frame(net, vertices = vertices)
    wt <- sort(E(g)$weight, decreasing = TRUE)
    if (useTop) {
        thresh <- wt[length(wt) * pctX]
        g2 <- delete.edges(g, which(E(g)$weight < thresh))
    }
    else {
        thresh <- wt[length(wt) * (1 - pctX)]
        g2 <- delete.edges(g, which(E(g)$weight > thresh))
    }
    df <- as.data.frame(get.edgelist(g2))
    df[, 1] <- as.character(df[, 1])
    df[, 2] <- as.character(df[, 2])
    df$weight <- edge_attr(g2, name = "weight")
    colnames(df) <- c("AliasA", "AliasB", "weight")
    return(df)
}
pruneNet_pctX <-
function (net, vertices, pctX = 0.10000000000000001, useTop = TRUE) 
{
    g <- igraph::graph_from_data_frame(net, vertices = vertices)
    wt <- sort(E(g)$weight, decreasing = TRUE)
    if (useTop) {
        thresh <- wt[length(wt) * pctX]
        g2 <- delete.edges(g, which(E(g)$weight < thresh))
    }
    else {
        thresh <- wt[length(wt) * (1 - pctX)]
        g2 <- delete.edges(g, which(E(g)$weight > thresh))
    }
    df <- as.data.frame(get.edgelist(g2))
    df[, 1] <- as.character(df[, 1])
    df[, 2] <- as.character(df[, 2])
    df$weight <- edge_attr(g2, name = "weight")
    colnames(df) <- c("AliasA", "AliasB", "weight")
    return(df)
}
pruneNets <-
function (oldDir, newDir = tempdir(), filterNets = "*", filterIDs = "*", 
    netSfx = "_cont.txt$", verbose = TRUE) 
{
    if (length(filterNets) == 1) {
        if (filterNets == "*") {
            if (verbose) 
                message("* Including all networks\n")
            fList <- dir(path = oldDir, pattern = netSfx)
            filterNets <- fList
        }
    }
    if (verbose) 
        message(sprintf("Limiting to %i networks\n", length(filterNets)))
    if (!file.exists(newDir)) 
        dir.create(newDir)
    if (length(filterIDs) == 1) {
        if (filterIDs == "*") {
            message("* Including all patients\n")
            for (f in filterNets) {
                oldf <- paste(oldDir, f, sep = "/")
                newf <- paste(newDir, f, sep = "/")
                file.copy(oldf, newf)
            }
        }
    }
    else {
        if (verbose) 
            message(sprintf("Limiting to %i patients\n", length(filterIDs)))
        for (f in filterNets) {
            dat <- read.delim(paste(oldDir, f, sep = "/"), 
                sep = "\t", header = FALSE, as.is = TRUE)
            idx <- intersect(which(dat[, 1] %in% filterIDs), 
                which(dat[, 2] %in% filterIDs))
            write.table(dat[idx, ], file = paste(newDir, f, sep = "/"), 
                sep = "\t", col.names = FALSE, row.names = FALSE, 
                quote = FALSE)
        }
    }
}
randAlphanumString <-
function (numStrings = 1L) 
{
    a <- do.call(paste0, replicate(5, sample(LETTERS, numStrings, 
        TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, numStrings, TRUE)), 
        sample(LETTERS, numStrings, TRUE))
}
readPathways <-
function (fname, MIN_SIZE = 10L, MAX_SIZE = 200L, EXCLUDE_KEGG = TRUE, 
    IDasName = FALSE, verbose = TRUE, getOrigNames = FALSE) 
{
    oldLocale <- Sys.getlocale("LC_ALL")
    Sys.setlocale("LC_ALL", "C")
    out <- list()
    srcList <- list()
    if (verbose) 
        message("---------------------------------------\n")
    if (verbose) 
        message(sprintf("File: %s\n\n", basename(fname)))
    f <- file(fname, "r")
    ctr <- 0
    options(warn = 1)
    repeat {
        s <- scan(f, what = "character", nlines = 1, quiet = TRUE, 
            sep = "\t")
        if (length(s) == 0) 
            break
        pPos <- gregexpr("%", s[1])[[1]]
        src <- ""
        src_id <- ""
        if (pPos[1] == -1) {
            s[1] <- s[1]
        }
        else {
            src <- substr(s[1], pPos[1] + 1, pPos[2] - 1)
            src_id <- substr(s[1], pPos[2] + 1, nchar(s[1]))
            if (IDasName) 
                s[1] <- paste(src, src_id, sep = ":")
            else s[1] <- substr(s[1], 1, pPos[1] - 1)
        }
        if (!EXCLUDE_KEGG || (src != "KEGG")) {
            idx <- which(s == "")
            if (any(idx)) 
                s <- s[-idx]
            out[[s[1]]] <- s[3:length(s)]
            srcList[[s[1]]] <- src
        }
        ctr <- ctr + 1
    }
    close(f)
    if (verbose) {
        message(sprintf(paste("Read %i pathways in total, ", 
            "internal list has %i entries", sep = ""), ctr, length(out)))
        message(sprintf("\tFILTER: sets with num genes in [%i, %i]", 
            MIN_SIZE, MAX_SIZE))
    }
    ln <- unlist(lapply(out, length))
    idx <- which(ln < MIN_SIZE | ln >= MAX_SIZE)
    if (any(idx)) {
        out[idx] <- NULL
        srcList[idx] <- NULL
    }
    if (verbose) 
        message(sprintf("\t  => %i pathways excluded\n\t  => %i left", 
            length(idx), length(out)))
    nm <- suppressMessages(suppressWarnings(cleanPathwayName(names(out))))
    idx <- which(duplicated(nm))
    if (any(idx)) {
        message("Resolving duplicate pathway names by appending source...")
        nm[idx] <- sprintf("%s_%s", nm[idx], srcList[[idx]])
    }
    if (getOrigNames) {
        pnames <- cbind(names(out), nm)
        names(out) <- nm
        out <- list(geneSets = out, pNames = pnames)
    }
    else {
        names(out) <- nm
    }
    return(out)
}
RR_featureTally <-
function (netmat, phenoDF, TT_STATUS, predClass, pScore, outDir = tempdir(), 
    enrichLabels = TRUE, enrichedNets, maxScore = 30L, verbose = FALSE) 
{
    pTally <- list()
    for (k in seq_len(length(TT_STATUS))) {
        dat <- pScore[[k]]
        dat[, 1] <- as.character(dat[, 1])
        for (m in seq_len(nrow(dat))) {
            curp <- dat[m, 1]
            if (!curp %in% names(pTally)) 
                pTally[[curp]] <- 0
            pTally[[curp]] <- pTally[[curp]] + dat[m, 2]
        }
    }
    pathDF <- data.frame(PATHWAY_NAME = names(pTally), SCORE = unlist(pTally))
    pathDF[, 2] <- as.integer(as.character(pathDF[, 2]))
    pathDF <- pathDF[order(pathDF[, 2], decreasing = TRUE), ]
    tmpOut <- paste(outDir, "pathway_cumTally.txt", sep = "/")
    write.table(pathDF, file = tmpOut, sep = "\t", col.names = TRUE, 
        row.names = FALSE, quote = FALSE)
    out <- list()
    tmp <- pathDF
    tmp[, 1] <- sub("_cont.txt", "", tmp[, 1])
    out[["cumulativeFeatScores"]] <- tmp
    scoreColl <- seq_len(maxScore)
    outdf <- matrix(NA, nrow = length(scoreColl), ncol = 9 + 
        4)
    colnames(outdf) <- c("score", "numPathways", "pred_tot", 
        "pred_ol", "pred_pct", "other_tot", "other_ol", "other_pct", 
        "rr", "pred_pct_min", "pred_pct_max", "other_pct_min", 
        "other_pct_max")
    outdf_train <- matrix(NA, nrow = length(scoreColl), ncol = 9 + 
        4)
    colnames(outdf_train) <- colnames(outdf)
    if (enrichLabels) {
        outdf_enriched <- matrix(NA, nrow = length(scoreColl), 
            ncol = 9 + 4)
        colnames(outdf_enriched) <- colnames(outdf)
        outdf_enriched_tr <- matrix(NA, nrow = length(scoreColl), 
            ncol = 9 + 4)
        colnames(outdf_enriched_tr) <- colnames(outdf)
    }
    ctr <- 1
    predContr <- rep("", length(scoreColl))
    otherContr <- rep("", length(scoreColl))
    predContr_cl <- rep("", length(scoreColl))
    otherContr_cl <- rep("", length(scoreColl))
    resampPerf <- list(allNets = list(), enrichedNets = list())
    for (setScore in scoreColl) {
        selPath <- pathDF[which(pathDF[, 2] >= setScore), 1]
        if (verbose) 
            message(sprintf("Thresh = %i ; %i pathways", setScore, 
                length(selPath)))
        currmat <- matrix(NA, nrow = length(TT_STATUS), ncol = 7)
        currmat_enriched <- matrix(NA, nrow = length(TT_STATUS), 
            ncol = 7)
        currmat_train <- matrix(NA, nrow = length(TT_STATUS), 
            ncol = 7)
        currmat_enriched_tr <- matrix(NA, nrow = length(TT_STATUS), 
            ncol = 7)
        predCurr <- ""
        otherCurr <- ""
        predCurr_cl <- ""
        otherCurr_cl <- ""
        for (k in seq_len(length(TT_STATUS))) {
            if (verbose) 
                message(sprintf("\t(k = %i)", k))
            pheno_test <- phenoDF[which(TT_STATUS[[k]] %in% "TEST"), 
                ]
            p_test <- netmat[which(rownames(netmat) %in% pheno_test$ID), 
                ]
            tmp <- updateNets(p_test, pheno_test, writeNewNets = FALSE, 
                verbose = FALSE)
            p_test <- tmp[[1]]
            pheno_test <- tmp[[2]]
            tmp <- getOR(p_test, pheno_test, predClass, selPath, 
                verbose = FALSE)
            predCurr <- c(predCurr, intersect(tmp$OLsamps, pheno_test$ID[which(pheno_test$STATUS %in% 
                predClass)]))
            otherCurr <- c(otherCurr, intersect(tmp$OLsamps, 
                pheno_test$ID[which(!pheno_test$STATUS %in% predClass)]))
            x <- tmp$stats
            currmat[k, ] <- c(x[1, 1], x[1, 2], x[1, 3], x[2, 
                1], x[2, 2], x[2, 3], tmp$relEnr)
            rm(x, tmp, pheno_test, p_test)
            pheno_train <- phenoDF[which(TT_STATUS[[k]] %in% 
                "TRAIN"), ]
            p_train <- netmat[which(rownames(netmat) %in% pheno_train$ID), 
                ]
            tmp <- updateNets(p_train, pheno_train, writeNewNets = FALSE, 
                verbose = FALSE)
            p_train <- tmp[[1]]
            pheno_train <- tmp[[2]]
            tmp <- getOR(p_train, pheno_train, predClass, selPath, 
                verbose = FALSE)
            x <- tmp$stats
            currmat_train[k, ] <- c(x[1, 1], x[1, 2], x[1, 3], 
                x[2, 1], x[2, 2], x[2, 3], tmp$relEnr)
            rm(x, tmp, pheno_train, p_train)
            if (enrichLabels) {
                pheno_test <- phenoDF[which(TT_STATUS[[k]] %in% 
                  "TEST"), ]
                p_test <- netmat[which(rownames(netmat) %in% 
                  pheno_test$ID), ]
                p_test <- p_test[, which(colnames(p_test) %in% 
                  enrichedNets[[k]])]
                tmp <- updateNets(p_test, pheno_test, writeNewNets = FALSE, 
                  verbose = FALSE)
                p_test <- tmp[[1]]
                pheno_test <- tmp[[2]]
                tmp <- getOR(p_test, pheno_test, predClass, selPath, 
                  verbose = FALSE)
                x <- tmp$stats
                currmat_enriched[k, ] <- c(x[1, 1], x[1, 2], 
                  x[1, 3], x[2, 1], x[2, 2], x[2, 3], tmp$relEnr)
                predCurr_cl <- c(predCurr_cl, intersect(tmp$OLsamps, 
                  pheno_test$ID[which(pheno_test$STATUS %in% 
                    predClass)]))
                otherCurr_cl <- c(otherCurr_cl, intersect(tmp$OLsamps, 
                  pheno_test$ID[which(!pheno_test$STATUS %in% 
                    predClass)]))
                pheno_train <- phenoDF[which(TT_STATUS[[k]] %in% 
                  "TRAIN"), ]
                p_train <- netmat[which(rownames(netmat) %in% 
                  pheno_train$ID), ]
                p_train <- p_train[, which(colnames(p_train) %in% 
                  enrichedNets[[k]])]
                tmp <- updateNets(p_train, pheno_train, writeNewNets = FALSE, 
                  verbose = FALSE)
                p_train <- tmp[[1]]
                pheno_train <- tmp[[2]]
                tmp <- getOR(p_train, pheno_train, predClass, 
                  selPath, verbose = FALSE)
                x <- tmp$stats
                currmat_enriched_tr[k, ] <- c(x[1, 1], x[1, 2], 
                  x[1, 3], x[2, 1], x[2, 2], x[2, 3], tmp$relEnr)
            }
        }
        resampPerf[["allNets"]][[setScore]] <- currmat
        resampPerf[["enrichedNets"]][[setScore]] <- currmat_enriched
        predCurr <- unique(predCurr)
        otherCurr <- unique(otherCurr)
        if (verbose) {
            message(sprintf("\t# contrib: %i pred ; %i other", 
                length(predCurr), length(otherCurr)))
        }
        predContr[ctr] <- paste(predCurr, collapse = ",")
        otherContr[ctr] <- paste(otherCurr, collapse = ",")
        if (enrichLabels) {
            predCurr_cl <- unique(predCurr_cl)
            otherCurr_cl <- unique(otherCurr_cl)
            if (verbose) {
                message(paste("\tLABEL ENRICHMENT: ", sprintf("# contributing: %i pred ; %i other", 
                  length(predCurr_cl), length(otherCurr_cl)), 
                  sep = ""))
            }
            predContr_cl[ctr] <- paste(predCurr_cl, collapse = ",")
            otherContr_cl[ctr] <- paste(otherCurr_cl, collapse = ",")
        }
        outdf[ctr, ] <- c(setScore, length(selPath), colMeans(currmat), 
            min(currmat[, 3]), max(currmat[, 3]), min(currmat[, 
                6]), max(currmat[, 6]))
        outdf_train[ctr, ] <- c(setScore, length(selPath), colMeans(currmat_train), 
            min(currmat_train[, 3]), max(currmat_train[, 3]), 
            min(currmat_train[, 6]), max(currmat_train[, 6]))
        if (enrichLabels) {
            outdf_enriched[ctr, ] <- c(setScore, length(selPath), 
                colMeans(currmat_enriched), min(currmat_enriched[, 
                  3]), max(currmat_enriched[, 3]), min(currmat_enriched[, 
                  6]), max(currmat_enriched[, 6]))
            outdf_enriched_tr[ctr, ] <- c(setScore, length(selPath), 
                colMeans(currmat_enriched_tr), min(currmat_enriched_tr[, 
                  3]), max(currmat_enriched_tr[, 3]), min(currmat_enriched_tr[, 
                  6]), max(currmat_enriched_tr[, 6]))
        }
        ctr <- ctr + 1
    }
    numresamp <- nrow(resampPerf[[1]][[1]])
    for (k in seq_len(length(resampPerf))) {
        tmp <- lapply(resampPerf[[k]], function(x) {
            as.numeric(t(x))
        })
        tmp <- do.call(rbind, tmp)
        x <- rep(c("pred_total", "pred_OL", "pred_OL_pct", "other_total", 
            "other_OL", "other_OL_pct", "relEnr"), numresamp)
        colnames(tmp) <- paste(x, rep(seq_len(numresamp), each = 7), 
            sep = "_")
        rownames(tmp) <- scoreColl
        resampPerf[[k]] <- tmp
    }
    out[["resamplingPerformance"]] <- resampPerf
    save(resampPerf, file = paste(outDir, "resamplingPerf.Rdata", 
        sep = "/"))
    outdf <- data.frame(outdf)
    outdf <- cbind(outdf, CONTRIBUT_PRED = predContr, CONTRIBUT_OTHER = otherContr)
    outFile <- paste(outDir, "RR_changeNetSum_stats_denAllNets.txt", 
        sep = "/")
    write.table(outdf, file = outFile, sep = "\t", col.names = TRUE, 
        row.names = FALSE, quote = FALSE)
    out[["performance_denAllNets"]] <- outdf
    write.table(outdf_train, file = outFile, sep = "\t", col.names = TRUE, 
        row.names = FALSE, quote = FALSE)
    out[["performance_denAllNets_TrainingSamples"]] <- outdf_train
    if (enrichLabels) {
        outdf_enriched <- data.frame(outdf_enriched)
        outdf_enriched <- cbind(outdf_enriched, CONTRIBUT_PRED = predContr_cl, 
            CONTRIBUT_OTHER = otherContr_cl)
        outFile <- paste(outDir, "RR_changeNetSum_stats_denEnrichedNets.txt", 
            sep = "/")
        write.table(outdf_enriched, file = outFile, sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
        out[["performance_denEnrichedNets"]] <- outdf_enriched
    }
    return(out)
}
runFeatureSelection <-
function (trainID_pred, outDir, dbPath, numTrainSamps = NULL, 
    incNets = "all", orgName = "predictor", fileSfx = "CV", verbose = FALSE, 
    numCores = 2L, JavaMemory = 6L, verbose_runQuery = FALSE, 
    debugMode = FALSE, ...) 
{
    if (!file.exists(outDir)) 
        dir.create(outDir)
    if (verbose) 
        message("\tWriting queries:\n")
    qSamps <- makeQueries(trainID_pred, verbose = verbose, ...)
    for (m in seq_len(length(qSamps))) {
        qFile <- paste(outDir, sprintf("%s_%i.query", fileSfx, 
            m), sep = "/")
        if (is.null(numTrainSamps)) {
            numTrainSamps = 5
            message("Memory saver option: using 5 training samples for CV")
        }
        writeQueryFile(qSamps[[m]], incNets, numTrainSamps, qFile, 
            orgName)
    }
    qFiles <- list()
    for (m in seq_len(length(qSamps))) {
        qFile <- paste(outDir, sprintf("%s_%i.query", fileSfx, 
            m), sep = "/")
        qFiles <- append(qFiles, qFile)
    }
    runQuery(dbPath, qFiles, outDir, JavaMemory = JavaMemory, 
        verbose = verbose_runQuery, numCores = numCores, debugMode = debugMode)
}
runQuery <-
function (dbPath, queryFiles, resDir, verbose = TRUE, JavaMemory = 6L, 
    numCores = 1L, debugMode = FALSE) 
{
    GM_jar <- getGMjar_path()
    qBase <- basename(queryFiles[[1]][1])
    logFile <- paste(resDir, sprintf("%s.log", qBase))
    queryStrings <- paste(queryFiles, collapse = " ")
    args <- c()
    java_ver <- suppressWarnings(system2("java", args = "--version", 
        stdout = TRUE, stderr = NULL))
    if (any(grep(" 11", java_ver)) || any(grep(" 12", java_ver)) || 
        any(grep(" 13", java_ver)) || any(grep(" 14", java_ver)) || 
        any(grep(" 16", java_ver))) {
        if (verbose) 
            message("Java 11 or later detected")
    }
    else {
        if (verbose) 
            message("Java 8 detected")
        args <- c(args, "-d64")
    }
    args <- c(args, sprintf("-Xmx%iG", JavaMemory * numCores), 
        "-cp", GM_jar)
    args <- c(args, "org.genemania.plugin.apps.QueryRunner")
    args <- c(args, "--data", dbPath, "--in", "flat", "--out", 
        "flat")
    args <- c(args, "--threads", numCores, "--results", resDir, 
        unlist(queryFiles))
    args <- c(args, "--netdx-flag", "true")
    resFile <- paste(resDir, sprintf("%s-results.report.txt", 
        qBase), sep = "/")
    t0 <- Sys.time()
    if (debugMode) {
        message(sprintf("java %s", paste(args, collapse = " ")))
        system2("java", args, wait = TRUE)
    }
    else {
        system2("java", args, wait = TRUE, stdout = NULL, stderr = NULL)
    }
    if (verbose) 
        message(sprintf("QueryRunner time taken: %1.1f s", Sys.time() - 
            t0))
    Sys.sleep(3)
    return(resFile)
}
setupFeatureDB <-
function (pheno, prepDir = tempdir()) 
{
    curd <- getwd()
    setwd(prepDir)
    pheno$INTERNAL_ID <- seq_len(nrow(pheno))
    con <- file("ORGANISMS.txt", "w")
    write(paste("1", "predictor", "my_predictor", "my_predictor", 
        -1, 1339, sep = "\t"), file = con, append = FALSE)
    close(con)
    tmp <- pheno[, c("INTERNAL_ID", "ID", "INTERNAL_ID")]
    tmp$dummy <- 1
    write.table(tmp, file = "NODES.txt", sep = "\t", col.names = FALSE, 
        row.names = FALSE, quote = FALSE)
    tmp <- paste(pheno$INTERNAL_ID, pheno$ID, "N/A", 1, pheno$INTERNAL_ID, 
        1, 0, sep = "\t")
    write.table(tmp, file = "GENES.txt", sep = "\t", col.names = FALSE, 
        row.names = FALSE, quote = FALSE)
    write.table(pheno[, c("INTERNAL_ID", "ID")], file = "GENE_DATA.txt", 
        sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(cbind(pheno$INTERNAL_ID, pheno$INTERNAL_ID), 
        file = "1.synonyms", sep = "\t", col.names = FALSE, row.names = FALSE, 
        quote = FALSE)
    con <- file("GENE_NAMING_SOURCES.txt", "w")
    write("1\tOther\t1\tOther", file = con, append = TRUE)
    write("2\tEntrez Gene ID\t0\tEntrez Gene ID", file = con, 
        append = TRUE)
    close(con)
    file.create("TAGS.txt")
    file.create("NETWORK_TAG_ASSOC.txt")
    file.create("ONTOLOGIES.txt")
    file.create("ONTOLOGY_CATEGORIES.txt")
    con <- file("colours.txt", "w")
    write("geneset_1\tff00ff", file = con, append = TRUE)
    close(con)
    file.create("ATTRIBUTES.txt")
    file.create("ATTRIBUTE_GROUPS.txt")
    con <- file("db.cfg", "w")
    write("[FileLocations]", file = con, append = TRUE)
    write("generic_db_dir = .", file = con, append = TRUE)
    write("[Organisms]", file = con, append = TRUE)
    write("organisms = org_1", file = con, append = TRUE)
    write("[org_1]", file = con, append = TRUE)
    write("gm_organism_id = 1", file = con, append = TRUE)
    write("short_name = predictor", file = con, append = TRUE)
    write("common_name = my_predictor", file = con, append = TRUE)
    write("", file = con, append = TRUE)
    close(con)
    setwd(curd)
    return(pheno)
}
sim.eucscale <-
function (dat, K = 20, alpha = 0.5) 
{
    ztrans <- function(m) {
        m <- as.matrix(m)
        m2 <- apply(m, 1, function(x) {
            (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
        })
        m2
    }
    normalize <- function(X) {
        print(dim(X))
        row.sum.mdiag <- rowSums(X, na.rm = TRUE) - diag(X)
        row.sum.mdiag[row.sum.mdiag == 0] <- 1
        X <- X/(2 * (row.sum.mdiag))
        diag(X) <- 0.5
        X <- (X + t(X))/2
        return(X)
    }
    nnodata <- which(abs(colSums(dat, na.rm = TRUE)) < .Machine$double.eps)
    z1 <- ztrans(dat)
    euc <- as.matrix(dist(z1, method = "euclidean"))^(1/2)
    N <- nrow(euc)
    euc <- (euc + t(euc))/2
    sortedColumns <- as.matrix(t(apply(euc, 2, sort, na.last = TRUE)))
    print(dim(sortedColumns))
    finiteMean <- function(x) {
        return(mean(x[is.finite(x)], na.rm = TRUE))
    }
    means <- apply(sortedColumns[, seq_len(K) + 1], 1, finiteMean)
    means <- means + .Machine$double.eps
    avg <- function(x, y) {
        return((x + y)/2)
    }
    Sig <- outer(means, means, avg)/3 * 2 + euc/3 + .Machine$double.eps
    Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
    densities <- dnorm(euc, 0, alpha * Sig, log = FALSE)
    W <- (densities + t(densities))/2
    W <- normalize(W)
    idx <- which(rowSums(is.na(euc)) == ncol(W) - 1)
    if (any(idx)) {
        W <- W[-idx, ]
        idx <- which(colSums(is.na(euc)) == ncol(W) - 1)
        W <- W[, -idx]
    }
    return(W)
}
sim.pearscale <-
function (dat, K = 20, alpha = 0.5) 
{
    ztrans <- function(m) {
        m <- as.matrix(m)
        m2 <- apply(m, 1, function(x) {
            (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
        })
        m2
    }
    normalize <- function(X) {
        print(dim(X))
        row.sum.mdiag <- rowSums(X, na.rm = TRUE) - diag(X)
        row.sum.mdiag[row.sum.mdiag == 0] <- 1
        X <- X/(2 * (row.sum.mdiag))
        diag(X) <- 0.5
        X <- (X + t(X))/2
        return(X)
    }
    if (nrow(dat) < 6) {
        z1 <- ztrans(dat)
        euc <- as.matrix(dist(z1, method = "euclidean"))^(1/2)
    }
    else {
        euc <- as.matrix(1 - cor(dat, method = "pearson", use = "pairwise.complete.obs"))
    }
    N <- nrow(euc)
    euc <- (euc + t(euc))/2
    sortedColumns <- as.matrix(t(apply(euc, 2, sort, na.last = TRUE)))
    print(dim(sortedColumns))
    finiteMean <- function(x) {
        return(mean(x[is.finite(x)], na.rm = TRUE))
    }
    means <- apply(sortedColumns[, seq_len(K) + 1], 1, finiteMean)
    means <- means + .Machine$double.eps
    avg <- function(x, y) {
        return((x + y)/2)
    }
    Sig <- outer(means, means, avg)/3 * 2 + euc/3 + .Machine$double.eps
    Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
    densities <- dnorm(euc, 0, alpha * Sig, log = FALSE)
    W <- (densities + t(densities))/2
    W <- normalize(W)
    idx <- which(rowSums(is.na(euc)) == ncol(W) - 1)
    if (any(idx)) {
        W <- W[-idx, ]
        idx <- which(colSums(is.na(euc)) == ncol(W) - 1)
        W <- W[, -idx]
    }
    return(W)
}
simpleCap <-
function (x) 
{
    x <- tolower(x)
    s <- strsplit(x, " ")[[1]]
    x <- paste(toupper(substring(s, 1, 1)), substring(s, 2), 
        sep = "", collapse = " ")
    x
}
smoothMutations_LabelProp <-
function (mat, net, numCores = 1L) 
{
    if (class(mat) == "data.frame") 
        mat <- as.matrix(mat)
    if (class(net) == "data.frame") 
        net <- as.matrix(net)
    inds <- split(seq_len(ncol(mat)), sort(rep_len(seq_len(numCores), 
        ncol(mat))))
    res.l <- list()
    required <- c("scater", "clusterExperiment", "netSmooth")
    ctr <- 0
    for (cur in required) {
        if (!requireNamespace(cur, quietly = TRUE)) {
            message(sprintf("Package \"%s\" needed for smoothMutations_LabelProp() to work. Please install it."))
            ctr <- ctr + 1
        }
        if (ctr > 0) 
            stop("Please install needed packages before proceeding.", 
                call. = FALSE)
    }
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    k <- NULL
    res.l <- foreach(k = 1:length(inds), .packages = c("netSmooth", 
        "scater", "clusterExperiment")) %dopar% {
        nS.res = netSmooth::netSmooth(mat[, inds[[k]]], net, 
            alpha = 0.20000000000000001, verbose = "auto", normalizeAdjMatrix = c("columns"))
        return(nS.res)
    }
    stopCluster(cl)
    nS.res <- do.call(cbind, res.l)
    return(nS.res)
}
sparsify2 <-
function (W, outFile = paste(tempdir(), "tmp.txt", sep = "/"), 
    cutoff = 0.29999999999999999, maxInt = 50, EDGE_MAX = 1000, 
    includeAllNodes = TRUE, verbose = TRUE) 
{
    if (verbose) 
        message(sprintf(paste("sparsify2:maxInt=%i;EDGE_MAX=%1.2f;", 
            "cutoff=%1.2e;includeAllNodes=%s", sep = ""), maxInt, 
            EDGE_MAX, cutoff, includeAllNodes))
    if (maxInt > ncol(W)) 
        maxInt <- ncol(W)
    W[upper.tri(W, diag = TRUE)] <- NA
    W[W < cutoff] <- NA
    x <- list()
    for (i in seq_len(nrow(W))) {
        x[[i]] <- sort(W[i, ], decreasing = TRUE, na.last = TRUE)
    }
    names(x) <- rownames(W)
    for (k in seq_len(length(x))) {
        cur <- x[[k]]
        tokeep <- names(cur)[seq_len(min(length(cur), maxInt))]
        W[k, which(!colnames(W) %in% tokeep)] <- NA
    }
    mmat <- na.omit(melt(W))
    mmat <- mmat[order(mmat[, 3], decreasing = TRUE), ]
    if (!is.infinite(EDGE_MAX)) {
        maxEdge <- nrow(mmat)
        if (maxEdge > EDGE_MAX) 
            maxEdge <- EDGE_MAX
        mmat <- mmat[seq_len(maxEdge), ]
    }
    if (includeAllNodes) {
        mmat[, 1] <- as.character(mmat[, 1])
        mmat[, 2] <- as.character(mmat[, 2])
        univ <- c(mmat[, 1], mmat[, 2])
        missing <- setdiff(rownames(W), univ)
        if (length(missing) > 0) {
            message(sprintf(paste("Sparsify2: ", "found %i missing patients; adding strongest edge\n", 
                sep = ""), length(missing)))
            for (k in missing) {
                tmp <- x[[k]]
                if (is.na(tmp[1])) {
                  message(paste("\tMissing edge is below cutoff; ", 
                    "setting to cutoff\n", sep = ""))
                  tmp[1] <- cutoff
                }
                mmat <- rbind(mmat, c(k, names(tmp)[1], tmp[1]))
            }
        }
    }
    head(mmat)
    mmat[, 3] <- as.numeric(mmat[, 3])
    mmat[, 3] <- round(mmat[, 3], digits = 4)
    write.table(mmat, file = outFile, sep = "\t", col.names = FALSE, 
        row.names = FALSE, quote = FALSE)
    return(mmat)
}
sparsify3 <-
function (W, outFile = sprintf("%s/tmp.txt", tempdir()), cutoff = 0.29999999999999999, 
    maxInt = 50, EDGE_MAX = Inf, includeAllNodes = TRUE, verbose = TRUE) 
{
    if (maxInt > ncol(W)) 
        maxInt <- ncol(W)
    if (!is(W, "matrix")) 
        W <- as.matrix(W)
    W[which(is.na(W))] <- .Machine$double.eps
    diag(W) <- NA
    mytop <- cbind(colnames(W), colnames(W)[apply(W, 1, which.max)], 
        apply(W, 1, max, na.rm = TRUE))
    W[upper.tri(W, diag = TRUE)] <- NA
    W[W < cutoff] <- NA
    maxind <- min(ncol(W), maxInt)
    W_order <- t(apply(W, 1, order, decreasing = TRUE, na.last = TRUE))
    W_order[which(W_order > maxInt)] <- NA
    W_order[which(W_order <= maxInt)] <- .Machine$double.eps
    W2 <- W + W_order
    mmat <- na.omit(melt(W2, varnames = names(dimnames(W2))))
    maxEdge <- nrow(mmat)
    if (!is.infinite(EDGE_MAX)) {
        if (maxEdge > EDGE_MAX) 
            maxEdge <- EDGE_MAX
        mmat <- mmat[seq_len(maxEdge), ]
    }
    if (includeAllNodes) {
        mmat[, 1] <- as.character(mmat[, 1])
        mmat[, 2] <- as.character(mmat[, 2])
        univ <- c(mmat[, 1], mmat[, 2])
        missing <- setdiff(rownames(W), univ)
        if (length(missing) > 0) {
            message(sprintf(paste("Sparsify2: found %i missing patients; ", 
                "adding strongest edge\n", sep = ""), length(missing)))
            for (k in missing) {
                tmp <- mytop[which(mytop[, 1] %in% k), ]
                x <- as.numeric(tmp[3])
                if (x < cutoff) {
                  message(paste("\tMissing edge is below cutoff; ", 
                    "setting to cutoff\n"))
                  x <- cutoff
                }
                mmat <- rbind(mmat, c(k, tmp[2], x))
            }
        }
    }
    head(mmat)
    mmat <- na.omit(mmat)
    mmat[, 3] <- as.numeric(mmat[, 3])
    mmat[, 3] <- round(mmat[, 3], digits = 4)
    write.table(mmat, file = outFile, sep = "\t", col.names = FALSE, 
        row.names = FALSE, quote = FALSE)
    return(mmat)
}
splitTestTrain <-
function (pheno_DF, pctT = 0.69999999999999996, verbose = FALSE) 
{
    lvls <- unique(pheno_DF$STATUS)
    IS_TRAIN <- rep("TEST", nrow(pheno_DF))
    for (lv in lvls) {
        idx <- which(pheno_DF$STATUS %in% lv)
        IS_TRAIN[sample(idx, floor(pctT * length(idx)), FALSE)] <- "TRAIN"
    }
    IS_TRAIN <- factor(IS_TRAIN, levels = c("TRAIN", "TEST"))
    pheno_DF <- cbind(pheno_DF, IS_TRAIN = IS_TRAIN)
    if (verbose) 
        print(table(pheno_DF[, c("STATUS", "IS_TRAIN")]))
    return(IS_TRAIN)
}
splitTestTrain_resampling <-
function (pheno_DF, nFold = 3L, predClass, verbose = FALSE) 
{
    plus_idx <- which(pheno_DF$STATUS %in% predClass)
    other_idx <- setdiff(seq_len(nrow(pheno_DF)), plus_idx)
    plus_csize <- floor((1/nFold) * length(plus_idx))
    other_csize <- floor((1/nFold) * length(other_idx))
    plus_tsize <- length(plus_idx) - plus_csize
    other_tsize <- length(other_idx) - other_csize
    if (verbose) {
        message(sprintf("\t(+) %s : %i total ; %i train, %i held-out per\n", 
            predClass, length(plus_idx), plus_tsize, plus_csize))
        message(sprintf("\t(-) (!%s): %i total ; %i train, %i held-out per\n", 
            predClass, length(other_idx), other_tsize, other_csize))
    }
    plus_order <- sample(plus_idx, replace = FALSE)
    other_order <- sample(other_idx, replace = FALSE)
    out <- list()
    for (k in seq_len(nFold)) {
        status <- rep("TRAIN", nrow(pheno_DF))
        sidx <- ((k - 1) * plus_csize) + 1
        eidx <- k * plus_csize
        if (k == nFold) 
            eidx <- length(plus_idx)
        if (verbose) 
            message(sprintf("\t%i (+): %i test (%i-%i);\n", k, 
                eidx - sidx + 1, sidx, eidx))
        status[plus_order[sidx:eidx]] <- "TEST"
        sidx <- ((k - 1) * other_csize) + 1
        eidx <- k * other_csize
        if (k == nFold) 
            eidx <- length(other_idx)
        if (verbose) 
            message(sprintf("\t\t%i (-): %i test\n", k, eidx - 
                sidx + 1))
        status[other_order[sidx:eidx]] <- "TEST"
        out[[k]] <- status
    }
    out
}
subsampleValidationData <-
function (dataMAE, pctValidation = 0.20000000000000001, verbose = TRUE) 
{
    if (pctValidation < 0.050000000000000003 || pctValidation > 
        0.94999999999999996) 
        stop("pctValidation should be between 0.05 and 0.95.")
    if (missing(dataMAE)) 
        stop("Supply dataMAE.")
    if (class(dataMAE) != "MultiAssayExperiment") 
        stop("dataMAE must be an object of type MultiAssayExperiment.")
    pheno <- colData(dataMAE)
    st <- unique(pheno$STATUS)
    nsamp <- round(pctValidation/length(st) * nrow(pheno))
    idx_holdout <- c()
    for (k in unique(pheno$STATUS)) {
        idx_holdout <- c(idx_holdout, sample(which(pheno$STATUS == 
            k), nsamp, FALSE))
    }
    holdout <- dataMAE[, rownames(pheno)[idx_holdout]]
    colData(holdout)$ID <- as.character(colData(holdout)$ID)
    tokeep <- setdiff(1:nrow(pheno), idx_holdout)
    dataMAE <- dataMAE[, rownames(pheno)[tokeep]]
    return(list(trainMAE = dataMAE, validationMAE = holdout))
}
thresholdSmoothedMutations <-
function (smoothedMutProfile, unsmoothedMutProfile, nameDataset, 
    n_topXmuts = c(10)) 
{
    smoothedMutProfile = apply(-smoothedMutProfile, 2, rank)
    n_muts = colSums(unsmoothedMutProfile)
    smoothedMutProfiles_l = list()
    for (k_top in 1:length(n_topXmuts)) {
        name_prop = paste(nameDataset, "_x", n_topXmuts[k_top], 
            sep = "")
        n_new_muts = n_muts * n_topXmuts[k_top]
        for (i_col in 1:length(n_new_muts)) {
            smoothedMutProfile[smoothedMutProfile[, i_col] <= 
                n_new_muts[i_col], i_col] = 1
            smoothedMutProfile[smoothedMutProfile[, i_col] > 
                n_new_muts[i_col], i_col] = 0
        }
        smoothedMutProfiles_l[[name_prop]] = smoothedMutProfile
    }
    if (length(smoothedMutProfiles_l) != 1) {
        return(smoothedMutProfiles_l)
    }
    if (length(smoothedMutProfiles_l) == 1) {
        return(smoothedMutProfile)
    }
}
tSNEPlotter <-
function (psn, pheno, ...) 
{
    message("* Making symmetric matrix")
    symmForm <- suppressMessages(makeSymmetric(psn))
    symmForm[which(is.na(symmForm))] <- .Machine$double.eps
    message("* Running tSNE")
    x <- Rtsne(symmForm, ...)
    dat <- x$Y
    samps <- rownames(symmForm)
    idx <- match(samps, pheno$ID)
    if (all.equal(pheno$ID[idx], samps) != TRUE) {
        stop("pheno IDs not matching psn rownames")
    }
    st <- pheno$STATUS[idx]
    y <- status <- NULL
    message("* Plotting")
    colnames(dat) <- c("x", "y")
    dat <- as.data.frame(dat, stringsAsFactors = TRUE)
    dat$status <- as.factor(st)
    p <- ggplot2::ggplot(dat, aes(x, y)) + geom_point(aes(colour = status))
    p <- p + xlab("") + ylab("") + ggtitle("Integrated PSN - tSNE")
    print(p)
    return(x)
}
updateNets <-
function (p_net, pheno_DF, writeNewNets = TRUE, oldNetDir, newNetDir, 
    verbose = TRUE, ...) 
{
    idx <- which(colSums(p_net) >= 2)
    p_net <- p_net[, idx]
    idx <- which(rowSums(p_net) >= 1)
    p_net <- p_net[idx, ]
    if (verbose) {
        message("Update: (num patients) x (num networks)")
        print(dim(p_net))
    }
    pheno_DF <- pheno_DF[which(pheno_DF$ID %in% rownames(p_net)), 
        ]
    if (writeNewNets) {
        pruneNets(oldNetDir, newNetDir, filterNets = colnames(p_net), 
            filterIDs = rownames(p_net), ...)
    }
    return(list(p_net = p_net, pheno_DF = pheno_DF))
}
writeNetsSIF <-
function (netPath, outFile = paste(tempdir(), "out.sif", sep = "/"), 
    netSfx = "_cont.txt") 
{
    if (.Platform$OS.type == "unix") {
        if (file.exists(outFile)) 
            unlink(outFile)
        file.create(outFile)
    }
    for (n in netPath) {
        netName <- sub(netSfx, "", basename(n))
        message(sprintf("%s\n", netName))
        dat <- read.delim(n, sep = "\t", header = FALSE, as.is = TRUE)
        dat2 <- cbind(dat[, 1], netName, dat[, 2])
        write.table(dat2, file = outFile, append = TRUE, sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
}
writeQueryBatchFile <-
function (netDir, netList, outDir = tempdir(), idFile, orgName = "predictor", 
    orgDesc = "my_predictor", orgAlias = "my_predictor", taxID = 1339) 
{
    outF <- paste(outDir, "batch.txt", sep = "/")
    fileConn <- file(outF, "w")
    tmp <- c("#organism", "id", "file", "name", "description", 
        "alias", "taxonomyid")
    tmp2 <- c("organism", basename(idFile), orgName, orgDesc, 
        orgAlias, as.character(taxID))
    writeLines(sprintf("%s", paste(tmp, collapse = "\t")), con = fileConn)
    writeLines(sprintf("%s\n", paste(tmp2, collapse = "\t")), 
        con = fileConn)
    rm(tmp, tmp2)
    groupName <- "dummy_group"
    groupCode <- "geneset_1"
    groupDesc <- "dummy_group"
    tmp <- c("#group", "name", "code", "description", "RRGGBB colour", 
        "organism")
    tmp2 <- c("group", groupName, groupCode, groupDesc, "ff00ff", 
        orgName)
    writeLines(sprintf("%s", paste(tmp, collapse = "\t")), con = fileConn)
    writeLines(sprintf("%s\n", paste(tmp2, collapse = "\t")), 
        con = fileConn)
    rm(tmp, tmp2)
    tmp <- c("#network", "filename", "name", "description", "group code")
    writeLines(sprintf("%s", paste(tmp, collapse = "\t")), fileConn)
    rm(tmp)
    close(fileConn)
    net_DF <- data.frame(type = "network", filename = netList, 
        name = sub(".txt", "", netList), description = netList, 
        groupCode = groupCode)
    write.table(net_DF, file = outF, sep = "\t", col.names = FALSE, 
        row.names = FALSE, quote = FALSE, append = TRUE)
}
writeQueryFile <-
function (qSamps, incNets = "all", numReturn = 1L, outFile, orgName = "predictor") 
{
    fileConn <- file(outFile, "w")
    writeLines(sprintf("%s", orgName), con = fileConn)
    writeLines(sprintf("%s", paste(qSamps, collapse = "\t")), 
        con = fileConn)
    writeLines(sprintf("%s", paste(incNets, collapse = "\t")), 
        con = fileConn)
    writeLines(sprintf("%i", numReturn), con = fileConn)
    writeLines("automatic", con = fileConn)
    close(fileConn)
}
writeWeightedNets <-
function (patientIDs, netIDs, netDir, keepNets, filterEdgeWt = 0, 
    aggNetFunc = "MAX", limitToTop = 50L, plotEdgeDensity = FALSE, 
    verbose = FALSE) 
{
    pid <- patientIDs
    netid <- netIDs
    aggNetFunc <- toupper(aggNetFunc)
    if (!aggNetFunc %in% c("MEAN", "MAX")) {
        stop("aggNetFunc should be one of: MAX|MEAN\n")
    }
    if (missing(keepNets)) {
        keepNets <- netid$NETWORK
    }
    if (class(keepNets) == "character") {
        keepNets <- data.frame(NETWORK = keepNets, WEIGHT = 1)
        keepNets[, 1] <- as.character(keepNets[, 1])
    }
    simMode <- "normal"
    if (any(colnames(netid) %in% "isBinary")) {
        message("Binary status provided; switching to BinProp mode of similarity!")
        colnames(netid)[3] <- "isBinary"
        simMode <- "BinProp"
    }
    nets <- merge(x = netid, y = keepNets, by = "NETWORK")
    if (simMode == "normal") {
        nets <- nets[, c("NETWORK", "NET_ID", "WEIGHT")]
    }
    else {
        nets <- nets[, c("NETWORK", "NET_ID", "WEIGHT", "isBinary")]
    }
    nets$NET_ID <- as.character(nets$NET_ID)
    x <- sub(".profile$", "", nets$NETWORK)
    x <- sub("_cont.txt", "", x)
    nets$NETWORK_NAME <- x
    numPat <- nrow(pid)
    intColl <- matrix(0, nrow = numPat, ncol = numPat)
    numInt <- matrix(0, nrow = numPat, ncol = numPat)
    contNets <- 1:nrow(nets)
    if (simMode == "BinProp") {
        intColl <- matrix(0, nrow = numPat, ncol = numPat)
        binNets <- which(nets[, "isBinary"] > 0)
        message(sprintf("Got %i binary nets", length(binNets)))
        for (i in binNets) {
            nf <- sprintf(netDir, nets$NET_ID[i], sep = "/")
            ints <- read.delim(nf, sep = "\t", header = F, as.is = T)
            ints <- subset(ints, ints[, 3] >= filterEdgeWt)
            if (nrow(ints) >= 1) {
                midx <- rbind(as.matrix(ints[, c(1:2)]), as.matrix(ints[, 
                  c(2:1)]))
                intColl[midx] <- intColl[midx] + ints[, 3]
                numInt[midx] <- 1
            }
        }
        if (length(binNets) > 0) {
            intColl <- intColl/length(binNets)
            tmp <- qexp(intColl)
            midx <- which(numInt > 0)
            oldVal <- intColl[midx]
            intColl[midx] <- ((tmp[midx]/max(tmp[midx])) * (1 - 
                filterEdgeWt))
            intColl[midx] <- intColl[midx] + filterEdgeWt
        }
        contNets <- setdiff(contNets, which(nets[, "isBinary"] > 
            0))
        message(sprintf("%i continuous nets left", length(contNets)))
    }
    for (i in contNets) {
        nf <- paste(netDir, sprintf("1.%s.txt", nets$NET_ID[i]), 
            sep = "/")
        ints <- read.delim(nf, sep = "\t", header = F, as.is = T)
        oldcount <- nrow(ints)
        ints <- subset(ints, ints[, 3] >= filterEdgeWt)
        if (verbose) {
            message(sprintf("Edge wt filter: %i -> %i interactions", 
                oldcount, nrow(ints)))
        }
        if (nrow(ints) >= 1) {
            midx <- rbind(as.matrix(ints[, c(1:2)]), as.matrix(ints[, 
                c(2:1)]))
            if (aggNetFunc == "MEAN") {
                still_empty <- which(is.na(intColl[midx]))
                if (any(still_empty)) 
                  intColl[midx[still_empty]] <- 0
                intColl[midx] <- intColl[midx] + ints[, 3]
                if (plotEdgeDensity) {
                  tmp <- na.omit(as.numeric(ints[, 3]))
                  hist(tmp, main = nets$NETWORK[i])
                }
                numInt[midx] <- numInt[midx] + 1
            }
            else if (aggNetFunc == "MAX") {
                intColl[midx] <- pmax(intColl[midx], ints[, 3], 
                  na.rm = TRUE)
                if (plotEdgeDensity) {
                  plot(density(na.omit(as.numeric(ints[, 3]))), 
                    main = nets$NETWORK[i])
                }
                numInt[midx] <- numInt[midx] + 1
                if (verbose) 
                  message(sprintf("\t%s: %i: %i interactions added", 
                    basename(nf), i, nrow(midx)))
            }
        }
    }
    if (verbose) 
        message(sprintf("Total of %i nets merged", nrow(nets)))
    intColl[which(numInt < 1)] <- NA
    if (aggNetFunc == "MEAN") 
        tmp <- intColl/numInt
    else tmp <- intColl
    if (!is.infinite(limitToTop)) {
        if (limitToTop >= ncol(intColl)) 
            limitToTop <- Inf
    }
    if (!is.infinite(limitToTop)) {
        message(sprintf("* Limiting to top %i edges per patient", 
            limitToTop))
        for (k in 1:ncol(intColl)) {
            mytop <- order(tmp[k, ], decreasing = TRUE)
            if (limitToTop <= (length(mytop) - 1)) {
                tmp[k, mytop[(limitToTop + 1):length(mytop)]] <- NA
            }
        }
    }
    if (is.infinite(limitToTop)) {
        tmp[lower.tri(tmp, diag = TRUE)] <- NA
    }
    ints <- melt(tmp)
    message(sprintf("\n\t%i pairs have no edges (counts directed edges)", 
        sum(is.nan(ints$value)) + sum(is.na(ints$value))))
    ints <- na.omit(ints)
    if (!is.infinite(limitToTop)) {
        torm <- c()
        n <- nrow(ints)
        for (k in 1:(n - 1)) {
            dup <- which(ints[(k + 1):n, 2] == ints[k, 1] & ints[(k + 
                1):n, 1] == ints[k, 2])
            if (any(dup)) 
                torm <- c(torm, dup)
        }
        message(sprintf("\tRemoving %i duplicate edges", length(torm)))
        if (length(torm) > 0) 
            ints <- ints[-torm, ]
        x <- paste(ints[, 1], ints[, 2], sep = ".")
        y <- paste(ints[, 2], ints[, 1], sep = ".")
        z <- which(y %in% x)
        if (any(z)) {
            ints <- ints[-z, ]
            message(sprintf("\tSecond pass-through: removed %i more dups", 
                length(z)))
        }
        x <- paste(ints[, 1], ints[, 2], sep = ".")
        y <- paste(ints[, 2], ints[, 1], sep = ".")
        dup <- intersect(x, y)
        if (length(dup) > 0) {
            stop("still have duplicates")
        }
    }
    den <- choose(ncol(intColl), 2)
    message(sprintf("\tSparsity = %i/%i (%i %%)", nrow(ints), 
        den, round((nrow(ints)/den) * 100)))
    midx <- match(ints[, 1], pid$GM_ID)
    if (all.equal(pid$GM_ID[midx], ints[, 1]) != TRUE) {
        stop("column 1 doesn't match\n")
    }
    ints$SOURCE <- pid$ID[midx]
    rm(midx)
    midx <- match(ints[, 2], pid$GM_ID)
    if (all.equal(pid$GM_ID[midx], ints[, 2]) != TRUE) {
        stop("column 2 doesn't match\n")
    }
    ints$TARGET <- pid$ID[midx]
    ints <- ints[, c(4, 5, 3, 1, 2)]
    colnames(ints)[1:3] <- c("source", "target", "weights")
    out <- list(filterEdgeWt = filterEdgeWt, aggNetFunc = aggNetFunc, 
        limitToTop = limitToTop, aggNet = ints)
    return(out)
}
checkMakeNetFuncSims <-
function (makeNetFunc, sims, groupList) 
{
    if (is.null(makeNetFunc) && is.null(sims)) {
        stop("Provide either makeNetFunc or sims (preferred).")
    }
    if (!is.null(makeNetFunc) && !is.null(sims)) {
        stop("Provide either makeNetFunc or sims (preferred).")
    }
    if (!is.null(sims)) {
        if (class(sims) != "list") 
            stop("sims must be a list.")
        if (all.equal(sort(names(sims)), sort(names(groupList))) != 
            TRUE) 
            stop("names(sims) must match names(groupList).")
    }
    return(TRUE)
}
psn__corr <-
function (settings, verbose, ...) 
{
    if (verbose) 
        message(sprintf("Layer %s: PEARSON CORR", settings$name))
    nm <- settings$name
    netList <- makePSN_NamedMatrix(xpr = settings$dataList[[nm]], 
        nm = rownames(settings$dataList[[nm]]), namedSets = settings$groupList[[nm]], 
        outDir = settings$netDir, verbose = FALSE, writeProfiles = TRUE, 
        ...)
    return(netList)
}
.get_cache <-
function () 
{
    cache <- rappdirs::user_cache_dir(appname = "netDx")
    BiocFileCache::BiocFileCache(cache, ask = FALSE)
}
buildPredictor_modif <-
function (dataList, groupList, outDir = tempdir(), makeNetFunc = NULL, 
    sims = NULL, featScoreMax = 10L, trainProp = 0.80000000000000004, 
    numSplits = 10L, numCores, JavaMemory = 4L, featSelCutoff = 9L, 
    keepAllData = FALSE, startAt = 1L, preFilter = FALSE, impute = FALSE, 
    preFilterGroups = NULL, imputeGroups = NULL, logging = "default", 
    debugMode = FALSE,seedval=1122007) 
{
    verbose_default <- TRUE
    verbose_runQuery <- FALSE
    verbose_compileNets <- FALSE
    verbose_runFS <- TRUE
    verbose_predict <- FALSE
    verbose_compileFS <- FALSE
    verbose_makeFeatures <- FALSE
    if (logging == "all") {
        verbose_runQuery <- TRUE
        verbose_compileNets <- TRUE
        verbose_compileFS <- TRUE
        verbose_makeFeatures <- TRUE
    }
    else if (logging == "none") {
        verbose_runFS <- FALSE
        verbose_default <- FALSE
        verbose_predict <- FALSE
    }
    if (missing(dataList)) 
        stop("dataList must be supplied.\n")
    if (missing(groupList)) 
        stop("groupList must be supplied.\n")
    if (length(groupList) < 1) 
        stop("groupList must be of length 1+\n")
    tmp <- unlist(lapply(groupList, class))
    not_list <- sum(tmp == "list") < length(tmp)
    nm1 <- setdiff(names(groupList), "clinical")
    if (!is(dataList, "MultiAssayExperiment")) 
        stop("dataList must be a MultiAssayExperiment")
    names_nomatch <- any(!nm1 %in% names(dataList))
    if (!is(groupList, "list") || not_list || names_nomatch) {
        msg <- c("groupList must be a list of lists.", " Names must match those in dataList, and each entry should be a list", 
            " of networks for this group.")
        stop(paste(msg, sep = ""))
    }
    x <- checkMakeNetFuncSims(makeNetFunc = makeNetFunc, sims = sims, 
        groupList = groupList)
    if (!is(dataList, "MultiAssayExperiment")) 
        stop("dataList must be a MultiAssayExperiment")
    if (trainProp <= 0 | trainProp >= 1) 
        stop("trainProp must be greater than 0 and less than 1")
    if (startAt > numSplits) 
        stop("startAt should be between 1 and numSplits")
    megaDir <- outDir
    if (file.exists(megaDir)) {
        stop(paste("outDir seems to already exist!", "Please provide a new directory, as its contents will be overwritten", 
            sprintf("You provided: %s", outDir), sep = "\n"))
    }
    else {
        dir.create(megaDir, recursive = TRUE)
    }
    pheno_all <- colData(dataList)
    pheno_all <- as.data.frame(pheno_all)
    message("Predictor started at:")
    message(Sys.time())
    subtypes <- unique(pheno_all$STATUS)
    exprs <- experiments(dataList)
    datList2 <- list()
    for (k in seq_len(length(exprs))) {
        tmp <- exprs[[k]]
        df <- sampleMap(dataList)[which(sampleMap(dataList)$assay == 
            names(exprs)[k]), ]
        colnames(tmp) <- df$primary[match(df$colname, colnames(tmp))]
        if ("matrix" %in% class(tmp)) {
            datList2[[names(exprs)[k]]] <- tmp
        }
        else {
            tmp <- as.matrix(assays(tmp)[[1]])
            datList2[[names(exprs)[k]]] <- tmp
        }
    }
    if ("clinical" %in% names(groupList)) {
        tmp <- colData(dataList)
        vars <- unique(unlist(groupList[["clinical"]]))
        datList2[["clinical"]] <- t(as.matrix(tmp[, vars, drop = FALSE]))
    }
    dataList <- datList2
    rm(datList2)
    if (verbose_default) {
        message(sprintf("-------------------------------"))
        message(sprintf("# patients = %i", nrow(pheno_all)))
        message(sprintf("# classes = %i { %s }", length(subtypes), 
            paste(subtypes, collapse = ",")))
        message("Sample breakdown by class")
        message(table(pheno_all$STATUS))
        message(sprintf("%i train/test splits", numSplits))
        message(sprintf("Feature selection cutoff = %i of %i", 
            featSelCutoff, featScoreMax))
        message(sprintf("Datapoints:"))
        for (nm in names(dataList)) {
            message(sprintf("\t%s: %i units", nm, nrow(dataList[[nm]])))
        }
    }
    outList <- list()
    tmp <- list()
    for (nm in names(groupList)) {
        curNames <- names(groupList[[nm]])
        tmp[[nm]] <- cbind(rep(nm, length(curNames)), curNames)
    }
    tmp <- do.call("rbind", tmp)
    if (length(nm) < 2) 
        tmp <- as.matrix(tmp)
    colnames(tmp) <- c("NetType", "NetName")
    outList[["inputNets"]] <- tmp
    if (verbose_default) {
        if (!is.null(makeNetFunc)) {
            message("\n\nCustom function to generate input nets:")
            print(makeNetFunc)
        }
        else {
            message("Similarity metrics provided:")
            print(sims)
        }
        message(sprintf("-------------------------------\n"))
    }
    for (rngNum in startAt:numSplits) {
        curList <- list()
        if (verbose_default) {
            message(sprintf("-------------------------------"))
            message(sprintf("Train/test split # %i", rngNum))
            message(sprintf("-------------------------------"))
        }
        outDir <- paste(megaDir, sprintf("rng%i", rngNum), sep = "/")
        dir.create(outDir)
        pheno_all$TT_STATUS <- splitTestTrain2(pheno_all, pctT = trainProp, 
            verbose = verbose_default,rngNum,numSplits,seedval=seedval)
        pheno <- pheno_all[which(pheno_all$TT_STATUS %in% "TRAIN"), 
            ]
        dats_train <- lapply(dataList, function(x) x[, which(colnames(x) %in% 
            pheno$ID), drop = FALSE])
        if (impute) {
            if (verbose_default) 
                message("**** IMPUTING ****")
            if (is.null(imputeGroups)) 
                imputeGroups <- names(dats_train)
            if (!any(imputeGroups %in% names(dats_train))) 
                stop("imputeGroups must match names in dataList")
            nmset <- names(dats_train)
            dats_train <- lapply(names(dats_train), function(nm) {
                x <- dats_train[[nm]]
                print(class(x))
                if (nm %in% imputeGroups) {
                  missidx <- which(rowSums(is.na(x)) > 0)
                  for (i in missidx) {
                    na_idx <- which(is.na(x[i, ]))
                    x[i, na_idx] <- median(x[i, ], na.rm = TRUE)
                  }
                }
                x
            })
            names(dats_train) <- nmset
        }
        if (preFilter) {
            if (is.null(preFilterGroups)) 
                preFilterGroups <- names(dats_train)
            if (!any(preFilterGroups %in% names(dats_train))) {
                stop("preFilterGroups must match names in dataList")
            }
            message("Prefiltering enabled")
            for (nm in preFilterGroups) {
                message(sprintf("%s: %i variables", nm, nrow(dats_train[[nm]])))
                if (nrow(dats_train[[nm]]) < 2) 
                  vars <- rownames(dats_train[[nm]])
                else {
                  newx <- na.omit(dats_train[[nm]])
                  tmp <- pheno[which(pheno$ID %in% colnames(newx)), 
                    ]
                  tryCatch({
                    fit <- cv.glmnet(x = t(newx), y = factor(tmp$STATUS), 
                      family = "binomial", alpha = 1)
                  }, error = function(ex) {
                    print(ex)
                    message("*** You may need to set impute=TRUE for prefiltering ***")
                  }, finally = {
                  })
                  wt <- abs(coef(fit, s = "lambda.min")[, 1])
                  vars <- setdiff(names(wt)[which(wt > .Machine$double.eps)], 
                    "(Intercept)")
                }
                if (length(vars) > 0) {
                  tmp <- dats_train[[nm]]
                  tmp <- tmp[which(rownames(tmp) %in% vars), 
                    , drop = FALSE]
                  dats_train[[nm]] <- tmp
                }
                else {
                }
                message(sprintf("rngNum %i: %s: %s pruned", rngNum, 
                  nm, length(vars)))
            }
        }
        if (verbose_default) {
            message("# values per feature (training)")
            for (nm in names(dats_train)) {
                message(sprintf("\tGroup %s: %i values", nm, 
                  nrow(dats_train[[nm]])))
            }
        }
        netDir <- paste(outDir, "tmp", sep = "/")
        dir.create(netDir)
        pheno_id <- setupFeatureDB(pheno, netDir)
        if (verbose_default) 
            message("** Creating features")
        createPSN_MultiData(dataList = dats_train, groupList = groupList, 
            pheno = pheno_id, netDir = netDir, makeNetFunc = makeNetFunc, 
            sims = sims, numCores = numCores, verbose = verbose_makeFeatures)
        if (verbose_default) 
            message("** Compiling features")
        dbDir <- compileFeatures(netDir, outDir, numCores = numCores, 
            verbose = verbose_compileFS, debugMode = debugMode)
        if (verbose_default) 
            message("\n** Running feature selection")
        curList[["featureScores"]] <- list()
        for (g in subtypes) {
            pDir <- paste(outDir, g, sep = "/")
            if (file.exists(pDir)) 
                unlink(pDir, recursive = TRUE)
            dir.create(pDir)
            if (verbose_default) 
                message(sprintf("\tClass: %s", g))
            pheno_subtype <- pheno
            pheno_subtype$STATUS[which(!pheno_subtype$STATUS %in% 
                g)] <- "nonpred"
            trainPred <- pheno_subtype$ID[which(pheno_subtype$STATUS %in% 
                g)]
            if (verbose_default) {
                print(table(pheno_subtype$STATUS, useNA = "always"))
            }
            resDir <- paste(pDir, "GM_results", sep = "/")
            message(sprintf("\tScoring features"))
            runFeatureSelection(trainPred, outDir = resDir, dbPath = dbDir$dbDir, 
                nrow(pheno_subtype), verbose = verbose_runFS, 
                numCores = numCores, verbose_runQuery = TRUE, 
                featScoreMax = featScoreMax, JavaMemory = JavaMemory, 
                debugMode = debugMode)
            tmp <- dir(path = resDir, pattern = "RANK$")[1]
            tmp <- sprintf("%s/%s", resDir, tmp)
            if (sum(grepl(pattern = ",", readLines(tmp, n = 6)) > 
                0)) {
                replacePattern(path = resDir, fileType = "RANK$")
            }
            nrank <- dir(path = resDir, pattern = "NRANK$")
            if (verbose_default) 
                message("\tCompiling feature scores")
            pTally <- compileFeatureScores(paste(resDir, nrank, 
                sep = "/"), verbose = verbose_compileFS)
            tallyFile <- paste(resDir, sprintf("%s_pathway_CV_score.txt", 
                g), sep = "/")
            write.table(pTally, file = tallyFile, sep = "\t", 
                col.names = TRUE, row.names = FALSE, quote = FALSE)
            curList[["featureScores"]][[g]] <- pTally
            if (verbose_default) 
                message("")
        }
        if (verbose_default) 
            message("\n** Predicting labels for test")
        pheno <- pheno_all
        predRes <- list()
        curList[["featureSelected"]] <- list()
        for (g in subtypes) {
            if (verbose_default) 
                message(sprintf("%s", g))
            pDir <- paste(outDir, g, sep = "/")
            pTally <- read.delim(paste(pDir, "GM_results", sprintf("%s_pathway_CV_score.txt", 
                g), sep = "/"), sep = "\t", header = TRUE, 
                as.is = TRUE)
            idx <- which(pTally[, 2] >= featSelCutoff)
            pTally <- pTally[idx, 1]
            pTally <- sub(".profile", "", pTally)
            pTally <- sub("_cont.txt", "", pTally)
            curList[["featureSelected"]][[g]] <- pTally
            if (verbose_default) 
                message(sprintf("\t%i feature(s) selected", length(pTally)))
            netDir <- paste(pDir, "networks", sep = "/")
            dats_tmp <- list()
            for (nm in names(dataList)) {
                passed <- rownames(dats_train[[nm]])
                tmp <- dataList[[nm]]
                dats_tmp[[nm]] <- tmp[which(rownames(tmp) %in% 
                  passed), , drop = FALSE]
            }
            if (impute) {
                train_samp <- pheno_all$ID[which(pheno_all$TT_STATUS %in% 
                  "TRAIN")]
                test_samp <- pheno_all$ID[which(pheno_all$TT_STATUS %in% 
                  "TEST")]
                nmSet <- names(dats_tmp)
                dats_tmp <- lapply(names(dats_tmp), function(nm) {
                  x <- dats_tmp[[nm]]
                  if (nm %in% imputeGroups) {
                    missidx <- which(rowSums(is.na(x)) > 0)
                    train_idx <- which(colnames(x) %in% train_samp)
                    test_idx <- which(colnames(x) %in% test_samp)
                    for (i in missidx) {
                      na_idx <- intersect(which(is.na(x[i, ])), 
                        train_idx)
                      na_idx1 <- na_idx
                      x[i, na_idx] <- median(x[i, train_idx], 
                        na.rm = TRUE)
                      na_idx <- intersect(which(is.na(x[i, ])), 
                        test_idx)
                      na_idx2 <- na_idx
                      x[i, na_idx] <- median(x[i, test_idx], 
                        na.rm = TRUE)
                    }
                  }
                  x
                })
                names(dats_tmp) <- nmSet
            }
            if (verbose_default) 
                message(sprintf("\t%s: Create & compile features", 
                  g))
            if (length(pTally) >= 1) {
                netDir <- paste(pDir, "tmp", sep = "/")
                dir.create(netDir)
                pheno_id <- setupFeatureDB(pheno, netDir)
                createPSN_MultiData(dataList = dats_tmp, groupList = groupList, 
                  pheno = pheno_id, netDir = netDir, makeNetFunc = makeNetFunc, 
                  sims = sims, numCores = numCores, filterSet = pTally, 
                  verbose = verbose_default)
                dbDir <- compileFeatures(netDir, outDir = pDir, 
                  numCores = numCores, verbose = verbose_compileNets, 
                  debugMode = debugMode)
                qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
                  pheno$TT_STATUS %in% "TRAIN")]
                qFile <- paste(pDir, sprintf("%s_query", g), 
                  sep = "/")
                writeQueryFile(qSamps, "all", nrow(pheno), qFile)
                if (verbose_default) 
                  message(sprintf("\t** %s: Compute similarity", 
                    g))
                resFile <- runQuery(dbDir$dbDir, qFile, resDir = pDir, 
                  JavaMemory = JavaMemory, numCores = numCores, 
                  verbose = verbose_runQuery, debugMode = debugMode)
                predRes[[g]] <- getPatientRankings(sprintf("%s.PRANK", 
                  resFile), pheno, g)
            }
            else {
                predRes[[g]] <- NA
            }
        }
        if (verbose_default) 
            message("")
        if (sum(is.na(predRes)) > 0 & verbose_default) {
            str <- sprintf("RNG %i : One or more classes have no selected features.", 
                rngNum)
            str <- sprintf("%s Not classifying.", str)
            message(str)
        }
        else {
            if (verbose_default) 
                message("** Predict labels")
            predClass <- predictPatientLabels(predRes, verbose = verbose_predict)
            out <- merge(x = pheno_all, y = predClass, by = "ID")
            outFile <- paste(outDir, "predictionResults.txt", 
                sep = "/")
            acc <- sum(out$STATUS == out$PRED_CLASS)/nrow(out)
            if (verbose_default) 
                message(sprintf("Split %i: ACCURACY (N=%i test) = %2.1f%%", 
                  rngNum, nrow(out), acc * 100))
            curList[["predictions"]] <- out
            curList[["accuracy"]] <- acc
        }
        if (!keepAllData) {
            unlink(outDir, recursive = TRUE)
        }
        if (verbose_default) {
            message("\n----------------------------------------")
        }
        outList[[sprintf("Split%i", rngNum)]] <- curList
    }
    message("Predictor completed at:")
    message(Sys.time())
    return(outList)
}
splitTestTrain3<- function (pheno_DF, pctT = 0.69999999999999996, verbose = FALSE, rngNum, numSplits)
{
  set.seed(1122007)
  nsamples<-length(pheno_DF$STATUS)
  permuted = sample(1:nsamples)
  
  indexStart = nsamples/numSplits*(rngNum-1)+1
  indexEnd = nsamples/numSplits*rngNum
  indextest = indexStart:indexEnd
  #Uindextrain = setdiff(1:nsamples,indextest)
  test = permuted[indextest]
  #train = permuted[indextrain]
  
  IS_TRAIN <- rep("TRAIN", nrow(pheno_DF))
  IS_TRAIN[test]<-"TEST"
  
  IS_TRAIN <- factor(IS_TRAIN, levels = c("TRAIN", "TEST"))
  pheno_DF <- cbind(pheno_DF, IS_TRAIN = IS_TRAIN)
  if (verbose)
    print(table(pheno_DF[, c("STATUS", "IS_TRAIN")]))
  return(IS_TRAIN)
}
splitTestTrain2<- function (pheno_DF, pctT = 0.69999999999999996, verbose = FALSE, rngNum, numSplits,seedval=1122007)
{
  set.seed(seedval)
  nsamples<-length(pheno_DF$STATUS)
  permuted = sample(rep(1:numSplits,length.out=nsamples))
  
  IS_TRAIN <- rep("TRAIN", nsamples)
  IS_TRAIN[permuted==rngNum]<-"TEST"
  
  IS_TRAIN <- factor(IS_TRAIN, levels = c("TRAIN", "TEST"))
  pheno_DF <- cbind(pheno_DF, IS_TRAIN = IS_TRAIN)
  if (verbose)
    print(table(pheno_DF[, c("STATUS", "IS_TRAIN")]))
  return(IS_TRAIN)
}
