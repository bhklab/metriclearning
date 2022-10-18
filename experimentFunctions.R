library(CMAPToolkit)
library(cmapR)
library(torch)

# Analyze L1000 datasets
analyzeL1KData <- function(dspath, metapath, cell_id, outpath=".", method="xval", epochs=10, saveModel=TRUE){
  l1k_meta <- read_l1k_meta(metapath, version=2020)
  
  if (cell_id == "all"){
    mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$pert_type == "trt_cp",]
    # This is a bit of a hack to get the learning approach to consider only same-cell line pairs
    mysigs$pert_iname <- sprintf("%s_%s", mysigs$pert_iname, mysigs$cell_id)
    
    # Mysigs has 196k unique pert+cell combinations. For expediency, sample:
    mysigs <- mysigs[mysigs$pert_iname %in% sample(unique(mysigs$pert_iname), 20000), ]

    ds <- parse_gctx(fname=get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=l1k_meta$siginfo$sig_id[l1k_meta$siginfo$pert_type == "trt_cp"])
    ds <- subset_gct(ds, cid=mysigs$sig_id)
  } else {
    mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$cell_id == cell_id & l1k_meta$siginfo$pert_type == "trt_cp",]
    ds <- parse_gctx(fname=get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  }
  
  print(sprintf("DS dimensions: %d x %d", dim(ds@mat)[1], dim(ds@mat)[2]))
  
  if (method == "xval"){
    nFolds <- 5
    
    L1Kxval <- metricCrossValidate(t(ds@mat), mysigs$pert_iname, nFolds=nFolds, epochs=epochs)

    # Save the result
    saveRDS(L1Kxval$res_all, file=file.path(outpath, sprintf("L1Kxval_epch=%s_folds=%s_cell=%s.rds", epochs, nFolds, cell_id)))
    
    if(saveModel){
      for (ii in seq_along(L1Kxval$mymodels)){
        torch_save(L1Kxval$mymodel[[ii]], file.path(outpath, sprintf("L1Kxval_epch=%s_folds=%s_cell=%s_model%d.pt", epochs, nFolds, cell_id, ii)))
      }
    }
    
    return(L1Kxval)
    
  } else if (method == "allds"){
    L1Kmetric <- learnInnerProduct(t(ds@mat), mysigs$pert_iname, epochs=epochs)
    saveRDS(L1Kmetric, file=file.path(outpath, sprintf("L1Kmetric_epch=%s_cell=%s.rds", epochs, cell_id)))
    torch_save(L1Kmetric$model, file.path(outpath, sprintf("L1Kmetric_epch=%s_cell=%s_model.pt", epochs, cell_id)))
    return(L1Kmetric)
    
  } else if (method == "ntraining"){
    
    res <- metricNTraining(t(ds@mat), mysigs$pert_iname, epochs=epochs, reps=5)
    saveRDS(res, file=file.path(outpath, sprintf("L1Kntraining_epch=%s_cell=%s.rds", epochs, cell_id)))
    return(res)
  }
}


# Analyze L1000 cell line specificity
analyzeL1KCellLineSpecificity <- function(dspath, metapath, cell_ids=c("HEPG2", "HCC515", "NPC", "ASC", "HEK293", "YAPC"), outpath=".", ncpds=1000, epochs=10, iter=1){
  l1k_meta <- read_l1k_meta(metapath, version=2020)
  
  mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$pert_type == "trt_cp" & l1k_meta$siginfo$cell_id %in% cell_ids,]
  mysigs$allpert_iname <- sprintf("%s_%s", mysigs$pert_iname, mysigs$cell_id)

  
  res <- data.frame(iter=numeric(), cell_id=character(), trainSet=character(), trainLoss=numeric(), 
                    testLoss=numeric(), testAvgLoss=numeric(), testAURank=numeric())
  
  ds <- parse_gctx(get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  
  for (ii in seq_len(iter)){
    print(sprintf("ii = %s", ii))
    grpPerts <- sample(unique(mysigs$pert_iname), min(length(unique(mysigs$pert_iname)), ncpds))
    trainPerts <- sample(grpPerts, round(length(grpPerts)/2))
    testPerts <- setdiff(grpPerts, trainPerts)
    
    for (mycell in cell_ids){
      print(sprintf("mycell = %s", mycell))
      all_train <- subset_gct(ds, cid=mysigs$sig_id[mysigs$cell_id != mycell & mysigs$pert_iname %in% trainPerts])
      
      ds_train <- subset_gct(ds, cid=mysigs$sig_id[mysigs$cell_id == mycell & mysigs$pert_iname %in% trainPerts])
      ds_test <- subset_gct(ds, cid=mysigs$sig_id[mysigs$cell_id == mycell & mysigs$pert_iname %in% testPerts])
      
      cellModel <- learnInnerProduct(t(ds_train@mat), mysigs$pert_iname[match(ds_train@cid, mysigs$sig_id)], epochs=epochs)
      allModel <- learnInnerProduct(t(all_train@mat), mysigs$allpert_iname[match(all_train@cid, mysigs$sig_id)], epochs=epochs)
      
      cosineTrainSim <- innerProductGroups("cosine", t(ds_train@mat), mysigs$pert_iname[match(ds_train@cid, mysigs$sig_id)], compact=1)
      cosineGrpSim <- innerProductGroups("cosine", t(ds_test@mat), mysigs$pert_iname[match(ds_test@cid, mysigs$sig_id)], compact=1)
      cellGrpSim <- innerProductGroups(cellModel$model, t(ds_test@mat), mysigs$pert_iname[match(ds_test@cid, mysigs$sig_id)], compact=1)
      allGrpSim <- innerProductGroups(allModel$model, t(ds_test@mat), mysigs$pert_iname[match(ds_test@cid, mysigs$sig_id)], compact=1)
      
      res <- rbind(res, data.frame(iter=ii, 
                                   cell_id=mycell,
                                   trainSet="cosine",
                                   trainAvgLoss=mean((mean(cosineTrainSim$diff) - sapply(cosineTrainSim$same, mean))/sd(cosineTrainSim$diff)),
                                   testLoss=(mean(cosineGrpSim$diff - mean(unlist(cosineGrpSim$same)))/sd(cosineGrpSim$diff)),
                                   testAvgLoss=mean((mean(cosineGrpSim$diff) - sapply(cosineGrpSim$same, mean))/sd(cosineGrpSim$diff)), 
                                   testAURank=1 - mean(rankVectors(unlist(cosineGrpSim$same), cosineGrpSim$diff))))
      
      res <- rbind(res, data.frame(iter=ii, 
                                   cell_id=mycell,
                                   trainSet="cell",
                                   trainAvgLoss=cellModel$mean_train_losses[length(cellModel$mean_train_losses)],
                                   testLoss=(mean(cellGrpSim$diff - mean(unlist(cellGrpSim$same)))/sd(cellGrpSim$diff)),
                                   testAvgLoss=mean((mean(cellGrpSim$diff) - sapply(cellGrpSim$same, mean))/sd(cellGrpSim$diff)), 
                                   testAURank=1 - mean(rankVectors(unlist(cellGrpSim$same), cellGrpSim$diff))))
      
      res <- rbind(res, data.frame(iter=ii, 
                                   cell_id=mycell,
                                   trainSet="all",
                                   trainAvgLoss=allModel$mean_train_losses[length(allModel$mean_train_losses)],
                                   testLoss=(mean(allGrpSim$diff - mean(unlist(allGrpSim$same)))/sd(allGrpSim$diff)),
                                   testAvgLoss=mean((mean(allGrpSim$diff) - sapply(allGrpSim$same, mean))/sd(allGrpSim$diff)), 
                                   testAURank=1 - mean(rankVectors(unlist(allGrpSim$same), allGrpSim$diff))))
    }
  }
  
  saveRDS(res, file=file.path(outpath, sprintf("L1KCellLineSpec_nLines=%d_epoch=%s_npcds=%s.rds", length(cell_ids), epochs, ncpds)))
  return(res)
}


# Apply models to replicate and PCL prediction
benchmarkL1KCellModel <- function(modelpath, dspath, metapath, pclpath, cell_id, outpath=".", modelname="mymodel", mkplots=1){
  pclds <- readRDS(pclpath)
  model <- torch_load(modelpath)
  pertnames <- unlist(pclds$pertnames)
  
  # Read appropriate data
  # Refactoring, I could make this a function like loadL1KData:
  l1k_meta <- read_l1k_meta(metapath, version=2020)
  
  if (cell_id == "all"){
    mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$pert_type == "trt_cp",]
    # This is a bit of a hack to get the learning approach to consider only same-cell line pairs
    mysigs <- mysigs[sample(seq_len(dim(mysigs)[1]), 15000)]
    ds <- parse_gctx(fname=get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=l1k_meta$siginfo$sig_id[l1k_meta$siginfo$pert_type == "trt_cp"])
    ds <- subset_gct(ds, cid=mysigs$sig_id)
  } else {
    mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$cell_id == cell_id & l1k_meta$siginfo$pert_type == "trt_cp",]
    
    if ((dim(mysigs)[1]) > 10000){
      pclCount <- sum(mysigs$pert_iname %in% pertnames)
      ax <- which(mysigs$pert_iname %in% pertnames)
      bx <- sample(setdiff(seq_len(dim(mysigs)[1]), ax), pclCount)
      mysigs <- rbind(mysigs[ax, ], mysigs[bx,])
    }
    ds <- parse_gctx(fname=get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  }
  
  cossim <- innerProduct("cosine", t(ds@mat), t(ds@mat))
  MLsim <- innerProduct(model, t(ds@mat), t(ds@mat))
  
  # For each PCL, compute the similarity within different compounds in each PCL
  
  MLSims <- list()
  cosSims <- list()
  
  for (ii in seq_along(pclds$pertnames)){
    pclcpds <- pclds$pertnames[[ii]]
    ix <- which(mysigs$pert_iname %in% pclcpds)
    
    # Identify PCL members that are different compounds
    filt <- outer(mysigs$pert_iname[ix], mysigs$pert_iname[ix], '!=')
    filt <- filt[upper.tri(filt)]
    
    if (sum(filt) > 1){
      pclMLSim <- MLsim[ix,ix]
      pclCosSim  <- cossim[ix,ix]
      MLSims[[pclds$pcldf$pclid[ii]]] <- (pclMLSim[upper.tri(pclMLSim)])[filt]
      cosSims[[pclds$pcldf$pclid[ii]]] <- (pclCosSim[upper.tri(pclCosSim)])[filt]
    }
  }
  
  cosbg <- sample(cossim[upper.tri(cossim)], 1e6)
  MLbg <- sample(MLsim[upper.tri(MLsim)], 1e6)
  
  MLranks <- listify(MLSims, rankVectors(unlist(MLSims), MLbg))
  cosranks <- listify(cosSims, rankVectors(unlist(cosSims), cosbg))

  retdf <- list(modelname=modelname, MLSims=MLSims, cosSims=cosSims, MLRanks=MLranks, cosranks=cosranks, MLbg=MLbg, cosbg=cosbg)
  saveRDS(retdf, file.path(outpath, sprintf("%s_pcl_res.rds", modelname)))
  
  if (mkplots){
    pdf(file=file.path(outpath, sprintf("%s_pcl_figs.pdf", modelname)), width=8, height=6)
    # ML density
    plot(density(MLsim[upper.tri(MLsim)], bw=0.01), col="blue", xlim=c(-1,1), xlab="Metric Learning similarity", lwd=2, main=modelname)
    lines(density(unlist(MLSims), bw=0.01), col="red", lwd=2)
    legend(x="topright", legend=c("PCL pairs", "All pairs"), lwd=c(4,4), col=c("red", "blue"))
    
    # Cos Density
    plot(density(cossim[upper.tri(cossim)], bw=0.01), col="purple", xlim=c(-1,1), xlab="Cosine similarity", lwd=2, main=modelname)
    lines(density(unlist(cosSims), bw=0.01), col="forestgreen", lwd=2)
    legend(x="topright", legend=c("PCL pairs", "All pairs"), lwd=c(4,4), col=c("purple", "forestgreen"))
    
    # Rank plot - unbalanced
    plot(ecdf(unlist(MLranks)), col="red", lwd=2, xlim=c(0,1), xlab="Rank relative to background", ylab="Cumulative Density", main=modelname)
    lines(ecdf(unlist(cosranks)), col="forestgreen", lwd=2)
    lines(c(0,1), c(0,1), col="black", lwd=1, lty=2)
    legend(x="bottomright", legend=c(sprintf("Metric Learning: %0.3f", 1-mean(unlist(MLranks))), 
                                     sprintf("Cosine: %0.3f", 1-mean(unlist(cosranks)))), col=c("red", "forestgreen"), lwd=c(4,4))
    
    # Rank plot - balanced
    MLBalRanks <- unlist(sapply(MLranks, FUN=function(x) sample(x, min(200, length(x)))))
    cosBalRanks <- unlist(sapply(cosranks, FUN=function(x) sample(x, min(200, length(x)))))
    
    plot(ecdf(MLBalRanks), col="red", lwd=2, xlim=c(0,1), 
         xlab="Balanced Rank relative to background", ylab="Cumulative Density", main=modelname)
    lines(ecdf(cosBalRanks), col="forestgreen", lwd=2)
    lines(c(0,1), c(0,1), col="black", lwd=1, lty=2)
    legend(x="bottomright", legend=c(sprintf("Metric Learning: %0.3f", 1-mean(unlist(MLBalRanks))), 
                                     sprintf("Cosine: %0.3f", 1-mean(unlist(cosBalRanks)))), col=c("red", "forestgreen"), lwd=c(4,4))
    
    # Maybe add Z-scores
    dev.off()

  }
  
  return(retdf)
  print("End of benchmark")
}


# Analyze Bray dataset
analyzeBrayData <- function(braypath, outpath=".", method="xval", epochs=10, saveModel=TRUE, subsample=0.2){
  epochs <- as.numeric(epochs)
  subsample <- as.numeric(subsample)
  
  brayds <- loadBrayData(braypath, center=1, subsample = subsample)
  brayData <- brayds$ds
  brayMeta <- brayds$metads
  
  if (method == "xval"){
    nFolds <- 3
    brayxval <- metricCrossValidate(brayData, brayMeta$Metadata_pert_id, nFolds=nFolds, epochs=epochs) 
    saveRDS(brayxval$res_all, file=file.path(outpath, sprintf("brayxval_epch=%d_smp=%d_folds=%d.rds", epochs, round(100*subsample), nFolds)))
    
    if (saveModel){
      for (ii in seq_along(brayxval$mymodels)){
        torch_save(brayxval$mymodel[[ii]], file.path(outpath, sprintf("brayxval_epch=%d_smp=%d_folds=%d_model%d.pt", epochs, round(100*subsample), nFolds, ii)))
      }
    }
    return(brayxval)
    
  } else if (method == "allds"){
    braymetric <- learnInnerProduct(brayData, brayMeta$Metadata_pert_id, epochs=epochs)
    saveRDS(braymetric, file=file.path(outpath, sprintf("braymetric_epch=%d_smp=%d.rds", epochs, round(100*subsample))))
    torch_save(braymetric$model, file.path(outpath, sprintf("braymetric_epch=%d_smp=%d_model.pt", epochs, round(100*subsample))))
    return(braymetric)
    
  } else if (method == "ntraining"){
    # Study how decreasing the amount of training data affects outcomes
    
    res <- metricNTraining(brayData, brayMeta$Metadata_pert_id, epochs=epochs, reps=2)
    saveRDS(res, file=file.path(outpath, sprintf("brayntraining_epch=%s_smp=%d.rds", epochs, round(100*subsample))))
    return(res)
  }
}

analyzeLincsCPData <- function(lincspath, outpath=".", method="xval", epochs=10, saveModel=TRUE, dsname="dsA"){
  epochs <- as.numeric(epochs)
  subsample <- as.numeric(subsample)
  
  lincsds <- loadLincsData(lincspath, splitGrps = 1)
  lincsData <- lincsds$ds
  lincsMeta <- lincsds$metads
  
  if (method == "xval"){
    nFolds <- 3
    lincsxval <- metricCrossValidate(lincsData, lincsMeta$Metadata_pert_iname, nFolds=nFolds, epochs=epochs) 
    saveRDS(lincsxval$res_all, file=file.path(outpath, sprintf("lincsxval_%s_epch=%d_smp=%d_folds=%d.rds", dsname, epochs, round(100*subsample), nFolds)))
    
    if (saveModel){
      for (ii in seq_along(lincsxval$mymodels)){
        torch_save(lincsxval$mymodel[[ii]], file.path(outpath, sprintf("lincsxval_%s_epch=%d_smp=%d_folds=%d_model%d.pt", dsname, epochs, round(100*subsample), nFolds, ii)))
      }
    }
    return(lincsxval)
    
  } else if (method == "allds"){
    lincsmetric <- learnInnerProduct(lincsData, lincsMeta$Metadata_pert_iname, epochs=epochs)
    saveRDS(lincsmetric, file=file.path(outpath, sprintf("lincsmetric_%s_epch=%d.rds", dsname, epochs)))
    torch_save(lincsmetric$model, file.path(outpath, sprintf("lincsmetric_%s_epch=%d_smp=%d_model.pt", dsname, epochs, round(100*subsample))))
    return(lincsmetric)
    
  } else if (method == "ntraining"){
    # Study how decreasing the amount of training data affects outcomes
    
    res <- metricNTraining(lincsData, lincsMeta$Metadata_pert_iname, epochs=epochs, reps=2)
    saveRDS(res, file=file.path(outpath, sprintf("lincsntraining_%s_epch=%s_smp=%d.rds", dsname, epochs, round(100*subsample))))
    return(res)
  }
}



# analyze generic cell painting dataset
analyzeNormCPData <- function(cppath, outpath=".", method="xval", epochs=10, dsname="myds", saveModel=FALSE){
  epochs <- as.numeric(epochs)
  
  cpds <- loadLincsData(cppath, splitGrps = 1)
  cpData <- cpds$ds
  cpMeta <- cpds$metads
  
  # Safety check
  cpData[is.na(cpData)] <- 0
  
  if (method == "xval"){
    nFolds <- 3
    cpxval <- metricCrossValidate(cpData, cpMeta$Metadata_pert_id, nFolds=nFolds, epochs=epochs) 
    saveRDS(cpxval$res_all, file=file.path(outpath, sprintf("%scp_xval_epch=%d_folds=%d.rds", dsname, epochs, nFolds)))
    
    if (saveModel){
      for (ii in seq_along(cpxval$mymodels)){
        torch_save(cpxval$mymodel[[ii]], file.path(outpath, sprintf("%scp_xval_epch=%d_folds=%d_model%d.pt", dsname, epochs, nFolds, ii)))
      }
    }
    return(cpxval)
    
  } else if (method == "allds"){
    cpmetric <- learnInnerProduct(cpData, cpMeta$Metadata_pert_id, epochs=epochs)
    saveRDS(cpmetric, file=file.path(outpath, sprintf("%scp_metric_epch=%d.rds", dsname, epochs)))
    torch_save(cpmetric$model, file.path(outpath, sprintf("%scp_metric_epch=%d_model.pt", dsname, epochs)))
    return(cpmetric)
    
  } else if (method == "ntraining"){
    # Study how decreasing the amount of training data affects outcomes
    
    res <- metricNTraining(cpData, cpMeta$Metadata_pert_id, epochs=epochs, reps=2)
    saveRDS(res, file=file.path(outpath, sprintf("%scp_ntraining_epch=%s.rds", epochs)))
    return(res)
  }
  
}


### Load Bray Data; Move to separate dataloaders file probs:
loadBrayData <- function(braypath, center=1, subsample=1){
  
  combds <- readRDS(braypath)
  
  ### Split data 
  ### Unnormalized data
  metax <- grep("Meta", colnames(combds))
  combMeta <- combds[, metax]
  combData <- as.matrix(combds[, setdiff(seq(dim(combds)[2]), metax)])
  
  if (center){
    combData <- scale(combData, center=TRUE, scale=TRUE)
    # apparently some axes do not vary
    combData[is.na(combData)] <- 0
  }
  
  # Assign controls unique pert_ids to make sure torch isn't grouping them, which could create massive memory hurdles  
  # This doesn't seem to be working for some reason, revisit later
  ix <- which(combMeta$Metadata_pert_type == "control")
  combMeta$Metadata_pert_id[ix] <- sprintf("ctl_%d", seq_along(ix))
  
  if (subsample < 1 & subsample > 0){
    perts <- sample(unique(combMeta$Metadata_pert_id), round(length(unique(combMeta$Metadata_pert_id))*subsample))
    jx <- which(combMeta$Metadata_pert_id %in% perts)
    
    combData <- combData[jx, ]
    combMeta <- combMeta[jx, ]
  }  
  
  return(list(ds=combData, metads=combMeta))
}


loadLincsData <- function(lincspath, samplecp=0, splitGrps=1){
  
  cpds <- readRDS(lincspath)
  
  metax <- grep("Meta", colnames(cpds))
  cpMeta <- cpds[, metax]
  cpData <- as.matrix(cpds[, -metax])
  
  cpMeta$Metadata_pert_iname[cpMeta$Metadata_pert_iname == ""] <- "empty"
  
  # There are a large number of posiive and negative controls. To help with training, I break them up
  # into groups of 60, which is the largest set of compounds apart from the main compounds
  
  if (splitGrps){
    cpcounts <- table(cpMeta$Metadata_pert_iname)
    bigCPs <- names(cpcounts[cpcounts > 60])
    
    for (mycp in bigCPs){
      ixgrps <- ceiling(sample(sum(cpMeta$Metadata_pert_iname == mycp))/60)
      cpMeta$Metadata_pert_id[cpMeta$Metadata_pert_iname == mycp] <- sprintf("%s_%d", cpMeta$Metadata_pert_id[cpMeta$Metadata_pert_iname == mycp], ixgrps)
      cpMeta$Metadata_pert_iname[cpMeta$Metadata_pert_iname == mycp] <- sprintf("%s_%d", mycp, ixgrps)
      
    }
  }
  
  #if (samplecp > 0){
  #  cpcounts <- table(cpMeta$Metadata_pert_iname)
  #}
  
  return(list(ds=cpData, metads=cpMeta))
}



analyzePermutedL1KData <- function(dspath, metapath, cell_id, outpath=".", method="xvalscram", epochs=10, saveModel=TRUE){
  l1k_meta <- read_l1k_meta(metapath, version=2020)
  
  if (cell_id == "all"){
    mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$pert_type == "trt_cp",]
    # This is a bit of a hack to get the learning approach to consider only same-cell line pairs
    mysigs$pert_iname <- sprintf("%s_%s", mysigs$pert_iname, mysigs$cell_id)
    
    # Mysigs has 196k unique pert+cell combinations. For expediency, sample:
    mysigs <- mysigs[mysigs$pert_iname %in% sample(unique(mysigs$pert_iname), 20000), ]
    
    ds <- parse_gctx(fname=get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=l1k_meta$siginfo$sig_id[l1k_meta$siginfo$pert_type == "trt_cp"])
    ds <- subset_gct(ds, cid=mysigs$sig_id)
  } else {
    mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$cell_id == cell_id & l1k_meta$siginfo$pert_type == "trt_cp",]
    ds <- parse_gctx(fname=get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  }
  
  print(sprintf("DS dimensions: %d x %d", dim(ds@mat)[1], dim(ds@mat)[2]))
  
  if (method == "xvalscram"){
    nFolds <- 3
    
    L1Kxval <- PermutedCrossValidate(t(ds@mat), mysigs$pert_iname, nFolds=nFolds, epochs=epochs)
    
    # Save the result
    saveRDS(L1Kxval$res_all, file=file.path(outpath, sprintf("L1Kxval_epch=%s_folds=%s_cell=%s.rds", epochs, nFolds, cell_id)))
    
    if(saveModel){
      for (ii in seq_along(L1Kxval$mymodels)){
        torch_save(L1Kxval$mymodel[[ii]], file.path(outpath, sprintf("L1Kxval_epch=%s_folds=%s_cell=%s_model%d.pt", epochs, nFolds, cell_id, ii)))
      }
    }
    
    return(L1Kxval)
  }
}

PermutedCrossValidate <- function(mat1, classes, nFolds=5, metric="", epochs=100, loss=mycos_t_loss){
  # Revise to intelligently choose folds in the case of a wide range of class sizes
  # But for now, do it naively
  
  mygroups <- table(classes)
  myfolds <- createFolds(names(mygroups), k=nFolds)
  myfolds <- sapply(myfolds, FUN=function(x) names(mygroups[x]))
  
  res_all <- list(myfolds=myfolds)
  mymodels <- list()
  
  for (testset in myfolds){
    testix <- which(classes %in% testset)
    trainix <- setdiff(seq_along(classes), testix)
    
    print("Permuting Training Labels")
    scrambledClasses <- sample(classes[trainix])
    
    gends_train <- genDataset(mat1[trainix,], scrambledClasses, pca_first = FALSE, scale = FALSE, center = FALSE)
    train_dl <- dataloader(gends_train, batch_size = 1, shuffle = TRUE)
    
    gends_valid <- genDataset(mat1[testix,], classes[testix], pca_first=FALSE, scale=FALSE, center=FALSE)
    valid_dl <- dataloader(gends_valid, batch_size = 1, shuffle = TRUE)
    
    model <- OneLayerLinear(dim(mat1)[2], dim(mat1)[2])
    device <- if (cuda_is_available()) torch_device("cuda:0") else "cpu"
    model <- model$to(device = device)
    optimizer <- optim_adam(model$parameters, lr = 0.01)
    
    res <- train_function(model=model, 
                          train_dl=train_dl, 
                          valid_dl=valid_dl, 
                          myloss=loss, 
                          device=device, 
                          optimizer=optimizer, 
                          epochs = epochs)
    
    print("Computing baseline sim")
    trainGrpSim <- innerProductGroups("cosine", mat1[trainix,], scrambledClasses, compact=1)
    validGrpSim <- innerProductGroups("cosine", mat1[testix,], classes[testix], compact=1)
    
    res$baseline_train_loss <- mean((mean(trainGrpSim$diff) - sapply(trainGrpSim$same, mean))/sd(trainGrpSim$diff))
    res$baseline_valid_loss <- mean((mean(validGrpSim$diff) - sapply(validGrpSim$same, mean))/sd(validGrpSim$diff))
    
    # Separate the final model from each fold to save the results
    res_all[sprintf("fold%d", length(res_all))] <- list(res[names(res) != "model"])
    mymodels <- c(mymodels, res$model)
  }
  
  return(list(res_all=res_all, mymodels=mymodels))
}


# Takes a vector and a list and maps the vector into the same structure as the list
# equivalent to relisting and unlisted list.
listify <- function(mylist, myvec){
  mylengths <- sapply(mylist, length)
  mynames <- names(mylist)
  
  newlist <- list()
  # Need to vectorize, but for the lists of interest, this is not costly
  ix <- 1
  for (ii in seq_along(mylist)){
    newlist[[mynames[ii]]] <- myvec[ix:(ix + mylengths[ii]-1)]
    ix <- ix + mylengths[ii]
  }
  names(newlist) <- mynames
  return(newlist)
}

