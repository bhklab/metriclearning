library(CMAPToolkit)
library(cmapR)
library(torch)

# Analyze L1000 datasets
analyzeL1KData <- function(dspath, metapath, cell_id, outpath=".", method="xval", epochs=10, saveModel=FALSE){
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




# if (cell_id == "all"){
#   mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$pert_type == "trt_cp",]
#   # This is a bit of a hack to get the learning approach to consider only same-cell line pairs
#   mysigs$pert_iname <- sprintf("%s_%s", mysigs$pert_iname, mysigs$cell_id)
#   
#   # Mysigs has 196k unique pert+cell combinations. For expediency, sample:
#   mysigs <- mysigs[mysigs$pert_iname %in% sample(unique(mysigs$pert_iname), 20000), ]
#   
#   ds <- parse_gctx(fname=get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=l1k_meta$siginfo$sig_id[l1k_meta$siginfo$pert_type == "trt_cp"])
#   ds <- subset_gct(ds, cid=mysigs$sig_id)
# } else {
#   mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$cell_id == cell_id & l1k_meta$siginfo$pert_type == "trt_cp",]
#   ds <- parse_gctx(fname=get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)
# }



# Analyze Bray dataset
analyzeBrayData <- function(braypath, outpath=".", method="xval", epochs=10, saveModel=FALSE){
  epochs <- as.numeric(epochs)
  
  # brayds <- readRDS(braypath)
  # 
  # metax <- grep("Meta", colnames(brayds))
  # brayMeta <- brayds[, metax]
  # brayData <- as.matrix(brayds[, setdiff(seq(dim(brayds)[2]), metax)])
  # 
  # brayData <- scale(brayData, center=TRUE, scale=TRUE)
  # # apparently some axes do not vary
  # # Check for infinities
  # brayData[is.na(brayData)] <- 0
  
  brayds <- loadBradyData(braypath, center=1)
  brayData <- brayds$ds
  brayMeta <- brayds$metads
  
  if (method == "xval"){
    nFolds <- 3
    brayxval <- metricCrossValidate(brayData, brayMeta$Metadata_pert_id, nFolds=nFolds, epochs=epochs) 
    saveRDS(brayxval$res_all, file=file.path(outpath, sprintf("brayxval_epch=%d_folds=%d.rds", epochs, nFolds)))
    
    if (saveModel){
      for (ii in seq_along(brayxval$mymodels)){
        torch_save(brayxval$mymodel[[ii]], file.path(outpath, sprintf("brayxval_epch=%d_folds=%d_model%d.pt", epochs, nFolds, ii)))
      }
    }
    return(brayxval)
    
  } else if (method == "allds"){
    braymetric <- learnInnerProduct(brayData, brayMeta$Metadata_pert_id, epochs=epochs)
    saveRDS(braymetric, file=file.path(outpath, sprintf("braymetric_epch=%d.rds", epochs)))
    torch_save(braymetric$model, file.path(outpath, sprintf("braymetric_epch=%d_model.pt", epochs)))
    return(braymetric)
    
  } else if (method == "ntraining"){
    # Study how decreasing the amount of training data affects outcomes
    
    res <- metricNTraining(brayData, brayMeta$Metadata_pert_id, epochs=epochs, reps=2)
    saveRDS(res, file=file.path(outpath, sprintf("brayntraining_epch=%s.rds", epochs)))
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


### Load Move to separate dataloaders file probs:
loadBrayData <- function(braypath, center=1){
  
  combds <- readRDS(braypath)
  
  ### Split data 
  ### Unnormalized data
  metax <- grep("Meta", colnames(combds))
  combMeta <- combds[, metax]
  combData <- as.matrix(combds[, setdiff(seq(dim(combds)[2]), metax)])
  
  if (center){
    combDataNorm <- scale(combData, center=TRUE, scale=TRUE)
    # apparently some axes do not vary
    combDataNorm[is.na(combDataNorm)] <- 0
  }
  
  # Assign controls unique pert_ids to make sure torch isn't grouping them, which could create massive memory hurdles  
  # This doesn't seem to be working for some reason, revisit later
  ix <- which(combMeta$Metadata_pert_type == "control")
  combMeta$Metadata_pert_id[ix] <- sprintf("ctl_%d", seq_along(ix))
  
  #return(list(ds=combData[-ix,], metads=combMeta[-ix,]))
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
      cpMeta$Metadata_pert_iname[cpMeta$Metadata_pert_iname == mycp] <- sprintf("%s_%d", mycp, ixgrps)
    }
  }
  
  #if (samplecp > 0){
  #  cpcounts <- table(cpMeta$Metadata_pert_iname)
  #}
  
  return(list(ds=cpData, metads=cpMeta))
}
