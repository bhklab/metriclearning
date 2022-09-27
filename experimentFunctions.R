# Analyze L1000 datasets





# Analyze Bray dataset
analyzeBrayData <- function(braypath, outpath=".", method="xval"){
  brayds <- readRDS(braypath)
  
  metax <- grep("Meta", colnames(brayds))
  brayMeta <- brayds[, metax]
  brayData <- as.matrix(brayds[, setdiff(seq(dim(brayds)[2]), metax)])
  
  brayData <- scale(brayData, center=TRUE, scale=TRUE)
  # apparently some axes do not vary
  brayData[is.na(brayData)] <- 0
  
  if (method == "xval"){
    epochs <- 10
    nFolds <- 5
    brayxval <- metricCrossValidate(brayData, brayMeta$Metadata_pert_id, nFolds=nFolds, epochs=epochs) 
    saveRDS(brayxval, file=file.path(outpath, sprintf("brayxval_epch=%d_folds=%d.rds", epochs, nFolds)))
    return(brayxval)
  } else if (method == "allds"){
    epochs <- 250
    braymetric <- learnInnerProduct(brayData, brayMeta$Metadata_pert_id, epochs=epochs)
    saveRDS(braymetric, file=file.path(outpath, sprintf("braymetric_epch=%d.rds", epochs)))
    return(braymetric)
  } else if (method == "ntraining"){
    # Study how decreasing the amount of training data affects outcome
  }
  
}




# Make ROC plots, summary plots