library(CMAPToolkit)
library(cmapR)
library(torch)

# Analyze L1000 datasets
analyzeL1KData <- function(dspath, metapath, cell_id, outpath=".", method="xval", epochs=10){
  l1k_meta <- read_l1k_meta(metapath, version=2020)
  
  if (cell_id == "all"){
    mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$pert_type == "trt_cp",]
    # This is a bit of a hack to get the learning approach to consider only same-cell line pairs
    mysigs$pert_iname <- sprintf("%s_%s", mysigs$pert_iname, mysigs$cell_id)
  } else {
    mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$cell_id == cell_id & l1k_meta$siginfo$pert_type == "trt_cp",]
  }
  
  ds <- parse_gctx(fname=get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  
  if (method == "xval"){
    nFolds <- 5
    
    L1Kxval <- metricCrossValidate(t(ds@mat), mysigs$pert_iname, nFolds=nFolds, epochs=epochs)
    saveRDS(L1Kxval, file=file.path(outpath, sprintf("L1Kxval_epch=%s_folds=%s_cell=%s.rds", epochs, nFolds, cell_id)))
    torch_save(L1Kxval$model, file.path(outpath, sprintf("L1Kxval_epch=%s_folds=%s_cell=%s_model.pt", epochs, nFolds, cell_id)))
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



# Analyze Bray dataset
analyzeBrayData <- function(braypath, outpath=".", method="xval", epochs=10){
  epochs <- as.numeric(epochs)
  brayds <- readRDS(braypath)
  
  metax <- grep("Meta", colnames(brayds))
  brayMeta <- brayds[, metax]
  brayData <- as.matrix(brayds[, setdiff(seq(dim(brayds)[2]), metax)])
  
  brayData <- scale(brayData, center=TRUE, scale=TRUE)
  # apparently some axes do not vary
  brayData[is.na(brayData)] <- 0
  
  if (method == "xval"){
    nFolds <- 3
    brayxval <- metricCrossValidate(brayData, brayMeta$Metadata_pert_id, nFolds=nFolds, epochs=epochs) 
    saveRDS(brayxval, file=file.path(outpath, sprintf("brayxval_epch=%d_folds=%d.rds", epochs, nFolds)))
    torch_save(brayxval$model, file.path(outpath, sprintf("brayxval_epch=%d_folds=%d_model.pt", epochs, nFolds)))
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
