# Apply models to replicate and PCL prediction
benchmarkL1K_PCLs <- function(modelpath, dspath, metapath, pclpath, cell_id, outpath=".", modelname="mymodel", mkplots=1){
  pclds <- readRDS(pclpath)
  model <- torch_load(modelpath)
  pertnames <- unlist(pclds$pertnames)
  
  # Read appropriate data
  l1k_meta <- read_l1k_meta(metapath, version=2020)
  
  if (cell_id == "all"){
    mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$pert_type == "trt_cp",]

    ds <- parse_gctx(fname=get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=l1k_meta$siginfo$sig_id[l1k_meta$siginfo$pert_type == "trt_cp"])
    ds <- subset_gct(ds, cid=mysigs$sig_id)
  } else {
    mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$cell_id == cell_id & l1k_meta$siginfo$pert_type == "trt_cp",]
    ds <- parse_gctx(fname=get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  
  }
  
  mx <- sample(seq_len(dim(mysigs)[1]), min(dim(mysigs)[1], 2000))
  cosbg <- innerProduct("cosine", t(ds@mat[, mx]), t(ds@mat[,mx]))
  MLbg <- innerProduct(model, t(ds@mat[, mx]), t(ds@mat[,mx]))
  
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
      pclMLSim <- innerProduct(model, t(ds@mat[, ix]), t(ds@mat[,ix]))  
      pclCosSim  <- innerProduct("cosine", t(ds@mat[, ix]), t(ds@mat[,ix]))
      MLSims[[pclds$pcldf$pclid[ii]]] <- (pclMLSim[upper.tri(pclMLSim)])[filt]
      cosSims[[pclds$pcldf$pclid[ii]]] <- (pclCosSim[upper.tri(pclCosSim)])[filt]
    }
  }
  
  MLranks <- listify(MLSims, rankVectors(unlist(MLSims), MLbg))
  cosranks <- listify(cosSims, rankVectors(unlist(cosSims), cosbg))
  
  retdf <- list(modelname=modelname, MLSims=MLSims, cosSims=cosSims, MLRanks=MLranks, cosranks=cosranks, MLbg=MLbg, cosbg=cosbg)
  saveRDS(retdf, file.path(outpath, sprintf("%s_pcl_res.rds", modelname)))
  
  if (mkplots){
    pdf(file=file.path(outpath, sprintf("%s_pcl_figs.pdf", modelname)), width=8, height=6)
    # ML density
    plot(density(MLbg[upper.tri(MLbg)], bw=0.01), col="blue", xlim=c(-1,1), xlab="Metric Learning similarity", lwd=2, main=modelname)
    lines(density(unlist(MLSims), bw=0.01), col="red", lwd=2)
    legend(x="topright", legend=c("PCL pairs", "All pairs"), lwd=c(4,4), col=c("red", "blue"))
    
    # Cos Density
    plot(density(cosbg[upper.tri(cosbg)], bw=0.01), col="purple", xlim=c(-1,1), xlab="Cosine similarity", lwd=2, main=modelname)
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

