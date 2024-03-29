library(gridExtra)
library(ggplot2)
library(cmapR)
library(CMAPToolkit)
library(torch)


getL1KReps <- function(datapath, l1kmeta, cellid="A375", mymodel, datalabel="L1K Reps", 
                       mkfigs=FALSE, outfile="./l1kreps_myfig.pdf", dataMode=c("balancePDF", "allPDF")){
 
  dataMode <- match.arg(dataMode)
  
  if (cellid == "all"){
    l1kmeta$siginfo$pert_iname <- sprintf("%s_%s", l1kmeta$siginfo$pert_iname, l1kmeta$siginfo$cell_id)
    mysigs <- l1kmeta$siginfo[l1kmeta$siginfo$pert_type == "trt_cp",]
  } else {
    mysigs <- l1kmeta$siginfo[l1kmeta$siginfo$cell_id == cellid & l1kmeta$siginfo$pert_type == "trt_cp", ]
  }
  
  ds <- parse_gctx(CMAPToolkit::get_level5_ds(datapath), rid=l1kmeta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  
  # Necessary because JURKAT has 32 bad signatures where all values are 0
  if (sum(colSums(ds@mat^2) == 0) > 0){
    excl_cids <- ds@cid[which(colSums(ds@mat^2) == 0)]
    mysigs <- mysigs[!(mysigs$sig_id %in% excl_cids), ]
    ds <- parse_gctx(CMAPToolkit::get_level5_ds(datapath), rid=l1kmeta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  }
  
  # This probably isn't fair because these replicate classes were used in training. I think I need to cross validate, even for the illustrative example.
  inProdReps <- innerProductGroups(model=mymodel, t(ds@mat), mysigs$pert_iname, compact=1)
  cosReps <- innerProductGroups(model="cosine", t(ds@mat), mysigs$pert_iname, compact=1)
  
  if (mkfigs){
    if (dataMode == "balancePDF"){ 
      pdf(outfile, width=8, height=6)
      plot(density(inProdReps$diff, bw=0.01), col="blue", lwd=2, xlim=c(-1,1), 
           ylim=c(0, 1.1*max(max(density(inProdReps$diff, bw=0.01)$y), max(density(cosReps$diff, bw=0.01)$y))), 
           xlab="Similarity", ylab="Density", main=sprintf("%s - All Similarities", datalabel))
      lines(density(cosReps$diff, bw=0.01), col="forestgreen", lwd=2, lty=6)
      
      balSampML <- sapply(inProdReps$same, FUN=function(x) sample(x, min(length(x), 100)))
      balSampCos <- sapply(cosReps$same, FUN=function(x) sample(x, min(length(x), 100)))
      
      lines(density(unlist(balSampML), bw=0.01), col="red", lwd=2)
      lines(density(unlist(balSampCos), bw=0.01), col="purple", lwd=2, lty=6)
      legend(x="topleft", legend=c("Cosine replicates", "Cosine all pairs", "MetLearn replicates", "MetLearn all pairs"), 
             lwd=c(2,2,2,2), lty=c(6, 6, 1, 1), col=c("purple", "forestgreen", "red", "blue"))
      
      dev.off()
    } 
    else if (dataMode == "allPDF"){
      pdf(outfile, width=8, height=6)
      plot(density(inProdReps$diff, bw=0.01), col="blue", lwd=2, xlim=c(-1,1), 
           ylim=c(0, 1.1*max(max(density(inProdReps$diff, bw=0.01)$y), max(density(cosReps$diff, bw=0.01)$y))), 
           xlab="Similarity", ylab="Density", main=sprintf("%s - All Similarities", datalabel))
      lines(density(unlist(inProdReps$same), bw=0.01), col="red", lwd=2)
      lines(density(unlist(cosReps$same), bw=0.01), col="purple", lwd=2, lty=6)
      lines(density(cosReps$diff, bw=0.01), col="forestgreen", lwd=2, lty=6)
      legend(x="topleft", legend=c("Cosine replicates", "Cosine all pairs", "MetLearn replicates", "MetLearn all pairs"), 
             lwd=c(2,2,2,2), lty=c(6, 6, 1, 1), col=c("purple", "forestgreen", "red", "blue"))
      dev.off()
    }
    
    pdf(sprintf("%s_perthist.pdf", strsplit(outfile, ".pdf", fixed=TRUE)[[1]][1]), width=8, height=6)
    hist(log2(sapply(inProdReps$same, length)), xlab="Count of replicate pairs by pert", 
         main=sprintf("%s Pert Pair count", datalabel))
    dev.off()
    
  }
  
  #   MLranks <- listify(MLSims, rankVectors(unlist(MLSims), MLbg))
  #   cosranks <- listify(cosSims, rankVectors(unlist(cosSims), cosbg))
  
  mlRanks <- listify(inProdReps$same, rankVectors(unlist(inProdReps$same), inProdReps$diff))
  cosRanks <- listify(cosReps$same, rankVectors(unlist(cosReps$same), cosReps$diff))
  
  # Sample the balanced statistics N times
  mldf <- data.frame()
  cosdf <- data.frame()
  
  for (ii in seq(100)){
    print(ii)
    mlBal <- sapply(mlRanks, FUN=function(x) sample(x, min(length(x), 100)))
    cosBal <- sapply(cosRanks, FUN=function(x) sample(x, min(length(x), 100)))
    
    mldf <- rbind(mldf, data.frame(auRankBal=mean(unlist(mlBal)), 
                                   q01=mean(p.adjust(unlist(mlBal), method="fdr") < 0.01), 
                                   q05=mean(p.adjust(unlist(mlBal), method="fdr") < 0.05), 
                                   q25=mean(p.adjust(unlist(mlBal), method="fdr") < 0.25)))
    cosdf <- rbind(cosdf, data.frame(auRankBal=mean(unlist(cosBal)), 
                                     q01=mean(p.adjust(unlist(cosBal), method="fdr") < 0.01), 
                                     q05=mean(p.adjust(unlist(cosBal), method="fdr") < 0.05), 
                                     q25=mean(p.adjust(unlist(cosBal), method="fdr") < 0.25)))
  }
  
  
  # auRank - balanced
  # q-value less than 0.05
  # q-value less than 0.25
  # query stats? top 1%
  
  repstats <- data.frame(datalabel=datalabel, cellid=cellid, datatype="L1K Reps, All", 
                         model="MetLearn", 
                         auRankBal=1-mean(mldf$auRankBal), 
                         auRankBalse=sd(mldf$auRankBal), 
                         q01Bal=mean(mldf$q01), 
                         q01Balse=sd(mldf$q01),
                         q05Bal=mean(mldf$q05),
                         q05Balse=sd(mldf$q05), 
                         q25Bal=mean(mldf$q25), 
                         q25Balse=sd(mldf$q25))
  repstats <- rbind(repstats, data.frame(datalabel=datalabel, cellid=cellid, datatype="L1K Reps, All",
                        model="cosine",
                        auRankBal=1-mean(cosdf$auRankBal), 
                        auRankBalse=sd(cosdf$auRankBal), 
                        q01Bal=mean(cosdf$q01), 
                        q01Balse=sd(cosdf$q01),
                        q05Bal=mean(cosdf$q05),
                        q05Balse=sd(cosdf$q05), 
                        q25Bal=mean(cosdf$q25), 
                        q25Balse=sd(cosdf$q25)))
  return(repstats)
}


getL1KXValReps <- function(myfolds, models, datapath, l1kmeta, cellid="A375"){
  
  if (cellid == "all"){
    l1kmeta$siginfo$pert_iname <- sprintf("%s_%s", l1kmeta$siginfo$pert_iname, l1kmeta$siginfo$cell_id)
    mysigs <- l1kmeta$siginfo[l1kmeta$siginfo$pert_type == "trt_cp", ]
  } else {
    mysigs <- l1kmeta$siginfo[l1kmeta$siginfo$cell_id == cellid & l1kmeta$siginfo$pert_type == "trt_cp", ]
    
  }
  
  ds <- parse_gctx(CMAPToolkit::get_level5_ds(datapath), rid=l1kmeta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  
  # Necessary because JURKAT has 32 bad signatures where all values are 0
  if (sum(colSums(ds@mat^2) == 0) > 0){
    excl_cids <- ds@cid[which(colSums(ds@mat^2) == 0)]
    mysigs <- mysigs[!(mysigs$sig_id %in% excl_cids), ]
    ds <- parse_gctx(CMAPToolkit::get_level5_ds(datapath), rid=l1kmeta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  }
  
  inProdReps <- cosReps <- mlRanks <- cosRanks <- list()
  for (ii in seq_len(length(myfolds))) {
    mymodel <- torch::torch_load(models[ii])
    foldsigs <- mysigs[mysigs$pert_iname %in% myfolds[[ii]], ]
    ds_fold <- subset_gct(ds, cid=foldsigs$sig_id)
    
    inProdReps[[ii]] <- innerProductGroups(model=mymodel, t(ds_fold@mat), foldsigs$pert_iname, compact=1)
    cosReps[[ii]] <- innerProductGroups(model="cosine", t(ds_fold@mat), foldsigs$pert_iname, compact=1)
    
    mlRanks[[ii]] <- listify(inProdReps[[ii]]$same, rankVectors(unlist(inProdReps[[ii]]$same), inProdReps[[ii]]$diff))
    cosRanks[[ii]] <- listify(cosReps[[ii]]$same, rankVectors(unlist(cosReps[[ii]]$same), cosReps[[ii]]$diff))
  }
  
  # Don't actually do this; this is debugging. Reduce the data instead.
  # Actually, this is a fine object to return from which to generate all the plots. 
  # Note: this does not support query (column-wise) ranks, which can't use innerProductGroups
  # or innerProduct. Consider implementing with sampling for performance. 
  return(list(inProdReps=inProdReps, cosReps=cosReps, mlRanks=mlRanks, cosRanks=cosRanks))
}


getL1KMoA <- function(datapath, l1kmeta, cellid, mymodel, pclds, altMets=FALSE){

  if (cellid == "all"){
    l1kmeta$siginfo$pert_iname <- sprintf("%s_%s", l1kmeta$siginfo$pert_iname, l1kmeta$siginfo$cell_id)
    mysigs <- l1kmeta$siginfo[l1kmeta$siginfo$pert_type == "trt_cp",]
  } else {
    mysigs <- l1kmeta$siginfo[l1kmeta$siginfo$cell_id == cellid & l1kmeta$siginfo$pert_type == "trt_cp", ]
  }
  
  ds <- parse_gctx(CMAPToolkit::get_level5_ds(datapath), rid=l1kmeta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  
  # Necessary because JURKAT has 32 bad signatures where all values are 0
  if (sum(colSums(ds@mat^2) == 0) > 0){
    excl_cids <- ds@cid[which(colSums(ds@mat^2) == 0)]
    mysigs <- mysigs[!(mysigs$sig_id %in% excl_cids), ]
    ds <- parse_gctx(CMAPToolkit::get_level5_ds(datapath), rid=l1kmeta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  }
  
  if (altMets == FALSE){
    mymodel <- torch::torch_load(mymodel)
    
    mlPCLs <- innerProductPairwiseGroups(model=mymodel, mat1=t(ds@mat), classes=mysigs$pert_iname, 
                                             sets=pclds$pertnames, compact=1)
    cosPCLs <- innerProductPairwiseGroups(model="cosine", mat1=t(ds@mat), classes=mysigs$pert_iname, 
                                          sets=pclds$pertnames, compact=1)
    
    mlRanks <- listify(mlPCLs$setSims, rankVectors(unlist(mlPCLs$setSims), mlPCLs$allSims))
    cosRanks <- listify(cosPCLs$setSims, rankVectors(unlist(cosPCLs$setSims), cosPCLs$allSims))
    
    return(list(mlPCLs=mlPCLs, cosPCLs=cosPCLs, mlRanks=mlRanks, cosRanks=cosRanks))
  } else {
    print("Running alternate metrics")
    print("cosine")
    cosPCLs <- simPairwiseGroups(metric="cosine", mat1=t(ds@mat), classes=mysigs$pert_iname, sets=pclds$pertnames, compact=1) 
    print("wtcs")
    wtcsPCLs <- simPairwiseGroups(metric="fastwtcs", mat1=t(ds@mat), classes=mysigs$pert_iname, sets=pclds$pertnames, compact=1)
    print("pearson")
    pearsonPCLs <- simPairwiseGroups(metric="pearson", mat1=t(ds@mat), classes=mysigs$pert_iname, sets=pclds$pertnames, compact=1)
    print("spearman")
    spearPCLs <- simPairwiseGroups(metric="spearman", mat1=t(ds@mat), classes=mysigs$pert_iname, sets=pclds$pertnames, compact=1)
    
    cosRanks <- listify(cosPCLs$setSims, rankVectors(unlist(cosPCLs$setSims), cosPCLs$allSims))
    wtcsRanks <- listify(wtcsPCLs$setSims, rankVectors(unlist(wtcsPCLs$setSims), wtcsPCLs$allSims))
    pearsonRanks <- listify(pearsonPCLs$setSims, rankVectors(unlist(pearsonPCLs$setSims), pearsonPCLs$allSims))
    spearRanks <- listify(spearPCLs$setSims, rankVectors(unlist(spearPCLs$setSims), spearPCLs$allSims))
    
    return(list(cosPCLs=cosPCLs, wtcsPCLs=wtcsPCLs, pearsonPCLs=pearsonPCLs, spearPCLs=spearPCLs,
                cosRanks=cosRanks, wtcsRanks=wtcsRanks, pearsonRanks=pearsonRanks, spearRanks=spearRanks))
  }
}


getCPReps <- function(cpmat, mymodel, pertlabels, datalabel="datalabel"){
  
  inProdReps <- innerProductGroups(model=mymodel, cpmat, pertlabels, compact=1)
  cosReps <- innerProductGroups(model="cosine", cpmat, pertlabels, compact=1)
  
  mlRanks <- listify(inProdReps$same, rankVectors(unlist(inProdReps$same), inProdReps$diff))
  cosRanks <- listify(cosReps$same, rankVectors(unlist(cosReps$same), cosReps$diff))
  
  return(list(inProdReps=inProdReps, cosReps=cosReps, mlRanks=mlRanks, cosRanks=cosRanks, datalabel=datalabel))
}

getCPXvalReps <- function(cpmat, cpmeta, myfolds, models, datalabel="datalabel"){
  
  inProdReps <- cosReps <- mlRanks <- cosRanks <- list()
  
  for (ii in seq_len(length(myfolds))){
    print(sprintf("fold = %d", ii))
    mymodel <- torch::torch_load(models[ii])
    
    foldix <- which(cpmeta$Metadata_pert_id %in% myfolds[[ii]])
    
    # This is a hack
    if (length(foldix) == 0){
      print("using pert_iname")
      foldix <- which(cpmeta$Metadata_pert_iname %in% myfolds[[ii]])
    }
    
    foldmeta <- cpmeta[foldix, ]
    foldmat <- cpmat[foldix, ]
    
    inProdReps[[ii]] <- innerProductGroups(model=mymodel, foldmat, foldmeta$Metadata_pert_id, compact=1)
    cosReps[[ii]] <- innerProductGroups(model="cosine", foldmat, foldmeta$Metadata_pert_id, compact=1)
    
    mlRanks[[ii]] <- listify(inProdReps[[ii]]$same, rankVectors(unlist(inProdReps[[ii]]$same), inProdReps[[ii]]$diff))
    cosRanks[[ii]] <- listify(cosReps[[ii]]$same, rankVectors(unlist(cosReps[[ii]]$same), cosReps[[ii]]$diff))
  }
  
  return(list(inProdReps=inProdReps, cosReps=cosReps, mlRanks=mlRanks, cosRanks=cosRanks, datalabel=datalabel))
}


getCPMoA <- function(cpmat, cpmeta, mymodel, pertlabels, pclds, datalabel="datalabel"){
  if (length(unique(cpmeta$Metadata_cell_id)) > 1){
    retlist <- list()
    
    for (mycell in unique(cpmeta$Metadata_cell_id)){
      print(mycell)
      ix <- which(cpmeta$Metadata_cell_id == mycell)
      b <- getCPMoA(cpmat[ix, ], cpmeta[ix, ], mymodel, pertlabels[ix], pclds, datalabel=sprintf("%s %s", datalabel, mycell))
      retlist <- c(retlist, list(b))
    }
    
    names(retlist) <- unique(cpmeta$Metadata_cell_id)
    
    return(retlist)
    
  } else {
    mlPCLs <- innerProductPairwiseGroups(model=mymodel, mat1=cpmat, classes=pertlabels, 
                                         sets=pclds$pertids, compact=1)
    cosPCLs <- innerProductPairwiseGroups(model="cosine", mat1=cpmat, classes=pertlabels, 
                                          sets=pclds$pertids, compact=1)
    
    mlRanks <- listify(mlPCLs$setSims, rankVectors(unlist(mlPCLs$setSims), mlPCLs$allSims))
    cosRanks <- listify(cosPCLs$setSims, rankVectors(unlist(cosPCLs$setSims), cosPCLs$allSims))
    
    return(list(mlPCLs=mlPCLs, cosPCLs=cosPCLs, mlRanks=mlRanks, cosRanks=cosRanks, datalabel=datalabel))
  }
}



balancedSample <- function(mylist, k=100){
  
  return(sapply(mylist, FUN=function(x) sample(x, min(length(x), k))))
  
}


splitCPds <- function(cpds, ix){
  return(list(ds = cpds$ds[ix, ], metads=cpds$metads[ix,]))
}


summarizeCPMoA <- function(cpMoAds, dsname="", balSampT=100, iter=10){
  return(data.frame(dataset=c(dsname, dsname), 
             method=c("ML", "Cos"), 
             auROC=c(1-mean(sapply(seq(iter), FUN=function(x) mean(unlist(balancedSample(cpMoAds$mlRanks, k=balSampT))))),
                     1-mean(sapply(seq(iter), FUN=function(x) mean(unlist(balancedSample(cpMoAds$cosRanks, k=balSampT)))))), 
             fdr05=c(mean(sapply(seq(iter), FUN=function(x) mean(p.adjust(unlist(balancedSample(cpMoAds$mlRanks, k=balSampT)), method="fdr") < 0.05))),
                     mean(sapply(seq(iter), FUN=function(x) mean(p.adjust(unlist(balancedSample(cpMoAds$cosRanks, k=balSampT)), method="fdr") < 0.05)))), 
             fdr10=c(mean(sapply(seq(iter), FUN=function(x) mean(p.adjust(unlist(balancedSample(cpMoAds$mlRanks, k=balSampT)), method="fdr") < 0.1))),
                     mean(sapply(seq(iter), FUN=function(x) mean(p.adjust(unlist(balancedSample(cpMoAds$cosRanks, k=balSampT)), method="fdr") < 0.1)))), 
             fdr25=c(mean(sapply(seq(iter), FUN=function(x) mean(p.adjust(unlist(balancedSample(cpMoAds$mlRanks, k=balSampT)), method="fdr") < 0.25))),
                     mean(sapply(seq(iter), FUN=function(x) mean(p.adjust(unlist(balancedSample(cpMoAds$cosRanks, k=balSampT)), method="fdr") < 0.25))))))
}



summarizeMoAdf_ByMoA <- function(cpMOAds, pclLabels, dsname=""){
  
  # A hack to solve an omission on getCPMoA
  if (length(pclLabels) > length(cpMOAds$mlPCLs$setSims)){
    print(sprintf("Warning: truncating pclLabels to length %d", length(cpMOAds$mlPCLs$setSims)))
    pclLabels <- pclLabels[1: length(cpMOAds$mlPCLs$setSims)]
  }
  
  return(data.frame(dsname=dsname,
                    pcl=pclLabels,
                    mlAUC=sapply(cpMOAds$mlRanks, FUN=mean), 
                    cosAUC=sapply(cpMOAds$cosRanks, FUN=mean),
                    mlFDR05=sapply(cpMOAds$mlRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.05)),
                    cosFDR05=sapply(cpMOAds$cosRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.05)),
                    mlFDR10=sapply(cpMOAds$mlRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.1)),
                    cosFDR10=sapply(cpMOAds$cosRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.1)),
                    mlSNR=(sapply(cpMOAds$mlPCLs$setSims, mean) - mean(cpMOAds$mlPCLs$allSims))/sd(cpMOAds$mlPCLs$allSims),
                    cosSNR=(sapply(cpMOAds$cosPCLs$setSims, mean) - mean(cpMOAds$cosPCLs$allSims))/sd(cpMOAds$cosPCLs$allSims),
                    pclSize=sapply(cpMOAds$mlPCLs$setSims, length)))
  
}

