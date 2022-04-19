library(cmapR)
library(CMAPToolkit)
library(ggplot2)
library(parallel)

test_innerprod <- function(datapath, metapath, cell_id="HEPG2", pcdim=978, epochlim=1000, outdir="."){
 
  l1k_meta <- read_l1k_meta(metapath, version=2020)
  attach(l1k_meta)
  
  mysigs <- siginfo[siginfo$cell_id == cell_id & siginfo$pert_type == "trt_cp", ]
  
  ds <- parse_gctx(get_level5_ds(datapath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)
  
  pcspace <- prcomp(t(ds@mat), scale=FALSE, center=FALSE)
  
  cpcount <- table(mysigs$pert_iname)
  mycmpds <- names(cpcount)[which(cpcount >= 12)]
  
  wts <- numeric(pcdim)+1
  m <- numeric(pcdim)
  v <- numeric(pcdim)
  epoch <- 1
  
  wmat <- matrix(numeric(epochlim*pcdim), ncol=pcdim)
  
    while (epoch <= epochlim){
      if (epoch %% 25 == 0){
        print(sprintf("epoch = %d", epoch))
      }
      mycp <- sample(mycmpds, 1)
      
      ix <- which(mysigs$pert_iname == mycp)
      jx <- sample(setdiff(seq_len(dim(pcspace$x)[1]), ix), 2*length(ix))
      
      mymat <- cbind(t(pcspace$x[ix, 1:pcdim]), t(pcspace$x[jx, 1:pcdim]))
      myclass <- c(rep(1, length(ix)), rep(0, length(jx)))
      
      retvals <- update_step_par(mymat, myclass, wts, m, v, epoch)
      wmat[epoch,] <- retvals$w
      
      wts <- retvals$w
      m <- retvals$m
      v <- retvals$v
      epoch <- epoch + 1
    }
  
  wmat[wmat < 0] <- 0
  saveRDS(list(pcspace=pcspace, wmat=wmat, mysigs=mysigs), file=file.path(outdir, sprintf("%s_metric_N=%d.rds", cell_id, epochlim)))
  return(list(pcspace=pcspace, wmat=wmat, mysigs=mysigs))
}



train_metric <- function(mymat, classes, epochlim=1000){
  dim1 <- dim(mymat)[1]
  cpcount <- table(classes)
  mycmpds <- names(cpcount)[which(cpcount >= 12)]
  
  wts <- numeric(dim1)+1
  m <- numeric(dim1)
  v <- numeric(dim1)
  epoch <- 1
  
  wmat <- matrix(numeric(epochlim*dim1), ncol=dim1)

  while (epoch <= epochlim){
    if (epoch %% 25 == 0){
      print(sprintf("epoch = %d", epoch))
    }
    mycp <- sample(mycmpds, 1)
    
    ix <- which(classes == mycp)
    jx <- sample(setdiff(seq_len(dim(mymat)[2]), ix), 2*length(ix))
    
    amat <- cbind(mymat[,ix], mymat[,jx])
    myclass <- c(rep(1, length(ix)), rep(0, length(jx)))
    
    retvals <- update_step(amat, myclass, wts, m, v, epoch)
    wmat[epoch,] <- retvals$w
    
    wts <- retvals$w
    m <- retvals$m
    v <- retvals$v
    epoch <- epoch + 1
  }
  
  wmat[wmat < 0] <- 0
  return(list(wmat=wmat))
  
}


eval_func <- function(mymat, classes, rotmat, weights, mkfigs=1, outdir=".", name="myfunccompare"){
  mysimmat0 <- cosine(mymat, mymat)
  
  #mysimmat <- cosine(sqrt(weights)*mymat, sqrt(weights)*mymat)
  mysimmat <- calc_simmat(ds1, ds2, rotmat, wts)

  
  print("Similarity matrices computed")
  simres <- data.frame(rep1= numeric(1e6), rep2=numeric(1e6), compound=character(1e6), sim0=numeric(1e6), sim1=numeric(1e6), rank0=numeric(1e6), rank1=numeric(1e6))
  
  repclass <- table(classes)[table(classes) > 1]
  mypos <- 1
  
  for (myclass in names(repclass)){
    ix <- which(classes == myclass)
    tres <- expand.grid(rep1=ix, rep2=ix)
    tres <- tres[tres$rep1 != tres$rep2,]
    tres <- tres[order(tres$rep1),]
    
    tres$compound <- myclass
    tres[,c("sim0", "sim1", "rank0", "rank1")] <- numeric(dim(tres)[1])
    
    for (ii in seq(dim(tres)[1])){
      ix1 <- tres$rep1[ii]
      ix2 <- tres$rep2[ii]
      tres[ii, seq(4, dim(tres)[2])] <- c(mysimmat0[ix1, ix2], mysimmat[ix1, ix2], 
                                          mean(mysimmat0[,ix2] > mysimmat0[ix1,ix2]), mean(mysimmat[, ix2] > mysimmat[ix1,ix2]))
    }
    simres[mypos:(mypos + dim(tres)[1] - 1),] <- tres
    mypos <- mypos + dim(tres)[1]
  }
  
  simres <- simres[1:max(which(simres$rep1 > 0)),]
  
  mean0 <- mean(sample(mysimmat0, 1e6))
  sd0 <- sd(sample(mysimmat0, 1e6))
  mean1 <- mean(sample(mysimmat, 1e6))
  sd1 <- sd(sample(mysimmat, 1e6))
  
  simres$z0 <- (simres$sim0 - mean0)/sd0
  simres$z1 <- (simres$sim1 - mean1)/sd1
  
  if (mkfigs){
    pdf(file.path(outdir, sprintf("%s_base_sim_density.pdf", name)), width=10, height=8)
    plot(density(sample(mysimmat0, 1e6), bw=0.01), col="blue", lwd=2, xlab="Cosine Similarity", ylab="Density, bw=0.01", main=sprintf("Base Similarity for %s", name), xlim=c(-1,1))
    lines(density(simres$sim0, bw=0.01), col="red", lwd=2)
    legend(x="topleft", legend=c("All signatures", "Duplicates"), col=c("blue", "red"), lwd=c(3,3))
    dev.off()
    
    pdf(file.path(outdir, sprintf("%s_mod_sim_density.pdf", name)), width=10, height=8)
    plot(density(sample(mysimmat, 1e6), bw=0.01), col="blue", lwd=2, xlab="Cosine Similarity", ylab="Density, bw=0.01", main=sprintf("Modified Similarity for %s", name), xlim=c(-1,1))
    lines(density(simres$sim1, bw=0.01), col="red", lwd=2)
    legend(x="topleft", legend=c("All signatures", "Duplicates"), col=c("blue", "red"), lwd=c(3,3))
    dev.off()
    
    pdf(file.path(outdir, sprintf("%s_comparison_hexplot.pdf", name)), width=10, height=8)
    print(ggplot(simres, aes(x=z0, y=z1)) + stat_binhex(bins=100) + geom_abline(intercept=0, slope=1, color="red", linetype="dashed", size=1.5) + xlab("Base replicate similarity") + 
      ylab("Modified replicate similarity") + ggtitle(sprintf("Comparison of z-scored replicate similarities for %s", name)) + theme_minimal() + xlim(c(-1,1)))
    dev.off()
    
    pdf(file.path(outdir, sprintf("%s_snr_comparison.pdf", name)), width=10, height=8)
    plot(ecdf(simres$z0), col="blue", lwd=2, xlab="SNR Z-scores", ylab="Cumulative Fraction", main=sprintf("CDF of replicate similarity Z-scores, %s", name))
    lines(ecdf(simres$z1), col="red", lwd=2)
    legend(x="topleft", legend=c("All signatures", "Duplicates"), col=c("blue", "red"), lwd=c(3,3))
    dev.off()
    
    pdf(file.path(outdir, sprintf("%s_rank_comparison.pdf", name)), width=10, height=8)
    plot(ecdf(simres$rank0), col="blue", lwd=2, xlab="Replicate rank", ylab="Cumulative fraction", main=sprintf("CDF of replicate rank, %s", name))
    lines(ecdf(simres$rank1), col="red", lwd=2)
    legend(x="topleft", legend=c("All signatures", "Duplicates"), col=c("blue", "red"), lwd=c(3,3))
    dev.off()
    
    saveRDS(simres, file=file.path(outdir, sprintf("%s_simres.rds", name)))
  }
  
  return(simres)
}