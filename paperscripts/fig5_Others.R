source("paperscripts/figinit.R")

library(DescTools)


#### Figure 5: Biological Interpretation

# Eigenvalue distribution 
eigencells <- c("A375", "A549", "HA1E", "MCF10A", "MCF7", "PC3", "VCAP", "ASC", "BT20", "HCC515", "HEK293", "HEPG2", "HUVEC", "NPC")
eigencells <- sort(eigencells)
fmodels <- list.files(file.path(l1kdir, "models"), pattern="L1Kmetric_")

# Load the eigendata - eigenvalues of the covariance matrices of the spaces
if (!file.exists(file.path(outdir, "../figspaper_res", "eigenvaluedata.rds"))){
  eigendata <- list()
  
  for (mycell in eigencells){
    ds <- parse_gctx(get_level5_ds(datapath), cid=siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_type == "trt_cp"], rid = landmarks$pr_gene_id)
    
    print(sprintf("%s: %s", mycell, fmodels[grep(mycell, fmodels)]))
    
    mymodel <- torch_load(file.path(l1kdir, "models", fmodels[grep(mycell, fmodels)]))
    modelmat <- as.matrix(mymodel(torch_tensor(t(ds@mat), dtype=torch_float())))
    
    # Centering is inappropriate because we are taking inner products in the uncentered space. 
    system.time(pcBase <- prcomp(t(ds@mat), center=FALSE))
    pcML <- prcomp(modelmat, center=FALSE)
    
    pctvarBase <- pcBase$sdev^2/sum(pcBase$sdev^2)
    pctvarML <- pcML$sdev^2/sum(pcML$sdev^2)
    
    x <- t(ds@mat) %*% pcBase$rotation
    mlRotation <- as.matrix(mymodel(torch_tensor(t(pcBase$rotation), dtype=torch_float())))
    xML <- modelmat %*% t(mlRotation)
    
    # CHECK: Is this correct? Not sure we want to compute the sd of the projection of the datasets onto the PCs
    origPCVars <- apply(x, 2, sd)^2
    mlPCVars <- apply(xML, 2, sd)^2
    
    origPCVarPct <- origPCVars/sum(origPCVars)
    mlPCVarPct <- mlPCVars/sum(mlPCVars)
    
    eigendata[[mycell]] <- list(pcBase=pcBase[c("sdev", "rotation")], pcML=pcML[c("sdev", "rotation")], 
                                pctvarBase=pctvarBase, pctvarML=pctvarML,
                                origPCVars=origPCVars, mlPCVars=mlPCVars)
  }
  
  saveRDS(eigendata, file=file.path(outdir, "../figspaper_res", "eigenvaluedata.rds"))
  
} else {
  eigendata <- readRDS(file.path(outdir, "../figspaper_res", "eigenvaluedata.rds"))
  attach(eigendata)
}

# Load the eigendata for brayds
if (!file.exists(file.path(outdir, "../figspaper_res", "eigendataBray.rds"))){
  braydir <- file.path(cpdir, "bray/models")

  brayds <- loadBrayData(braypath)
  
  # Bray PCA is going crazy because of a small number (< 1%) of crazy outliers. 
  # It's unclear how much these would affect metric learning, as the loss function
  # is insensitive under rescaling. For PCA, I am removing them from both base and 
  # embedded computations because they contribute massively to the variance
  
  xsum <- rowSums(abs(brayds$ds))
  jx <- which((xsum - mean(xsum))/sd(xsum) > 10)
  
  braymodel <- torch::torch_load(file.path(braydir, "braymetric_epch=10_smp=50_model.pt"))
  braymodelmat <- as.matrix(braymodel(torch_tensor(brayds$ds, dtype=torch_float())))
  
  ix <- sample(setdiff(seq(dim(brayds$ds)[1]), jx), 10000)
  
  # Normalize - project onto the unit ball 
  baseMat <- brayds$ds[ix,]/sqrt(rowSums(brayds$ds[ix,]^2))
  mlMat <- braymodelmat[ix,]/sqrt(rowSums(braymodelmat[ix,]^2))
  
  #system.time(pcBrayBase <- prcomp(brayds$ds[ix,], center=FALSE))  # The base space is already scaled
  #system.time(pcBrayML <- prcomp(braymodelmat[ix,], center=FALSE))
  
  system.time(pcBrayBase <- prcomp(baseMat, center=FALSE))
  system.time(pcBrayML <- prcomp(mlMat, center=FALSE))
  
  pctvarBaseBray <- pcBrayBase$sdev^2/sum(pcBrayBase$sdev^2)
  pctvarMLBray <- pcBrayML$sdev^2/sum(pcBrayML$sdev^2)
  
  eigenBray <- list(pcBrayBase=pcBrayBase[c("sdev", "rotation")], 
                    pcBrayML=pcBrayML[c("sdev", "rotation")], 
                    pctvarBaseBray=pctvarBaseBray, 
                    pctvarMLBray=pctvarMLBray, 
                    sample_ix=ix)
  saveRDS(eigenBray, file=file.path(outdir, "../figspaper_res", "eigendataBray.rds"))
  
} else {
  eigenBray <- readRDS(file.path(outdir, "../figspaper_res", "eigendataBray.rds"))
}

eigenBraydf <- rbind(data.frame(dataset="base", pctVar=pctvarBaseBray, cumVar=cumsum(pctvarBaseBray), index=seq_along(pctvarBaseBray)), 
                     data.frame(dataset="ML", pctVar=pctvarMLBray, cumVar=cumsum(pctvarMLBray), index=seq_along(pctvarMLBray)))

#### Get hallmark variance - a supervised approach to a biological interpretation of changes in the gene space 
# hmarks loaded in figinit

if (!file.exists(file.path(outdir, "../figspaper_res", "hallmarkVariance.rds"))){
  hmarkVar <- list()
  hmarkVarRand <- list()
  
  # Consider pct variance along unit vector pointing in direction of each hallmark geneset.
  # This is necessarily extraordinarily limited; better probably would be to consider the *subspace*
  # corresponding to the genes in the hallmark geneset
  for (mycell in eigencells){
    hmarkCell <- data.frame()
    hmarkCellRand <- data.frame()
    ds <- parse_gctx(get_level5_ds(datapath), cid=siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_type == "trt_cp"], rid = landmarks$pr_gene_id)
    
    print(sprintf("%s: %s", mycell, fmodels[grep(mycell, fmodels)]))
    
    mymodel <- torch_load(file.path(l1kdir, "models", fmodels[grep(mycell, fmodels)]))
    modelmat <- as.matrix(mymodel(torch_tensor(t(ds@mat), dtype=torch_float())))
    
    pcML <- eigendata[[mycell]]$pcML
    
    lengthBase <- matLength(t(ds@mat))
    lengthML <- matLength(modelmat)
    
    for (ii in seq_along(hmarks)){
      print(ii)
      # unit vector
      hvec <- as.numeric(ds@rid %in% hmarks[[ii]]$entry)/sqrt(sum(ds@rid %in% hmarks[[ii]]$entry))
      
      rvec <- as.numeric(ds@rid %in% sample(ds@rid, sum(ds@rid %in% hmarks[[ii]]$entry)))/sqrt(sum(sum(ds@rid %in% hmarks[[ii]]$entry)))
      
      hvecRot <- as.matrix(mymodel(torch_tensor(hvec, dtype=torch_float())))
      hvecScale <- sqrt(sum(hvecRot^2))
      hvecRot <- hvecRot/sqrt(sum(hvecRot^2))
      
      rvecRot <- as.matrix(mymodel(torch_tensor(rvec, dtype=torch_float())))
      rvecScale <- sqrt(sum(rvecRot^2))
      rvecRot <- rvecRot/sqrt(sum(rvecRot^2))
      
      hprodBase <- t(ds@mat) %*% hvec
      hprodML <- modelmat %*% hvecRot
      
      rprodBase <- t(ds@mat) %*% rvec
      rprodML <- modelmat %*% rvecRot
      
      hmarkCell <- rbind(hmarkCell, data.frame(hmarkSet=hmarks[[ii]]$head, 
                                               cellid=mycell, 
                                               baseVar=sum(hprodBase^2)/sum(lengthBase^2), 
                                               mlVar=sum(hprodML^2)/sum(lengthML^2),
                                               hvecScale=hvecScale, 
                                               maxeigen=pcML$sdev[1]))
      
      hmarkCellRand <- rbind(hmarkCellRand, data.frame(hmarkSet=sprintf("%s:%s", hmarks[[ii]]$head, "random"), 
                                                       cellid=mycell,
                                                       baseVar=sum(rprodBase^2)/sum(lengthBase^2), 
                                                       mlVar=sum(rprodML^2)/sum(lengthML^2),
                                                       rvecScale=rvecScale, 
                                                       maxeigen=pcML$sdev[1]))
  
    }
    hmarkVar[[mycell]] <- hmarkCell
    hmarkVarRand[[mycell]] <- hmarkCellRand
  }
  
  # Revised embedding along unit vector pointing in direction of each hallmark geneset
  hVarEmbed <- list()
  
  for (mycell in eigencells){
    hmarkCell <- data.frame()
    
    ds <- parse_gctx(get_level5_ds(datapath), cid=siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_type == "trt_cp"], rid = landmarks$pr_gene_id)
    
    print(sprintf("%s: %s", mycell, fmodels[grep(mycell, fmodels)]))
    
    mymodel <- torch_load(file.path(l1kdir, "models", fmodels[grep(mycell, fmodels)]))
    modelmat <- as.matrix(mymodel(torch_tensor(t(ds@mat), dtype=torch_float())))
    
    svdmat <- svd(as.matrix(mymodel$fc1$weight))
    
    embedmat <- svdmat$v %*% diag(svdmat$d) %*% t(svdmat$v) %*% ds@mat
    
    lengthBase <- matLength(t(ds@mat))
    lengthML <- matLength(t(embedmat))
    
    for (ii in seq_along(hmarks)){
      print(ii)
      # unit vector
      hvec <- as.numeric(ds@rid %in% hmarks[[ii]]$entry)/sqrt(sum(ds@rid %in% hmarks[[ii]]$entry))
      hsize <- sum(ds@rid %in% hmarks[[ii]]$entry)
      
      hprodBase <- t(ds@mat) %*% hvec
      hprodML <- t(embedmat) %*% hvec
      
      if (hsize > 1){
        hvarBase <- colSums(ds@mat[which(hvec > 0), ]^2)
        hvarML <- colSums(embedmat[which(hvec > 0), ]^2)
      }
      
      hmarkCell <- rbind(hmarkCell,
                         data.frame(hmarkSet=hmarks[[ii]]$head, 
                                    size=hsize, 
                                    cellid=mycell, 
                                    baseVar=sum(hprodBase^2)/sum(lengthBase^2), 
                                    mlVar=sum(hprodML^2)/sum(lengthML^2), 
                                    subspVarbase=sum(hvarBase)/sum(lengthBase^2),
                                    subspVarML=sum(hvarML)/sum(lengthML^2)))
      
    }
    hVarEmbed[[mycell]] <- hmarkCell
  }
  
  saveRDS(file=file.path(outdir, "../figspaper_res", "hallmarkVariance.rds"), 
          list(hmarkVar=hmarkVar, hmarkVarRand=hmarkVarRand, hVarEmbed=hVarEmbed))
} else {
  hmarkVards <- readRDS(file.path(outdir, "../figspaper_res", "hallmarkVariance.rds"))
  attach(hmarkVards)
}


#### Unsupervised interpretation of biological changes in the embedded space

if (!file.exists(file.path(outdir, "../figspaper_res", "eigenOrigEigenVars.rds"))){
  origEigenVars <- list()
    
  for (mycell in eigencells){
    ds <- parse_gctx(get_level5_ds(datapath), cid=siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_type == "trt_cp"], rid = landmarks$pr_gene_id)
    
    print(sprintf("%s: %s", mycell, fmodels[grep(mycell, fmodels)]))
    
    mymodel <- torch_load(file.path(l1kdir, "models", fmodels[grep(mycell, fmodels)]))
    modelmat <- as.matrix(mymodel(torch_tensor(t(ds@mat), dtype=torch_float())))
    
    svdmat <- svd(as.matrix(mymodel$fc1$weight))
    
    embedmat <- svdmat$v %*% diag(svdmat$d) %*% t(svdmat$v) %*% ds@mat
    
    lengthBase <- matLength(t(ds@mat))
    lengthML <- matLength(t(embedmat))
    
    # Get variance in each direction of the base PCs:
    baseMatRot <- t(ds@mat) %*% eigendata[[mycell]]$pcBase$rotation
    embedmatRot <- t(embedmat) %*% eigendata[[mycell]]$pcBase$rotation
    baseVar <- colMeans(baseMatRot^2)
    embedVar <- colMeans(embedmatRot^2)
    
    # Have to normalize the variance, as scalar multiplication affects the computed variance:
    baseVar <- baseVar/sum(baseVar)
    embedVar <- embedVar/sum(embedVar)
    
    origEigenVars[[mycell]] <- list(baseVar=baseVar, embedVar=embedVar)
    # Sanity check: baseVar should correlate very close to 1 with the calculated pcvar:
    print(sprintf("Correlation of base variances (should be nearly 1): %0.4f", cor(baseVar, eigendata[[mycell]]$pcBase$sdev^2)))
  }
  
  saveRDS(origEigenVars, file=file.path(outdir, "../figspaper_res", "eigenOrigEigenVars.rds"))
} else {
  origEigenVars <- readRDS(file.path(outdir, "../figspaper_res", "eigenOrigEigenVars.rds"))
}

#### Gini coefficients of eigenvalues

if (!file.exists(file.path(outdir, "../figspaper_res", "eigendataGini.rds"))){
  baseVar <- sapply(eigendata, FUN=function(x) x$pctvarBase)
  MLVar <- sapply(eigendata, FUN=function(x) x$pctvarML)
  
  giniCoefs <- rbind(data.frame(ds="L1000 ML", gini=sapply(seq(14), FUN=function(x) Gini(MLVar[1:200,x])), cellid=names(eigendata)),
                     data.frame(ds="L1000 base", gini=sapply(seq(14), FUN=function(x) Gini(baseVar[1:200,x])), cellid=names(eigendata)),
                     data.frame(ds="CDRP ML", gini=Gini(eigenBray$pctvarMLBray[1:100]), cellid="CDRP"), 
                     data.frame(ds="CDRP base", gini=Gini(eigenBray$pctvarBaseBray[1:100]), cellid="CDRP"))
  
  
  saveRDS(giniCoefs, file.path(outdir, "../figspaper_res", "eigendataGini.rds"))
} else {
  giniCoefs <- readRDS(file.path(outdir, "../figspaper_res", "eigendataGini.rds"))
}


# Ghetto plots
pdf(file.path(outdir, "fig5_eigenvalueDistributionFigures.pdf"), width=8, height=6)
plot(-100, -100, xlim=c(0, 978), ylim=c(-5, 0), xlab="Eigenvalue index", ylab="Percent Variance explained", main="Eigenvalues of L1000 data distribution")
for (ii in seq_along(eigendata)){
  lines(log10(eigendata[[ii]]$pctvarBase), col="blue", type="l", lwd=2)
  lines(log10(eigendata[[ii]]$pctvarML), col="red", type="l", lwd=2)
}
legend(x="bottomleft", legend=c("Base (cosine)", "Embedding (ML)"), col=c("blue", "red"), lwd=c(3,3))

plot(-100, -100, xlim=c(0,978), ylim=c(-3, 2), xlab="Eigenvalue index", ylab="log10 Rescaling", main="Rescaling of Eigenvalues of L1000 data distribution")
for (ii in seq_along(eigendata)){
  lines(log10(eigendata[[ii]]$pctvarML/eigendata[[ii]]$pctvarBase), col="black", type="l", lwd=2)
}
grid()

plot(-10, -10, xlim=c(0, 978), ylim=c(0,1), xlab="Eigenvalue index", ylab="Cumulative pct variance explained", main="Cumulative Percentage of Variance of L1000 data distribution")
for (ii in seq_along(eigendata)){
  lines(cumsum(eigendata[[ii]]$pctvarBase), col="blue", type="l", lwd=2)
  lines(cumsum(eigendata[[ii]]$pctvarML), col="red", type="l", lwd=2)
}
legend(x="bottomright", legend=c("Base (cosine)", "Embedding (ML)"), col=c("blue", "red"), lwd=c(3,3))
dev.off()


# GGplot version of the cumulative variance distribution. 
eigendf <- data.frame(cellid=character(), dataset=character(), x=numeric(), cumpctvar=numeric())
for (ii in seq_along(eigendata)){
  eigendf <- rbind(eigendf, data.frame(cellid=names(eigendata)[ii], x=seq(978), dataset="base", cumpctvar=cumsum(eigendata[[ii]]$pctvarBase)))
  eigendf <- rbind(eigendf, data.frame(cellid=names(eigendata)[ii], x=seq(978), dataset="ML", cumpctvar=cumsum(eigendata[[ii]]$pctvarML)))
}

pdf(file.path(outdir, "fig5_eigenvalueDistributionFiguresGGP.pdf"), width=8, height=6)
ggplot(eigendf, aes(x=x, y=cumpctvar, color=dataset, lty=cellid)) + geom_line(lwd=1) + theme_minimal() + 
  coord_cartesian(xlim=c(0, 250)) + xlab("Eigenvalue Index") + ylab("Cumulative percent variance") + ggtitle("L1000 Pct Variance")
dev.off()


pcdf <- data.frame(cellid=names(eigendata), 
                   pc25=c(sapply(eigendata, FUN=function(x) min(which(cumsum(x$pcBase$sdev^2/sum(x$pcBase$sdev^2)) > 0.25))), 
                          sapply(eigendata, FUN=function(x) min(which(cumsum(x$pcML$sdev^2/sum(x$pcML$sdev^2)) > 0.25)))), 
                   pc50=c(sapply(eigendata, FUN=function(x) min(which(cumsum(x$pcBase$sdev^2/sum(x$pcBase$sdev^2)) > 0.50))), 
                          sapply(eigendata, FUN=function(x) min(which(cumsum(x$pcML$sdev^2/sum(x$pcML$sdev^2)) > 0.50)))), 
                   pc95=c(sapply(eigendata, FUN=function(x) min(which(cumsum(x$pcBase$sdev^2/sum(x$pcBase$sdev^2)) > 0.95))), 
                          sapply(eigendata, FUN=function(x) min(which(cumsum(x$pcML$sdev^2/sum(x$pcML$sdev^2)) > 0.95)))), 
                   pc5e3=c(sapply(eigendata, FUN=function(x) sum((x$pcBase$sdev^2/sum(x$pcBase$sdev^2) > 0.005))), 
                           sapply(eigendata, FUN=function(x) sum((x$pcML$sdev^2/sum(x$pcML$sdev^2) > 0.005)))),
                   pc1e2=c(sapply(eigendata, FUN=function(x) sum((x$pcBase$sdev^2/sum(x$pcBase$sdev^2) > 0.01))), 
                           sapply(eigendata, FUN=function(x) sum((x$pcML$sdev^2/sum(x$pcML$sdev^2) > 0.01)))),
                   space=c(rep("Base", 14), rep("ML", 14)))


pdf(file.path(outdir, "fig5_eigenvaluePlots2.pdf"), width=8, height=6)
ggplot(pcdf, aes(x=cellid, y=pc95, fill=space)) + geom_bar(stat="identity", position="dodge") + 
  theme_minimal() + ggtitle("Number of PCs needed to account for 95% of variance")
ggplot(pcdf, aes(x=cellid, y=pc5e3, fill=space)) + geom_bar(stat="identity", position="dodge") + 
  theme_minimal() + ggtitle("Number of PCs with at least 0.5% of variance")

plot(-10, -10, xlim=c(0, 250), ylim=c(0,1), xlab="Eigenvalue index", ylab="Cumulative pct variance explained", main="Cumulative Percentage of Variance of L1000 data distribution")
  for (ii in seq_along(eigendata)){
    lines(cumsum(eigendata[[ii]]$pctvarBase), col="coral2", type="l", lwd=2)
    lines(cumsum(eigendata[[ii]]$pctvarML), col="cyan4", type="l", lwd=2)
  }
  legend(x="bottomright", legend=c("Base (cosine)", "Embedding (ML)"), col=c("coral2", "cyan4"), lwd=c(3,3))

ggplot(giniCoefs, aes(x=factor(ds, levels=c("L1000 base", "L1000 ML", "CDRP base", "CDRP ML")), y=gini)) + 
         geom_violin(aes(fill=giniCoefs$ds)) + geom_jitter(width = 0.25) + theme_minimal() + 
  xlab("Embedding") + ylab("Gini Coefficient") + ggtitle("Gini coefficients of leading eigenvalues, Wilcox p = 3.4e-3") + 
  ylim(c(0.4, 1)) + theme(legend.position="none")


ggplot(rbind(data.frame(space="ML", eigenvalue=as.numeric(MLVar)), 
             data.frame(space="base", eigenvalue=as.numeric(baseVar))), 
       aes(x=eigenvalue, color=space, fill=space)) + geom_density(alpha=0.6) + 
      coord_cartesian(xlim=c(1e-10, 1)) + scale_x_continuous(trans="log10") + theme_minimal() + 
  xlab("Eigenvalue (pct variance)") + ylab("Density") + ggtitle("L1000 Eigenvalue distributions")
dev.off()


pdf(file.path(outdir, "fig5_CDRPEigenvaluePlots.pdf"), width=8, height=6)
ggplot(eigenBraydf, aes(x=index, y=cumVar, color=dataset)) + geom_line(lwd=2) + coord_cartesian(xlim=c(0,100)) + theme_minimal() + 
  xlab("Eigenvalue Index") + ylab("Cumulative pct variance explained") + ggtitle("CDRP Cell Painting Pct Variance explained")

ggplot(eigenBraydf, aes(x=pctVar, color=dataset, fill=dataset)) + geom_density(alpha=0.4) + 
  scale_x_continuous(trans="log10") + coord_cartesian(xlim=c(1e-12, 1)) + theme_minimal() + 
  xlab("Eigenvalue (pct variance)") + ylab("Density") + ggtitle("CDRP Eigenvalue distributions")
dev.off()

pdf(file.path(outdir, "Sfig_ginicoefs.pdf"), width=8, height=6)
ggplot(giniCoefs, aes(x=factor(ds, levels=c("L1000 base", "L1000 ML", "CDRP base", "CDRP ML")), y=gini)) + 
  geom_violin(aes(fill=giniCoefs$ds)) + geom_jitter(width = 0.25, size=3) + theme_minimal() + 
  xlab("Embedding") + ylab("Gini Coefficient") + ggtitle("Gini coefficients of leading eigenvalues, Wilcox p = 1.7e-3") + 
  ylim(c(0.4, 1)) + theme(legend.position="none")
dev.off()


# Summarize the variance of the two spaces along the base eigenvectors
pdf(file.path(outdir, "../biointerp/baseEigenvectorVariance.pdf"))
for (ii in seq_along(origEigenVars)){
  myds <- origEigenVars[[ii]]
  plot(log10(myds$baseVar), log10(myds$embedVar), pch=16, xlim=c(-5,0), ylim=c(-5,0), xlab="log10 Base pct variance", ylab="log10 Embedded pct variance", 
       main=sprintf("%s percent variance along base eigenvectors", names(origEigenVars)[ii]), col=rainbow(14)[ii])
  lines(c(-10,10), c(-10,10), col="grey", lty=2)
}
dev.off()


biominq <- sapply(bioRes, FUN=function(y) sapply(y, FUN=function(x) -log10(min(x$padj))))

biomindf <- rbind(data.frame(cellid=eigencells, count = colSums(biominq > 2), condition="q < 1e-2"), 
                  data.frame(cellid=eigencells, count = colSums(biominq > 1), condition="q < 1e-1"))

pdf(file.path(outdir, "../biointerp/baseEigenvectorQCount.pdf"), width=8, height=6)
ggplot(biomindf, aes(x=cellid, y=count, fill=condition)) + geom_bar(position="dodge", stat="identity") + 
  theme_minimal() + xlab("Cell Line") + ylab("N PCs with at least 1 significant gene set") + ggtitle("300 base PC q-value count among GO c5 genesets")
dev.off()