source("figinit.R")




#### Figure 5: Biological Interpretation

# Eigenvalue distribution 
eigencells <- c("A375", "A549", "HA1E", "MCF10A", "MCF7", "PC3", "VCAP", "ASC", "BT20", "HCC515", "HEK293", "HEPG2", "HUVEC", "NPC")
eigencells <- sort(eigencells)
fmodels <- list.files(file.path(l1kdir, "models"), pattern="L1Kmetric_")

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
}


#### Get hallmark variance
# hmarks loaded in figinit
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
  
  pcML <- eigendata[[mycell]]
  
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


pcdf <- data.frame(cellid=names(eigendata), 
                   pc95=c(sapply(eigendata, FUN=function(x) min(which(cumsum(x$pcBase$sdev^2/sum(x$pcBase$sdev^2)) > 0.95))), 
                          sapply(eigendata, FUN=function(x) min(which(cumsum(x$pcML$sdev^2/sum(x$pcML$sdev^2)) > 0.95)))), 
                   pc5e3=c(sapply(eigendata, FUN=function(x) sum((x$pcBase$sdev^2/sum(x$pcBase$sdev^2) > 0.005))), 
                           sapply(eigendata, FUN=function(x) sum((x$pcML$sdev^2/sum(x$pcML$sdev^2) > 0.005)))),
                   pc1e2=c(sapply(eigendata, FUN=function(x) sum((x$pcBase$sdev^2/sum(x$pcBase$sdev^2) > 0.01))), 
                           sapply(eigendata, FUN=function(x) sum((x$pcML$sdev^2/sum(x$pcML$sdev^2) > 0.01)))),
                   space=c(rep("Base", 14), rep("ML", 14)))

ggplot(pcdf, aes(x=cellid, y=pc95, fill=space)) + geom_bar(stat="identity", position="dodge") + theme_minimal() + ggtitle("Number of PCs needed to account for 95% of variance")
ggplot(pcdf, aes(x=cellid, y=pc5e3, fill=space)) + geom_bar(stat="identity", position="dodge") + theme_minimal() + ggtitle("Number of PCs with at least 0.5% of variance")


# Centering is inappropriate because we are taking inner products in the uncentered space. 

pdf(file.path(outdir, "fig6_VarianceExplained_HEPG2.pdf"), width=8, height=6)
ggplot(data.frame(index=seq(978), native=pctvar, embedded=pctvarML), aes(x=native, y=embedded)) + geom_point() + scale_x_log10() + scale_y_log10() + 
  geom_abline(intercept=0, slope=1, col="blue", lty=2) + theme_minimal() + coord_cartesian(xlim=c(1e-4, 3e-1), ylim=c(1e-4, 3e-1)) + xlab("Native Pct Variance") +
  ylab("Embedded Pct Variance") + ggtitle("Variance Explained by Metric Learning in HEPG2")
dev.off()

pdf(file.path(outdir, "fig6_EigenvalueRatio_HEPG2.pdf"), width=8, height=6)
ggplot(data.frame(index=seq(978), native=pctvar, embedded=pctvarML, ratio=pctvarML/pctvar), aes(x=index, y=ratio)) + geom_point() + theme_minimal() + 
  xlab("Eigenvalue index") + ylab("Ratio of eigenvalues of embedding to native") + ggtitle("Ratio of eigenvalues of HEPG2 Embedding vs native")
dev.off()


# I'm not sure why these aren't identical:
x <- t(ds@mat) %*% pcHEPG2$rotation

mlRotation <- as.matrix(mymodel(torch_tensor(t(pcHEPG2$rotation), dtype=torch_float())))
xML <- modelmat %*% t(mlRotation)

plot(apply(x, 2, sd), pcHEPG2$sdev)

origPCVars <- apply(x, 2, sd)^2
mlPCVars <- apply(xML, 2, sd)^2

origPCVarPct <- origPCVars/sum(origPCVars)
mlPCVarPct <- mlPCVars/sum(mlPCVars)


scaledf <- data.frame(logScaleFactor = log10(mlPCVarPct/origPCVarPct), pcIX = seq(978), cellid=mycell)
pcVars <- data.frame(cellid=mycell, pcix=seq(978), origPCVarPct = origPCVarPct, mlPCVarPct = mlPCVarPct)

pdf(file.path(outdir, "fig6_VarRescaled_HEPG2.pdf"), width=8, height=6)
ggplot(scaledf, aes(x=pcIX, y=logScaleFactor, color=cellid)) + geom_point() + theme_minimal() + xlab("PC Index") + ylab("Log10 rescaling factor") + 
  ggtitle("Rescaling of original principal components in embedded space")
dev.off()


#### Figure 7: Cosine Theory (MOVE TO METASIG PAPER) ####
# Figure 6b: see bhk/code/cosine/cosineScript.R

#### Fig 6a:
# pdf(file.path(outdir, "fig2a_cosineVsDim_dists.pdf"), width=7, height=5)
# plot(density(x1000[upper.tri(x1000)], bw=0.01), col="forestgreen", lwd=3, xlim=c(-1,1), 
#      main="Distribution of cosine similarities for standard normals of dimension N")
# lines(density(x10[upper.tri(x10)], bw=0.01), col="red", lwd=3)
# lines(density(x100[upper.tri(x100)], bw=0.01), col="blue", lwd=3)
# legend(x="topleft", legend=c("N=20", "N=100", "N=500"), lwd=c(4,4,4), col=c("red", "darkblue", "forestgreen"))
# dev.off()

#### Fig 6a Revised:
datapath <- "~/Work/bhk/data/l1k/2020/"
l1kmeta <- CMAPToolkit::read_l1k_meta(datapath, version=2020)
attach(l1kmeta)

mycells <- c("A375", "A549", "ASC", "HCC515", "HEK293", "HEPG2", "MCF7", "NPC", "PC3")
l1kCellDf <- data.frame(cellid = character(), mean = numeric(), sd = numeric())


pdf(file.path(outdir, "fig6_L1K_landmark_cos.pdf"), width=7, height=6)
x10 <- getCosineDist(mydim=20, N=1000, mkfig=0)
x100 <- getCosineDist(mydim=100, N=1000, mkfig=0)
cos978 <- getCosineDist(mydim=978, N=1000, mkfig=0)
plot(density(cos978[upper.tri(cos978)], bw=0.01), col="black", lwd=2, xlim=c(-0.5, 0.5), xlab="Cosine Similarity", ylab="Density", main="L1000 Cosine Similarities, landmark")
lines(density(x10[upper.tri(x10)], bw=0.01), col="black", lty=2, lwd=2)
lines(density(x100[upper.tri(x100)], bw=0.01), col="black", lty=4, lwd=2)

for (cella in mycells){
  print(cella)
  ds <- parse_gctx(get_level5_ds(datapath), rid = geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], 
                   cid=siginfo$sig_id[siginfo$cell_id == cella & siginfo$pert_type == "trt_cp"])
  ix <- sample(length(ds@cid), 2000)
  dscos <- CMAPToolkit::cosine(ds@mat[, ix], ds@mat[, ix])
  
  l1kCellDf <- rbind(l1kCellDf, data.frame(cellid = cella, mean = mean(dscos[upper.tri(dscos)]), sd = sd(dscos[upper.tri(dscos)])))
  lines(density(dscos[upper.tri(dscos)], bw=0.01), col=rainbow(10)[match(cella, mycells)], lwd=2)
}
# All cell lines
ds <- parse_gctx(get_level5_ds(datapath), rid = geneinfo$pr_gene_id[geneinfo$pr_is_lm == 1], 
                 cid=siginfo$sig_id[siginfo$pert_type == "trt_cp" & siginfo$cell_id != "JURKAT"])   # Jurkat has a handful of signatures that are all 0s for some reason.
ix <- sample(length(ds@cid), 4000)
dscos <- CMAPToolkit::cosine(ds@mat[, ix], ds@mat[, ix])
l1kCellDf <- rbind(l1kCellDf, data.frame(cellid = "all", mean = mean(dscos[upper.tri(dscos)]), sd = sd(dscos[upper.tri(dscos)])))
lines(density(dscos[upper.tri(dscos)], bw=0.01), col=rainbow(10)[10], lwd=2)

legend(x="topright", legend=c("Null: 978 genes", "Null: 100", "Null: 10", mycells, "all"), col=c("black", "black", "black", rainbow(10)), lwd=3, 
       lty=c(1,2,4,rep(1,10)))

dev.off()

l1kCellDf$effectDim <- (1/l1kCellDf$sd)^2
saveRDS(l1kCellDf, file=file.path("~/Work/bhk/analysis/cosine/l1kCellDf_cosineDist.RDS"))


printdf <- rbind(data.frame(cellid="Std Norm 978", mean=mean(cos978[upper.tri(cos978)]), sd=sd(cos978[upper.tri(cos978)]), effectDim=978), l1kCellDf[1:8,])
pdf(file.path(outdir, "L1K_landmark_table.pdf"), width=4, height=6)
grid.table(printdf %>% mutate(across(where(is.numeric), ~ round(., 5))), rows=NULL)
dev.off()

# Supp Fig: 
pdf(file.path(outdir, "fig6_cosineFit_supp.pdf"), width=7, height=5)
x100 <- getCosineDist(mydim=50, N=1000, mkfig=1)
dev.off()


# L1000 SD vs ngenes
l1kdf <- readRDS("~/Work/bhk/analysis/cosine/L1K_HEPG2_geneSubsample_MeanVsDim.rds")
l1kdf$effectDim <- (1/l1kdf$sd)^2

# Add jitter to the 978 to avoid tripping up geom_violin
# l1kdf$sd[l1kdf$ngenes %in% c(978, 12328)] <- l1kdf$sd[l1kdf$ngenes %in% c(978, 12328)] + rnorm(sum(l1kdf$ngenes %in% c(978, 12328)), 0, 0.003)
# l1kdf$effectDim[l1kdf$ngenes %in% c(978, 12328)] <- l1kdf$effectDim[l1kdf$ngenes %in% c(978, 12328)] + rnorm(sum(l1kdf$ngenes %in% c(978, 12328)), 0, 1)
# pdf(file.path(outdir, "fig6de_cosineVsNGenesL1K.pdf"), width=7, height=6, onefile = TRUE)
# topplot <- ggplot(l1kdf, aes(x=ngenes, y=sd, group=factor(ngenes))) + geom_violin(fill="lightblue") + theme_minimal() + 
#   scale_x_log10(breaks=c(30, 100, 300, 1000, 3000, 10000), labels=c("30", "100", "300", "1000", "3000", "10000")) + 
#   xlab("Number of Genes") + ylab("Stdev of cosine similarity") #+ geom_smooth(method="loess") 
# botplot <- ggplot(l1kdf, aes(x=ngenes, y=effectDim, group=factor(ngenes))) + geom_violin(fill="lightgreen") + theme_minimal() + 
#   scale_x_log10(breaks=c(30, 100, 300, 1000, 3000, 10000), labels=c("30", "100", "300", "1000", "3000", "10000")) +
#   xlab("Number of Genes") + ylab("Dimensional Equivalent") #+ geom_smooth(method="loess") 
# grid.arrange(topplot, botplot)
# dev.off()

pdf(file.path(outdir, "fig6_cosineVsNGenesL1K.pdf"), width=7, height=6, onefile = TRUE)
topplot <- ggplot(l1kdf, aes(x=ngenes, y=sd)) + geom_point(position="jitter", alpha=0.8, colour="darkblue") + theme_minimal() + 
  scale_x_log10(breaks=c(30, 100, 300, 1000, 3000, 10000), labels=c("30", "100", "300", "1000", "3000", "10000")) + 
  xlab("Number of Genes") + ylab("Stdev of cosine similarity") + geom_smooth(method="loess") 
botplot <- ggplot(l1kdf, aes(x=ngenes, y=effectDim)) + geom_point(position="jitter", alpha=0.8, colour="darkgreen") + theme_minimal() +
  scale_x_log10(breaks=c(30, 100, 300, 1000, 3000, 10000), labels=c("30", "100", "300", "1000", "3000", "10000")) +
  xlab("Number of Genes") + ylab("Dimensional Equivalent") + geom_smooth(method="loess")
grid.arrange(topplot, botplot)
dev.off()

