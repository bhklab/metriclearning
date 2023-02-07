library(gridExtra)
library(ggplot2)
library(cmapR)
library(CMAPToolkit)
#library(dplyr)
library(torch)
library(reshape2)
library(hexbin)
library(RColorBrewer)

### This is a script to generate figures for the manuscript.


topdir <- "~/Work/bhk/analysis/metric_learning/2022_L1K/"
cpdir <- "~/Work/bhk/analysis/metric_learning/2022_cellpaint"
outdir <- file.path(topdir, "figspaper")

datapath <- "~/Work/bhk/data/l1k/2020/"
l1kmeta <- CMAPToolkit::read_l1k_meta(datapath, version=2020)
attach(l1kmeta)

braypath <- "~/Work/bhk/data/cellpainting/bray-2017/bray_2017_combined_profiles.rds"
lincspath <- "~/Work/bhk/data/cellpainting/lincs-cell-painting/lincs-cell-painting/spherized_profiles/profiles"
lincs1 <- file.path(lincspath, "2016_04_01_a549_48hr_batch1_dmso_spherized_profiles_with_input_normalized_by_dmso.rds")
lincs2 <- file.path(lincspath, "2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_dmso.rds")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7") #, "#F0E442"


#### Fig 1 workflow - matrix ####
heatthemes <- theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), 
                    legend.position="none", panel.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), plot.background=element_blank())


myds <- data.frame(expand.grid(X=seq(150), Y=seq(120)), val=rbeta(150*120,1,1))
png(file=file.path(outdir, "fig1_matrix.png"), width=400, height=320)
ggplot(myds, aes(x=X, y=Y, fill=val)) + geom_tile() + scale_fill_distiller(palette = "RdBu") + theme_minimal() + xlab("") + ylab("") + 
  heatthemes
dev.off()

png(file=file.path(outdir, "fig1_matrixLg.png"), width=400, height=320)
ggplot(myds, aes(x=X, y=Y, fill=val)) + geom_tile() + scale_fill_distiller(palette = "RdBu") + theme_minimal() + xlab("") + ylab("") + 
  heatthemes
dev.off()

png(file=file.path(outdir, "fig1_matrixmed.png"), width=400, height=180)
ggplot(data.frame(expand.grid(X=seq(120), Y=seq(60)), val=rbeta(120*60,1,1)), aes(x=X, y=Y, fill=val)) + geom_tile() + scale_fill_distiller(palette = "RdBu") + theme_minimal() + xlab("") + ylab("") + 
  heatthemes
dev.off()

png(file=file.path(outdir, "fig1_matrixsm.png"), width=400, height=110)
ggplot(data.frame(expand.grid(X=seq(120), Y=seq(30)), val=rbeta(30*120,1,1)), aes(x=X, y=Y, fill=val)) + geom_tile() + scale_fill_distiller(palette = "RdBu") + theme_minimal() + xlab("") + ylab("") + 
  heatthemes
dev.off()



# Fig 1 Cosine loss
png(file=file.path(outdir, "fig1_cos.png"), width=400, height=320)
#pdf(file=file.path(outdir, "fig1_cos.pdf"), width=6, height=6)
plot(seq(-1,1, 0.02), 0.5*dbeta(seq(0,1,0.01),25,25), col="black", lwd=3, type="l", xlab="Similarity", ylab="", main="", yaxt="n")
lines(seq(-1,1,0.01), 0.7*dnorm(seq(-1,1,0.01), mean = 0.4, sd=0.15) + 0.3*dnorm(seq(-1,1,0.01), mean=0, sd=0.143), col="magenta", lwd=3)
dev.off()


# Fig 1 - signature
sigA <- rbeta(20, 1,1)
sigB <- rbeta(20, 1,1)
sigA2 <- 0.8*sigA + 0.2*rbeta(20,1,1)
png(file=file.path(outdir, "fig1_sigA.png"), width=200, height=30)
ggplot(data.frame(X=seq(20), Y=rep(1,20), val=sigA), aes(x=X, y=Y, fill=val)) + geom_tile() + scale_fill_distiller(palette = "RdBu") + heatthemes
dev.off()
png(file=file.path(outdir, "fig1_sigA2.png"), width=200, height=30)
ggplot(data.frame(X=seq(20), Y=rep(1,20), val=sigA2), aes(x=X, y=Y, fill=val)) + geom_tile() + scale_fill_distiller(palette = "RdBu") + heatthemes
dev.off()
png(file=file.path(outdir, "fig1_sigB.png"), width=200, height=30)
ggplot(data.frame(X=seq(20), Y=rep(1,20), val=sigB), aes(x=X, y=Y, fill=val)) + geom_tile() + scale_fill_distiller(palette = "RdBu") + heatthemes
dev.off()

# Fig 1 - CDF
a <- c(rnorm(10000, 0.4, 0.25), rnorm(0, 0.143))
b <- 2*rbeta(10000, 5, 5) - 1
png(file=file.path(outdir, "fig1_rankplot.png"), width=400, height=320)
plot(seq(0,1,1e-4), c(0,sort(1-rankVectors(a,b))), type="l", lwd=4, col="blue", ylim=c(0,1), xlab="Rank", ylab="Cumulative Fraction")
dev.off()


#### Figure 2: L1000 performance ####

l1kdir <- file.path(topdir, "modelRuns")

### Fig 2a - cosine similarity for replicates vs all
# Pick a representative cell line, say A375

##mymod <- torch::torch_load(file.path(l1kdir, "models", "L1Kmetric_epch=20_cell=A375_model.pt"))
#xvalds <- readRDS(file.path(l1kdir, "L1Kxval_epch=10_folds=5_cell=A375.rds"))
#f <- list.files(file.path(l1kdir, "models"), pattern="L1Kxval.*A375")

mycell <- "A375"
xvalds <- readRDS(file.path(l1kdir, list.files(l1kdir, pattern=sprintf("L1Kxval.*%s.rds", mycell))))
f <- list.files(file.path(l1kdir, "models"), pattern=sprintf("L1Kxval.*%s", mycell))

xvalreps <- getL1KXValReps(myfolds=xvalds$myfolds, 
                           models=file.path(l1kdir, "models", f), 
                           datapath=datapath,
                           l1kmeta=l1kmeta, 
                           cellid=mycell)

mlsame <- sapply(xvalreps$inProdReps, FUN=function(x) balancedSample(x$same, k=100))
mldiff <- sapply(xvalreps$inProdReps, FUN=function(x) x$diff)

cossame <- sapply(xvalreps$cosReps, FUN=function(x) balancedSample(x$same, k=100))
cosdiff <- sapply(xvalreps$cosReps, FUN=function(x) x$diff)

pdf(file.path(outdir, sprintf("fig2a_%s_pdf.pdf", mycell)), width=8, height=6)
plot(density(unlist(cosdiff), bw=0.01), col="forestgreen", lwd=1.5, lty=5, xlim=c(-0.5,1), 
     xlab="Similarity", ylab="Density", main=sprintf("%s Balanced Similarity", mycell))
lines(density(unlist(cossame), bw=0.01), col="orange", lwd=1.5, lty=5)
lines(density(unlist(mldiff), bw=0.01), col="blue", lwd=1.5)
lines(density(unlist(mlsame), bw=0.01), col="firebrick3", lwd=1.5)
lines(c(0, 0), c(0, 100), col="grey", lty=2)
legend(x="topright", legend=c("Metric learning replicates", "Metric Learning all pairs", 
    "Cosine replicates", "Cosine all pairs"), lwd=2.5, col=c("firebrick3", "blue", "purple", "forestgreen"), 
    lty=c(1,1,5,5))
dev.off()

### Fig 2b - one example of cross validated replicate recall curves

pdf(file.path(outdir, sprintf("fig2b_%s_xvalBalRankCDF.pdf", mycell)), width=8, height=6)
plot(c(-10), c(-10), xlim=c(0,1), ylim=c(0,1), xlab="Balanced Rank", ylab="Cumulative Fraction", 
         main=sprintf("%s balanced rank", mycell))
for (ii in seq(5)){
  lines(ecdf(unlist(balancedSample(xvalreps$cosRanks[[ii]], 100))), col="orange", lty=5, lwd=1.5)
}
for (ii in seq(5)){
  lines(ecdf(unlist(balancedSample(xvalreps$mlRanks[[ii]], 100))), col="firebrick3", lwd=1.5)
}
legend(x="bottomright", legend=c("Metric Learning", "Cosine"), col=c("firebrick3", "orange"), lwd=2.5,
       lty=c(1,5), bg="white")
dev.off()

# Fig 2c - summary cross validated replicate recall (bar plots, auRank?)

# Get cell names
mycells <- sapply(strsplit(sapply(strsplit(list.files(l1kdir, pattern="L1Kxval"), "cell="), 
                                  FUN=function(x) x[[2]]), ".rds"), FUN=function(x) x[[1]])

# !!!!!! Fix this once you have the all_xval models
mycells <- setdiff(mycells, c("all", "MCF10A", "SKL"))
xvalres <- list()

for (ii in seq_along(mycells)){
  mycell <- mycells[ii]
  print(sprintf("%d: %s", ii, mycell))
  
  xvalds <- readRDS(file.path(l1kdir, list.files(l1kdir, pattern=sprintf("L1Kxval.*%s.rds", mycell))))
  f <- list.files(file.path(l1kdir, "models"), pattern=sprintf("L1Kxval.*%s", mycell))
  
  xvalres[[mycell]] <- getL1KXValReps(myfolds=xvalds$myfolds, 
                             models=file.path(l1kdir, "models", f), 
                             datapath=datapath,
                             l1kmeta=l1kmeta, 
                             cellid=mycell)
  
}

saveRDS(xvalres, file=file.path(outdir, "../figspaper_res/L1K_xvalres.rds"))

### mean balanced AUC
# Consider also plotting all 100 balanced AUCs rather than taking the mean for each fold
# The results are stored as ranks on [0,1], where 0 is the best rank. 
# To convert to AUC, take 1 - mean balanced rank.

if (!file.exists(file.path(outdir, "../figspaper_res/L1K_MeanBalAUC.rds"))){
  L1KPertBalAUCML <- 1 - sapply(mycells, FUN=function(x) 
    sapply(xvalres[[x]]$mlRanks, FUN=function(y) 
      mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(y)))))))
  
  L1KPertBalAUCCos <- 1 - sapply(mycells, FUN=function(x) 
    sapply(xvalres[[x]]$cosRanks, FUN=function(y) 
      mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(y)))))))
  
  saveRDS(list(L1KPertBalAUCML=L1KPertBalAUCML, L1KPertBalAUCCos=L1KPertBalAUCCos), 
          file=file.path(outdir, "../figspaper_res/L1K_MeanBalAUC.rds"))
} else {
  balAucList <- readRDS(file.path(outdir, "../figspaper_res/L1K_MeanBalAUC.rds"))
  attach(balAucList)
}

#ggplot(df, aes(x=Var2, y=value, fill=method)) + geom_boxplot() + facet_wrap(~Var2, scale="free")
pdf(file=file.path(outdir, "fig2c_L1KMeanBalAUC.pdf"), width=8, height=6)
ggplot(rbind(data.frame(melt(L1KPertBalAUCCos), method="Cosine"), data.frame(melt(L1KPertBalAUCML), method="Metric Learning")), aes(x=Var2, y=value, fill=method)) + 
  geom_boxplot() + geom_point(position=position_jitterdodge()) + 
  theme_minimal() + ylim(c(0.5, 1)) + geom_hline(yintercept=0.5, col="grey", linetype=2) + 
  xlab("Cell Line") + ylab("Mean Balanced AUC") + ggtitle("L1K Cross-validated Mean Balanced AUC") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
dev.off()

### mean compound AUC - avoids sampling
L1KPertMeanAUCML <- 1 - sapply(mycells, FUN=function(x) 
  mean(sapply(xvalres[[x]]$mlRanks, FUN=function(y)
    mean(sapply(y, FUN=mean)))))

L1KPertMeanAUCCos <- 1 - sapply(mycells, FUN=function(x) 
  mean(sapply(xvalres[[x]]$cosRanks, FUN=function(y)
    mean(sapply(y, FUN=mean)))))

pdf(file=file.path(outdir, "fig2c_L1KMeanPertAUC.pdf"), width=8, height=6)
ggplot(data.frame(ML=L1KPertMeanAUCML, Cos=L1KPertMeanAUCCos, cell=names(L1KPertMeanAUCML)), 
       aes(x=Cos, y=ML, color=cell)) + geom_point() + xlim(c(0.6, 1)) + theme_minimal() + 
  ylim(c(0.6,1)) + geom_abline(intercept=0, slope=1, linetype=2, col="grey") + xlab("Cosine Mean AUC") +
  ylab("Metric Learning Mean AUC") + ggtitle("L1K Cross Validated Mean Compound AUC")
dev.off()

# FDR fractions - used balanced 
L1KPertBalFDRML <- data.frame(fdr01=sapply(mycells, FUN=function(x) 
  mean(sapply(xvalres[[x]]$mlRanks, FUN=function(y) 
    mean(sapply(seq(100), FUN=function(z) mean(p.adjust(unlist(balancedSample(y)), "fdr") < 0.01)))))), 
  fdr05=sapply(mycells, FUN=function(x) 
    mean(sapply(xvalres[[x]]$mlRanks, FUN=function(y) 
      mean(sapply(seq(100), FUN=function(z) mean(p.adjust(unlist(balancedSample(y)), "fdr") < 0.05)))))),
  fdr10=sapply(mycells, FUN=function(x) 
    mean(sapply(xvalres[[x]]$mlRanks, FUN=function(y) 
      mean(sapply(seq(100), FUN=function(z) mean(p.adjust(unlist(balancedSample(y)), "fdr") < 0.10))))))
)
  
L1KPertBalFDRCos <- data.frame(fdr01=sapply(mycells, FUN=function(x) 
  mean(sapply(xvalres[[x]]$cosRanks, FUN=function(y) 
    mean(sapply(seq(100), FUN=function(z) mean(p.adjust(unlist(balancedSample(y)), "fdr") < 0.01)))))), 
  fdr05=sapply(mycells, FUN=function(x) 
    mean(sapply(xvalres[[x]]$cosRanks, FUN=function(y) 
      mean(sapply(seq(100), FUN=function(z) mean(p.adjust(unlist(balancedSample(y)), "fdr") < 0.05)))))),
  fdr10=sapply(mycells, FUN=function(x) 
    mean(sapply(xvalres[[x]]$cosRanks, FUN=function(y) 
      mean(sapply(seq(100), FUN=function(z) mean(p.adjust(unlist(balancedSample(y)), "fdr") < 0.10))))))
)

L1KPertBalFDRCos$method <- "cosine"
L1KPertBalFDRML$method <- "metric learning"
L1KPertBalFDRCos$cell <- row.names(L1KPertBalFDRCos)
L1KPertBalFDRML$cell <- row.names(L1KPertBalFDRML)

pdf(file=file.path(outdir, "fig2d_L1Kpertfdr05.pdf"), width=8, height=6)
ggplot(rbind(L1KPertBalFDRCos, L1KPertBalFDRML), aes(x=cell, y=fdr05, fill=method)) + 
  geom_bar(stat="identity", position="dodge") + theme_minimal() + ylab("Balanced FDR < 0.05") + 
  ggtitle("L1K Replicate Balanced FDR < 0.05") + ylim(c(0, 0.6))
dev.off()

# signal-to-noise

L1KsnrML <- sapply(xvalres, FUN=function(x) 
  sapply(x$inProdReps, FUN=function(y) 
    (sapply(y$same, mean) - mean(y$diff))/sd(y$diff)))
L1KsnrCos <- sapply(xvalres, FUN=function(x) 
  sapply(x$cosReps, FUN=function(y) 
    (sapply(y$same, mean) - mean(y$diff))/sd(y$diff)))

#bin <- hexbin(unlist(L1KsnrCos), unlist(L1KsnrML), xbins=50)
#my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
#plot(bin, colramp = my_colors, xlab="Cosine Signal-to-Noise ratio", ylab="Metric Learning SNR", 
#     main="Signal-to-noise ratio comparison by compound, L1000")

pdf(file=file.path(outdir, "Sfig_L1KSNR_density.pdf"), width=10, height=10)
ggplot(data.frame(cos=unlist(L1KsnrCos), ml=unlist(L1KsnrML)), aes(x=cos, y=ml)) + geom_bin_2d(bins=100) + 
  scale_fill_continuous(type="viridis") + ylim(c(-5, 10)) + xlim(c(-5, 10)) + theme_minimal() + 
  geom_abline(intercept=0, slope=1, col="red", linetype=2) + xlab("Cosine Signal-to-noise ratio") + 
  ylab("Metric Learning SNR") + ggtitle("Signal-to-noise ratio comparison by compound, L1000")
dev.off()

pdf(file=file.path(outdir, "Sfig_L1KSNR_pdf.pdf"), width=8, height=6)  
plot(density(unlist(L1KsnrCos), bw=0.03), col="blue", lwd=2, xlab="SNR", ylab="Density", 
     main="L1K compound signal-to-noise ratio")
lines(density(unlist(L1KsnrML), bw=0.03), col="red", lwd=2)
legend(x="topright", legend=c("Cosine", "Metric Learning"), lwd=4, col=c("blue", "red"))
dev.off()

pdf(file=file.path(outdir, "Sfig_L1KSNR_ecdf.pdf"), width=8, height=6)
plot(ecdf(unlist(L1KsnrCos)), col="blue", lwd=2, xlab="SNR", ylab="Density", 
     main="L1K compound SNR ECDF", xlim=c(-2, 5))
lines(ecdf(unlist(L1KsnrML)), col="red", lwd=2)
legend(x="bottomright", legend=c("Cosine", "Metric Learning"), lwd=4, col=c("blue", "red"))
grid()
dev.off()


# mean average precision? (Moshkov)
# folds of enrichment (Moshkov)



#### Fig 2e - summary PCL recall (auRank?)
pclds <- readRDS("data/pclds.rds")

pclcells <- sapply(strsplit(sapply(strsplit(list.files(l1kdir, pattern="L1Kmetric"), "cell="), 
                                  FUN=function(x) x[[2]]), ".rds"), FUN=function(x) x[[1]])
pclcells <- setdiff(pclcells, "all")
pclres <- list()

for (ii in seq_along(pclcells)){
  acell <- pclcells[ii]
  print(acell)
  modelds <- readRDS(file.path(l1kdir, list.files(l1kdir, pattern=sprintf("L1Kmetric.*%s.rds", acell))))
  f <- list.files(file.path(l1kdir, "models"), pattern=sprintf("L1Kmetric.*%s", acell))
  
  pclres[[acell]] <- getL1KMoA(datapath, l1kmeta, acell, mymodel=file.path(l1kdir, "models", f[1]), pclds=pclds)
}

L1KMOABalAUCML <- 1 - sapply(pclcells, FUN=function(x) 
    mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(pclres[[x]]$mlRanks, 1000))))))

L1KMOABalAUCCos <- 1 - sapply(pclcells, FUN=function(x) 
    mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(pclres[[x]]$cosRanks, 1000))))))

pdf(file.path(outdir, "fig2e_L1KMoAMeanBalAUC.pdf"), width=8, height=6)
ggplot(rbind(data.frame(method="ml", auc=L1KMOABalAUCML, cell=names(L1KMOABalAUCML)), 
      data.frame(method="cos", auc=L1KMOABalAUCCos, cell=names(L1KMOABalAUCCos))),
      aes(x=cell, y=auc, fill=method)) + geom_bar(stat="identity", position="dodge") + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  coord_cartesian(ylim=c(0.5,0.8)) + ylab("MoA Mean Balanced AUC") + ggtitle("L1K MoA Mean Balanced AUC")
dev.off()

pdf(file.path(outdir, "fig2e_L1KMoAMeanBalAUC_scatter.pdf"), width=8, height=6)
ggplot(data.frame(cell=names(L1KMOABalAUCML), cosAUC=L1KMOABalAUCCos, mlAUC=L1KMOABalAUCML), 
       aes(x=cosAUC, y=mlAUC, color=cell)) + geom_point(aes(shape=cell, stroke=1.5, size=1.5)) + xlim(c(0.5, 0.8)) + ylim(c(0.5, 0.8)) + 
  theme_minimal() + ggtitle("L1K MoA Mean Balanced AUC") + geom_abline(intercept=0, slope=1, linetype=2, col="black") + 
  scale_shape_manual(values=rep(c(3, 4, 16, 17, 18),5)) + 
  scale_color_manual(values=rep(cbbPalette, 3)[1:15]) + guides(size="none")
dev.off()


# Fig 2d - PCL recall example


# PCL SNR  

L1KMoASNRML <- sapply(pclres, FUN=function(x) (sapply(x$mlPCLs$setSims, mean) - mean(x$mlPCLs$allSims))/sd(x$mlPCLs$allSims))
L1KMoASNRCos <- sapply(pclres, FUN=function(x) (sapply(x$cosPCLs$setSims, mean) - mean(x$cosPCLs$allSims))/sd(x$cosPCLs$allSims))

# Need to fix the size of the points
pdf(file.path(outdir, "Sfig_L1KSNR_MoA_scatter.pdf"), width=8, height=6)
ggplot(data.frame(ml=unlist(L1KMoASNRML), cos=unlist(L1KMoASNRCos), 
        cell=unlist(sapply(seq(15), FUN=function(x) rep(names(L1KMoASNRML)[x], sapply(L1KMoASNRML, length)[x])))), 
       aes(x=cos, y=ml, color=cell)) + geom_point(aes(shape=cell, size=0.001, alpha=0.8)) + theme_minimal() + 
  scale_color_manual(values=sort(rep(cbbPalette, 3))[1:15]) + 
  scale_shape_manual(values=rep(c(16, 17, 18),5)) + guides(size="none", alpha="none") + 
  geom_abline(intercept=0, slope=1, col="red", linetype=2) + xlim(c(-2.5, 10)) + ylim(c(-2.5, 10)) + 
  ylab("Cosine MoA SNR") + ylab("Metric Learning MoA SNR") + ggtitle("Comparison of MoA SNR")
dev.off()

pdf(file.path(outdir, "Sfig_L1KSNR_MoA_DeltaDensity.pdf"), width=8, height=6)
ggplot(data.frame(ml=unlist(L1KMoASNRML), cos=unlist(L1KMoASNRCos), 
                  cell=unlist(sapply(seq(15), FUN=function(x) rep(names(L1KMoASNRML)[x], sapply(L1KMoASNRML, length)[x])))), 
       aes(x=ml-cos, color=cell, fill=cell)) + geom_density(alpha=0.6, position="stack") + theme_minimal() + 
  geom_vline(xintercept=0, linetype=2, col="black") + xlab("MoA SNR: Metric Learning - Cosine") +
  ggtitle("Change in MoA Mean SNR with metric learning vs cosine")
dev.off()

# Fig 2f - FDR plots for replicates, PCLs
L1KMoABalFDRML <- data.frame(fdr01=sapply(pclres, FUN=function(x) 
    mean(sapply(seq(10), FUN=function(y) mean(p.adjust(unlist(balancedSample(x$mlRanks, k=1000)), "fdr") < 0.01)))), 
  fdr05=sapply(pclres, FUN=function(x) 
    mean(sapply(seq(10), FUN=function(y) mean(p.adjust(unlist(balancedSample(x$mlRanks, k=1000)), "fdr") < 0.05)))), 
  fdr10=sapply(pclres, FUN=function(x) 
    mean(sapply(seq(10), FUN=function(y) mean(p.adjust(unlist(balancedSample(x$mlRanks, k=1000)), "fdr") < 0.10)))) 
)

L1KMoABalFDRCos <- data.frame(fdr01=sapply(pclres, FUN=function(x) 
  mean(sapply(seq(10), FUN=function(y) mean(p.adjust(unlist(balancedSample(x$cosRanks, k=1000)), "fdr") < 0.01)))), 
  fdr05=sapply(pclres, FUN=function(x) 
    mean(sapply(seq(10), FUN=function(y) mean(p.adjust(unlist(balancedSample(x$cosRanks, k=1000)), "fdr") < 0.05)))), 
  fdr10=sapply(pclres, FUN=function(x) 
    mean(sapply(seq(10), FUN=function(y) mean(p.adjust(unlist(balancedSample(x$cosRanks, k=1000)), "fdr") < 0.10)))) 
)

L1KMoABalFDRML$cell <- row.names(L1KMoABalFDRML)
L1KMoABalFDRML$method <- "metric learning"
L1KMoABalFDRCos$cell <- row.names(L1KMoABalFDRCos)
L1KMoABalFDRCos$method <- "cosine"

pdf(file=file.path(outdir, "fig2f_L1KMoAfdr10.pdf"), width=8, height=6)
ggplot(rbind(L1KMoABalFDRCos, L1KMoABalFDRML), aes(x=cell, y=fdr10, fill=method)) + 
  geom_bar(stat="identity", position="dodge") + theme_minimal() + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  ylab("Balanced FDR < 0.1") + ggtitle("L1K PCL Balanced FDR < 0.1") + ylim(c(0, 0.5))
dev.off()

pdf(file=file.path(outdir, "fig2f_L1KMoAfdr05.pdf"), width=8, height=6)
ggplot(rbind(L1KMoABalFDRCos, L1KMoABalFDRML), aes(x=cell, y=fdr05, fill=method)) + 
  geom_bar(stat="identity", position="dodge") + theme_minimal() + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  ylab("Balanced FDR < 0.05") + ggtitle("L1K PCL Balanced FDR < 0.05") + ylim(c(0, 0.5))
dev.off()


#### Figure 3: Cell Painting performance ####

# Fig 3a - cosine similarity example
# Fig 3b - cross validated replicate recall curves
# Fig 3c - summary cross validated replicate recall 
# Fig 3d - PCL recall example
# Fig 3e - summary PCL recall (auRank)
# Fig 3f - FDR plots for replicates, PCLs

braydir <- file.path(cpdir, "bray/models")
pclds <- readRDS("data/pclds.rds")

brayds <- loadBrayData(braypath)
lincsds1 <- loadLincsData(file.path(lincspath, "2016_04_01_a549_48hr_batch1_dmso_spherized_profiles_with_input_normalized_by_dmso.rds"), byCell = 1, splitGrps = 0)
lincsds2 <- loadLincsData(file.path(lincspath, "2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_dmso.rds"), byCell = 0, splitGrps = 0)

lincsds2_hidose$ds <- lincsds2$ds[lincsds2$metads$Metadata_mmoles_per_liter >= 10, ]
lincsds2_hidose$metads <- lincsds2$metads[lincsds2$metads$Metadata_mmoles_per_liter >= 10, ]

l1m <- table(lincsds1$metads$Metadata_moa)
l2m <- table(lincsds2$metads$Metadata_moa)

# The first element is "", i.e. empty and unassigned, so exclude this.
# It is very important that lincsds2 be read with byCell = 0 and splitGrps = 0, otherwise cell IDs and groupIDs get appended to pert identifiers. 
# The by Cell and splitGrps flags are important for training the metric and defining batch sizes, not for evaluation. 
lincsMoADS <- list(moaNames = c(names(l1m)[-1], names(l2m)[-1]), 
                   pertids = c(sapply(names(l1m)[-1], FUN=function(x) unique(lincsds1$metads$Metadata_pert_id[lincsds1$metads$Metadata_moa == x])), 
                               sapply(names(l2m)[-1], FUN=function(x) unique(lincsds2$metads$Metadata_pert_id[lincsds2$metads$Metadata_moa == x]))), 
                   pertnames = c(sapply(names(l1m)[-1], FUN=function(x) unique(lincsds1$metads$Metadata_pert_iname[lincsds1$metads$Metadata_moa == x])), 
                                 sapply(names(l2m)[-1], FUN=function(x) unique(lincsds2$metads$Metadata_pert_iname[lincsds2$metads$Metadata_moa == x]))))
                   

braymodel <- torch::torch_load(file.path(braydir, "braymetric_epch=10_smp=50_model.pt"))
lincs1model <- torch::torch_load(file.path(cpdir, "lincs/models/lincsmetric_a549_2016_dmso_epch=10_model.pt"))
lincs2model <- torch::torch_load(file.path(cpdir, "lincs/models/lincsmetric_batch2_dmso_bycell_epch=10_model.pt"))


brayReps <- getCPReps(brayds$ds, mymodel=braymodel, pertlabels=brayds$metads$Metadata_pert_id, datalabel="Bray DS centered full")
lincs1Reps <- getCPReps(lincsds1$ds, mymodel=lincs1model, pertlabels=lincsds1$metads$Metadata_pert_id, datalabel="LINCS A549 DMSO full model")
lincs2Reps <- getCPReps(lincsds2$ds, mymodel=lincs2model, pertlabels=lincsds2$metads$Metadata_pert_id, datalabel="LINCS Batch2 DMSO full model")

# These MoA results use the L1000 PCLs (pclds) as the definitive labels for what constitutes a MoA class
brayMoA <- getCPMoA(brayds$ds, brayds$metads, mymodel=braymodel, pertlabels=brayds$metads$Metadata_pert_id, 
                    pclds=pclds, datalabel="Bray MoA, centered full")
lincs1MoA <- getCPMoA(lincsds1$ds, lincsds1$metads, mymodel=lincs1model, pertlabels=lincsds1$metads$Metadata_pert_id, 
                      pclds=pclds, datalabel="LINCS A549 MoA DMSO full model")
lincs2MoA <- getCPMoA(lincsds2$ds, lincsds2$metads, mymodel=lincs2model, pertlabels=lincsds2$metads$Metadata_pert_id, 
                      pclds=pclds, datalabel="LINCS Batch2 MoA DMSO full model")
lincs2MoAByCell <- getCPMoA(lincsds2$ds, lincsds2$metads, mymodel=lincs2model, pertlabels=lincsds2$metads$Metadata_pert_id, 
                           pclds=pclds, datalabel="LINCS Batch2 MoA DMSO full model")

lincs2MoAByCellHighDose <- getCPMoA(lincsds2_hidose$ds, lincsds2_hidose$metads, mymodel=lincs2model, pertlabels=lincsds2_hidose$metads$Metadata_pert_id, 
                            pclds=pclds, datalabel="LINCS Batch2 MoA DMSO full model High Dose")

# The 'Native' MoA results use mechanism of action annotations directly from the LINCS Cell Painting datasets.
brayMoANative <- getCPMoA(brayds$ds, brayds$metads, mymodel=braymodel, pertlabels=brayds$metads$Metadata_pert_id, 
                    pclds=lincsMoADS, datalabel="Bray Native MoA, centered full")
lincs1MoANative <- getCPMoA(lincsds1$ds, lincsds1$metads, mymodel=lincs1model, pertlabels=lincsds1$metads$Metadata_pert_id, 
                       pclds=lincsMoADS, datalabel="LINCS A549 Native MoAs DMSO full model")
lincs2MoABCNative <- getCPMoA(lincsds2$ds, lincsds2$metads, mymodel=lincs2model, pertlabels=lincsds2$metads$Metadata_pert_id, 
                           pclds=lincsMoADS, datalabel="LINCS Batch2 Native MoAs DMSO full model by cell")


# Fig 3a - Cell Painting Cosine

pdf(file=file.path(outdir, "fig3a_bray_simpdfs.pdf"), width=8, height=6)
plot(density(unlist(brayReps$inProdReps$diff), bw=0.01), col="blue", lwd=1.5, xlim=c(-1, 1), 
     xlab="Similarity", ylab="Density", main=sprintf("CDRP similarity density"))
lines(density(unlist(brayReps$inProdReps$same), bw=0.01), col="firebrick3", lwd=1.5)
lines(density(unlist(brayReps$cosReps$diff), bw=0.01), col="forestgreen", lwd=1.5, lty=5)
lines(density(unlist(brayReps$cosReps$same), bw=0.01), col="orange", lwd=1.5, lty=5)
lines(c(0, 0), c(0, 100), col="grey", lty=2)
legend(x="topright", legend=c("Metric learning replicates", "Metric Learning all pairs", 
                              "Cosine replicates", "Cosine all pairs"), lwd=2.5, col=c("firebrick3", "blue", "orange", "forestgreen"), 
       lty=c(1,1,5,5))
dev.off()


L1Perts <- intersect(names(lincs1Reps$inProdReps$same), unique(lincsds1$metads$Metadata_pert_id[lincsds1$metads$Metadata_pert_type == "trt"]))
L2Perts <- intersect(names(lincs2Reps$inProdReps$same), unique(lincsds2$metads$Metadata_pert_id[lincsds2$metads$Metadata_pert_type == "trt"]))

pdf(file=file.path(outdir, "fig3a2_lincs1_simpdfs.pdf"), width=8, height=6)
plot(density(unlist(lincs1Reps$cosReps$diff), bw=0.01), col="forestgreen", lwd=1.5, lty=5, xlab="Similarity", ylab="Density", main=sprintf("Cell Health A549 similarity density"))
lines(density(unlist(lincs1Reps$cosReps$same[L1Perts]), bw=0.01), col="orange", lwd=1.5, lty=5)
lines(density(unlist(lincs1Reps$inProdReps$diff), bw=0.01), col="blue", lwd=1.5)
lines(density(unlist(lincs1Reps$inProdReps$same[L1Perts]), bw=0.01), col="firebrick3", lwd=1.5)
legend(x="topright", legend=c("Metric learning replicates", "Metric Learning all pairs", 
                              "Cosine replicates", "Cosine all pairs"), lwd=2.5, col=c("firebrick3", "blue", "orange", "forestgreen"), 
       lty=c(1,1,5,5))
dev.off()

pdf(file=file.path(outdir, "fig3a3_lincs2_simpdfs.pdf"), width=8, height=6)
plot(density(unlist(lincs2Reps$cosReps$diff), bw=0.01), col="forestgreen", lwd=1.5, lty=5, xlab="Similarity", ylab="Density", main=sprintf("Cell Health Batch 2 similarity density"))
lines(density(unlist(lincs2Reps$cosReps$same[L2Perts]), bw=0.01), col="orange", lwd=1.5, lty=5)
lines(density(unlist(lincs2Reps$inProdReps$diff), bw=0.01), col="blue", lwd=1.5)
lines(density(unlist(lincs2Reps$inProdReps$same[L2Perts]), bw=0.01), col="firebrick3", lwd=1.5)
legend(x="topright", legend=c("Metric learning replicates", "Metric Learning all pairs", 
                              "Cosine replicates", "Cosine all pairs"), lwd=2.5, col=c("firebrick3", "blue", "orange", "forestgreen"), 
       lty=c(1,1,5,5))
dev.off()


# Fig 3b: Cross validation


# Fig 3d, 3e: MoA

pdf(file=file.path(outdir, "fig3_MoABrayAUC.pdf"), width=8, height=6)
plot(ecdf(unlist(balancedSample(brayMoA$mlRanks, k=500))), col="red", xlim=c(0,1), xlab="Balanced Rank", ylab="Cumulative Fraction", 
     main="CDRP MoA ROC", lwd=2)
lines(ecdf(unlist(balancedSample(brayMoA$cosRanks, k=500))), col="blue", lwd=2)
legend(x="bottomright", legend=c(sprintf("Metric Learning, AUC = %0.4f", 1-mean(sapply(seq(100), FUN=function(x) mean(unlist(balancedSample(brayMoA$mlRanks, k=500)))))), 
                                 sprintf("Cosine, AUC = %0.4f", 1-mean(sapply(seq(100), FUN=function(x) mean(unlist(balancedSample(brayMoA$cosRanks, k=500))))))),
       col=c("red", "blue"), lwd=2.5)
dev.off()

pdf(file=file.path(outdir, "fig3_MoALincs1AUC.pdf"), width=8, height=6)
plot(ecdf(unlist(balancedSample(lincs1MoA$mlRanks, k=500))), col="red", xlim=c(0,1), xlab="Balanced Rank", ylab="Cumulative Fraction", 
     main="Cell Health 1 MoA ROC", lwd=2)
lines(ecdf(unlist(balancedSample(lincs1MoA$cosRanks, k=500))), col="blue", lwd=2)
legend(x="bottomright", legend=c(sprintf("Metric Learning, AUC = %0.4f", 1-mean(sapply(seq(100), FUN=function(x) mean(unlist(balancedSample(lincs1MoA$mlRanks, k=500)))))), 
                                 sprintf("Cosine, AUC = %0.4f", 1-mean(sapply(seq(100), FUN=function(x) mean(unlist(balancedSample(lincs1MoA$cosRanks, k=500))))))),
       col=c("red", "blue"), lwd=2.5)
dev.off()

pdf(file=file.path(outdir, "fig3_MoALincs2AUC.pdf"), width=8, height=6)
plot(ecdf(unlist(balancedSample(lincs2MoAByCell$A549$mlRanks, k=500))), col=cbbPalette[2], xlim=c(0,1), xlab="Balanced Rank", ylab="Cumulative Fraction", 
     main="Cell Health 2 MoA ROC")
lines(ecdf(unlist(balancedSample(lincs2MoAByCell$A549$cosRanks, k=500))), col=cbbPalette[3], lwd=2)
lines(ecdf(unlist(balancedSample(lincs2MoAByCell$MCF7$mlRanks, k=500))), col=cbbPalette[4], lwd=2)
lines(ecdf(unlist(balancedSample(lincs2MoAByCell$MCF7$cosRanks, k=500))), col=cbbPalette[5], lwd=2)
lines(ecdf(unlist(balancedSample(lincs2MoAByCell$U2OS$mlRanks, k=500))), col=cbbPalette[6], lwd=2)
lines(ecdf(unlist(balancedSample(lincs2MoAByCell$U2OS$cosRanks, k=500))), col=cbbPalette[7], lwd=2)
legend(x="bottomright", legend=c(sprintf("A549 Metric Learning, AUC = %0.4f", 1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$A549$mlRanks, k=500)))))), 
                                 sprintf("A549 Cosine, AUC = %0.4f", 1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$A549$cosRanks, k=500)))))),
                                 sprintf("MCF7 Metric Learning, AUC = %0.4f", 1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$MCF7$mlRanks, k=500)))))), 
                                 sprintf("MCF7 Cosine, AUC = %0.4f", 1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$MCF7$cosRanks, k=500)))))),
                                 sprintf("U2OS Metric Learning, AUC = %0.4f", 1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$U2OS$mlRanks, k=500)))))), 
                                 sprintf("U2OS Cosine, AUC = %0.4f", 1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$U2OS$cosRanks, k=500))))))),
       col=cbbPalette[2:7], lwd=2.5)
dev.off()

pdf(file=file.path(outdir, "fig3_MoALincs2AUCHighDose.pdf"), width=8, height=6)
plot(ecdf(unlist(balancedSample(lincs2MoAByCellHighDose$A549$mlRanks, k=500))), col=cbbPalette[2], xlim=c(0,1), xlab="Balanced Rank", ylab="Cumulative Fraction", 
     main="Cell Health 2 MoA ROC")
lines(ecdf(unlist(balancedSample(lincs2MoAByCellHighDose$A549$cosRanks, k=500))), col=cbbPalette[3], lwd=2)
lines(ecdf(unlist(balancedSample(lincs2MoAByCellHighDose$MCF7$mlRanks, k=500))), col=cbbPalette[4], lwd=2)
lines(ecdf(unlist(balancedSample(lincs2MoAByCellHighDose$MCF7$cosRanks, k=500))), col=cbbPalette[5], lwd=2)
lines(ecdf(unlist(balancedSample(lincs2MoAByCellHighDose$U2OS$mlRanks, k=500))), col=cbbPalette[6], lwd=2)
lines(ecdf(unlist(balancedSample(lincs2MoAByCellHighDose$U2OS$cosRanks, k=500))), col=cbbPalette[7], lwd=2)
legend(x="bottomright", legend=c(sprintf("A549 Metric Learning, AUC = %0.4f", 1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$A549$mlRanks, k=500)))))), 
                                 sprintf("A549 Cosine, AUC = %0.4f", 1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$A549$cosRanks, k=500)))))),
                                 sprintf("MCF7 Metric Learning, AUC = %0.4f", 1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$MCF7$mlRanks, k=500)))))), 
                                 sprintf("MCF7 Cosine, AUC = %0.4f", 1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$MCF7$cosRanks, k=500)))))),
                                 sprintf("U2OS Metric Learning, AUC = %0.4f", 1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$U2OS$mlRanks, k=500)))))), 
                                 sprintf("U2OS Cosine, AUC = %0.4f", 1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$U2OS$cosRanks, k=500))))))),
       col=cbbPalette[2:7], lwd=2.5)
dev.off()

##### auRankCPdf #####
auRankCPdf <- data.frame(dataset=character(), auROC=numeric(), method=character(), fdr10=numeric())
auRankCPdf <- rbind(auRankCPdf, data.frame(dataset=c("CDRP MoA", "CDRP MoA"), 
                                           method=c("ML", "Cos"), 
                                           auROC=c(1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(brayMoA$mlRanks, k=500))))),
                                                   1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(brayMoA$cosRanks, k=500)))))), 
                                           fdr10=c(mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(brayMoA$mlRanks, k=500)), method="fdr") < 0.1))),
                                                   mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(brayMoA$cosRanks, k=500)), method="fdr") < 0.1))))))
auRankCPdf <- rbind(auRankCPdf, data.frame(dataset=c("Lincs1 MoA", "Lincs1 MoA"), 
                                           method=c("ML", "Cos"), 
                                           auROC=c(1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs1MoA$mlRanks, k=500))))),
                                                   1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs1MoA$cosRanks, k=500)))))), 
                                           fdr10=c(mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs1MoA$mlRanks, k=500)), method="fdr") < 0.1))),
                                                   mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs1MoA$cosRanks, k=500)), method="fdr") < 0.1))))))
auRankCPdf <- rbind(auRankCPdf, data.frame(dataset=c("Lincs1 Native MoA", "Lincs1 Native MoA"), 
                                           method=c("ML", "Cos"), 
                                           auROC=c(1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs1MoANative$mlRanks, k=500))))),
                                                   1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs1MoANative$cosRanks, k=500)))))),
                                           fdr10=c(mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs1MoANative$mlRanks, k=500)), method="fdr") < 0.1))),
                                                   mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs1MoANative$cosRanks, k=500)), method="fdr") < 0.1))))))
auRankCPdf <- rbind(auRankCPdf, data.frame(dataset=c("Lincs2 MoA A549", "Lincs2 MoA A549"), 
                                           method=c("ML", "Cos"), 
                                           auROC=c(1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$A549$mlRanks, k=500))))),
                                                   1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$A549$cosRanks, k=500)))))),
                                           fdr10=c(mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoAByCell$A549$mlRanks, k=500)), method="fdr") < 0.1))),
                                                   mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoAByCell$A549$cosRanks, k=500)), method="fdr") < 0.1))))))
auRankCPdf <- rbind(auRankCPdf, data.frame(dataset=c("Lincs2 MoA MCF7", "Lincs2 MoA MCF7"), 
                                           method=c("ML", "Cos"), 
                                           auROC=c(1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$MCF7$mlRanks, k=500))))),
                                                   1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$MCF7$cosRanks, k=500)))))), 
                                           fdr10=c(mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoAByCell$MCF7$mlRanks, k=500)), method="fdr") < 0.1))),
                                                   mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoAByCell$MCF7$cosRanks, k=500)), method="fdr") < 0.1))))))
auRankCPdf <- rbind(auRankCPdf, data.frame(dataset=c("Lincs2 MoA U2OS", "Lincs2 MoA U2OS"), 
                                           method=c("ML", "Cos"), 
                                           auROC=c(1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$U2OS$mlRanks, k=500))))),
                                                   1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCell$U2OS$cosRanks, k=500)))))),
                                           fdr10=c(mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoAByCell$U2OS$mlRanks, k=500)), method="fdr") < 0.1))),
                                                   mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoAByCell$U2OS$cosRanks, k=500)), method="fdr") < 0.1))))))
auRankCPdf <- rbind(auRankCPdf, data.frame(dataset=c("Lincs2 MoA - High Dose A549", "Lincs2 MoA - High Dose A549"), 
                                           method=c("ML", "Cos"), 
                                           auROC=c(1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCellHighDose$A549$mlRanks, k=500))))),
                                                   1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCellHighDose$A549$cosRanks, k=500)))))), 
                                           fdr10=c(mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoAByCellHighDose$A549$mlRanks, k=500)), method="fdr") < 0.1))),
                                                   mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoAByCellHighDose$A549$cosRanks, k=500)), method="fdr") < 0.1))))))
auRankCPdf <- rbind(auRankCPdf, data.frame(dataset=c("Lincs2 MoA - High Dose MCF7", "Lincs2 MoA - High Dose MCF7"), 
                                           method=c("ML", "Cos"), 
                                           auROC=c(1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCellHighDose$MCF7$mlRanks, k=500))))),
                                                   1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCellHighDose$MCF7$cosRanks, k=500)))))), 
                                           fdr10=c(mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoAByCellHighDose$MCF7$mlRanks, k=500)), method="fdr") < 0.1))),
                                                   mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoAByCellHighDose$MCF7$cosRanks, k=500)), method="fdr") < 0.1))))))
auRankCPdf <- rbind(auRankCPdf, data.frame(dataset=c("Lincs2 MoA - High Dose U2OS", "Lincs2 MoA - High Dose U2OS"), 
                                           method=c("ML", "Cos"), 
                                           auROC=c(1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCellHighDose$U2OS$mlRanks, k=500))))),
                                                   1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoAByCellHighDose$U2OS$cosRanks, k=500)))))),
                                           fdr10=c(mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoAByCellHighDose$U2OS$mlRanks, k=500)), method="fdr") < 0.1))),
                                                   mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoAByCellHighDose$U2OS$cosRanks, k=500)), method="fdr") < 0.1))))))
auRankCPdf <- rbind(auRankCPdf, data.frame(dataset=c("Lincs2 MoA - Native A549", "Lincs2 MoA - Native A549"), 
                                           method=c("ML", "Cos"), 
                                           auROC=c(1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoABCNative$A549$mlRanks, k=500))))),
                                                   1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoABCNative$A549$cosRanks, k=500)))))), 
                                           fdr10=c(mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoABCNative$A549$mlRanks, k=500)), method="fdr") < 0.1))),
                                                   mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoABCNative$A549$cosRanks, k=500)), method="fdr") < 0.1))))))
auRankCPdf <- rbind(auRankCPdf, data.frame(dataset=c("Lincs2 MoA - Native MCF7", "Lincs2 MoA - Native MCF7"), 
                                           method=c("ML", "Cos"), 
                                           auROC=c(1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoABCNative$MCF7$mlRanks, k=500))))),
                                                   1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoABCNative$MCF7$cosRanks, k=500)))))), 
                                           fdr10=c(mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoABCNative$MCF7$mlRanks, k=500)), method="fdr") < 0.1))),
                                                   mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoABCNative$MCF7$cosRanks, k=500)), method="fdr") < 0.1))))))
auRankCPdf <- rbind(auRankCPdf, data.frame(dataset=c("Lincs2 MoA - Native U2OS", "Lincs2 MoA - Native U2OS"), 
                                           method=c("ML", "Cos"), 
                                           auROC=c(1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoABCNative$U2OS$mlRanks, k=500))))),
                                                   1-mean(sapply(seq(10), FUN=function(x) mean(unlist(balancedSample(lincs2MoABCNative$U2OS$cosRanks, k=500)))))), 
                                           fdr10=c(mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoABCNative$U2OS$mlRanks, k=500)), method="fdr") < 0.1))),
                                                   mean(sapply(seq(10), FUN=function(x) mean(p.adjust(unlist(balancedSample(lincs2MoABCNative$U2OS$cosRanks, k=500)), method="fdr") < 0.1))))))

##### auRankCPdf Plot #####
pdf(file.path(outdir, "fig3S_CPMoAAUROC_bar.pdf"), width=8, height=6)
ggplot(auRankCPdf, aes(x=dataset, y=auROC, fill=method)) + geom_bar(stat="identity", position="dodge") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + xlab("Cell Painting dataset") + ylab("AUROC") + 
  ggtitle("Mechanism of Action AUROC for cell painting datasets") + coord_cartesian(ylim=c(0.4,1))
dev.off()

pdf(file.path(outdir, "fig3S_CPMoAfdr10_bar.pdf"), width=8, height=6)
ggplot(auRankCPdf, aes(x=dataset, y=fdr10, fill=method)) + geom_bar(stat="identity", position="dodge") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + xlab("Cell Painting dataset") + ylab("Balanced FDR < 0.1") + 
  ggtitle("Mechanism of Action FDR < 0.1 for cell painting datasets")
dev.off()

#### MoA Scatter for Cell Painting Datasets, Mean Rank and Mean FDR < 0.1
L2MoA <- data.frame(cell=c(rep("MCF7", 88), rep("U2OS", 88), rep("A549", 88)), 
                    mlL2vals = as.numeric(sapply(lincs2MoAByCell, FUN=function(x) sapply(x$mlRanks, mean))), 
                    cosL2vals = as.numeric(sapply(lincs2MoAByCell, FUN=function(x) sapply(x$cosRanks, mean))),
                    mlL2FDR = as.numeric(sapply(lincs2MoAByCell, FUN=function(x) sapply(x$mlRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.1)))), 
                    cosL2FDR = as.numeric(sapply(lincs2MoAByCell, FUN=function(x) sapply(x$cosRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.1)))), 
                    dose="all")

L2MoAHD <- data.frame(cell=c(rep("MCF7", 88), rep("U2OS", 88), rep("A549", 88)), 
                      mlL2vals = as.numeric(sapply(lincs2MoAByCellHighDose, FUN=function(x) sapply(x$mlRanks, mean))), 
                      cosL2vals = as.numeric(sapply(lincs2MoAByCellHighDose, FUN=function(x) sapply(x$cosRanks, mean))),
                      mlL2FDR = as.numeric(sapply(lincs2MoAByCellHighDose, FUN=function(x) sapply(x$mlRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.1)))), 
                      cosL2FDR = as.numeric(sapply(lincs2MoAByCellHighDose, FUN=function(x) sapply(x$cosRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.1)))), 
                      dose = "10 uM")

LMoA <- data.frame(cell="A549", 
                   mlL1vals = sapply(lincs1MoA$mlRanks, FUN=mean), 
                   cosL1vals = sapply(lincs1MoA$cosRanks, FUN=mean),
                   mlL1FDR = sapply(lincs1MoA$mlRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.1)), 
                   cosL1FDR = sapply(lincs1MoA$cosRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.1)))


BrayMoADF <- data.frame(cell="U2OS", 
                        mlBrayvals = sapply(brayMoA$mlRanks, FUN=mean), 
                        cosBrayvals = sapply(brayMoA$cosRanks, FUN=mean), 
                        mlBrayFDR = sapply(brayMoA$mlRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.1)), 
                        cosBrayFDR = sapply(brayMoA$cosRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.1)))

pdf(file.path(outdir, "fig3S_CPLincs2_MoAMean.pdf"), width=8, height=6)
ggplot(rbind(L2MoA, L2MoAHD), aes(x=cosL2vals, y=mlL2vals, color=cell, shape=dose)) + geom_point(size=3) + theme_minimal() + xlim(c(0,0.7)) + ylim(c(0,0.7)) + 
  ggtitle("Cell Health Batch 2 MoA Mean Rank") + geom_abline(intercept=0, slope=1, col="black", lty=2) + xlab("Cosine Mean Rank") + ylab("Metric Learning Mean Rank")
dev.off()

pdf(file.path(outdir, "fig3S_CPLincs2_MoaFDR.pdf"))
ggplot(rbind(L2MoA, L2MoAHD), aes(x=cosL2FDR, y=mlL2FDR, color=cell, shape=dose)) + geom_point(size=3) + theme_minimal() + xlim(c(0,1)) + ylim(c(0,1)) + 
  ggtitle("Cell Health Batch 2 MoA FDR < 0.1") + geom_abline(intercept=0, slope=1, col="black", lty=2) + xlab("Cosine Fraction FDR < 0.1") + ylab("Metric Learning Fraction < 0.1")
dev.off()

pdf(file.path(outdir, "fig3S_CPLincs1_MoAMean.pdf"), width=8, height=6)
ggplot(LMoA, aes(x=cosL1vals, y=mlL1vals)) + geom_point(size=3, color=cbbPalette[2]) + theme_minimal() + geom_abline(intercept=0, slope=1, col="red", lty=2) + xlab("Cosine Mean Rank") + 
  ylab("Metric Learning Mean Rank") + xlim(c(0, 0.8)) + ylim(c(0, 0.8)) + ggtitle("Cell Health Batch 1 MoA Mean Rank")
dev.off()

pdf(file.path(outdir, "fig3S_CPLincs1_MoAFDR.pdf"))
ggplot(LMoA, aes(x=cosL1FDR, y=mlL1FDR)) + geom_point(size=3, color=cbbPalette[2]) + theme_minimal() + xlim(c(0,1)) + ylim(c(0,1)) + 
  ggtitle("Cell Health Batch 1 MoA FDR < 0.1") + geom_abline(intercept=0, slope=1, col="black", lty=2) + xlab("Cosine Fraction FDR < 0.1") + ylab("Metric Learning Fraction < 0.1")
dev.off()

pdf(file.path(outdir, "fig3S_CPBray_MoAMean.pdf"), width=8, height=6)
ggplot(BrayMoADF, aes(x=cosBrayvals, y=mlBrayvals)) + geom_point(size=3, color="forestgreen") + theme_minimal() + xlim(c(0, 0.7)) + ylim(c(0, 0.7)) + 
  ggtitle("CDRP MoA Mean Rank") + xlab("Cosine Mean Rank") + ylab("Metric Learning Mean Rank") + geom_abline(intercept=0, slope=1, col="black", lty=2)
dev.off()

pdf(file.path(outdir, "fig3S_CPBray_MoAFDR.pdf"))
ggplot(BrayMoADF, aes(x=cosBrayFDR, y=mlBrayFDR)) + geom_point(size=3, color="forestgreen") + theme_minimal() + xlim(c(0,1)) + ylim(c(0,1)) + 
  ggtitle("CDRP MoA FDR < 0.1") + geom_abline(intercept=0, slope=1, col="black", lty=2) + xlab("Cosine Fraction FDR < 0.1") + ylab("Metric Learning Fraction < 0.1")
dev.off()



# Supp Fig: doses
pdf(file=file.path(outdir, "fig3S_CPDoses.pdf"), width=8, height=6)
plot(density(log10(brayds$metads$Metadata_mmoles_per_liter + 1e-3), bw=0.03), col=cbbPalette[2], lwd=3, xlab="Log10 Dose uM", ylab="Density", main="Cell Painting Doses, uM")
lines(density(log10(lincsds1$metads$Metadata_mmoles_per_liter + 1e-3), bw=0.03), col=cbbPalette[4], lwd=3)
lines(density(log10(lincsds2$metads$Metadata_mmoles_per_liter + 1e-3), bw=0.03), col=cbbPalette[5], lwd=3)
legend(x="topleft", legend=c("CDRP 2017", "Cell Health 1", "Cell Health 2"), col=cbbPalette[c(2, 4, 5)], lwd=3)
dev.off()


#### Figure 4: How much data do you need to learn a useful representation? ####
# L1000
ntraindir <- "~/Work/bhk/analysis/metric_learning/2022_L1K/modelRuns"
f <- list.files(ntraindir, pattern="ntraining")
fcells <- sapply(strsplit(f, split="cell="), FUN=function(x) strsplit(x[2], ".rds")[[1]])

ntrainL1K <- lapply(seq_along(f), FUN=function(x) 
  makeNTrainingPlots(ntrainingResPath = file.path(ntraindir, f[x]), saveFiles=1, fname=fcells[x], mytitle=fcells[x], outpath=file.path(ntraindir, "nTrainingFigs")))

resNorm <- bind_rows(lapply(seq(6), FUN=function(x) ntrainL1K[[x]]$resNorm))

pdf(file=file.path(outdir, "fig4a_normValidLoss.pdf"), width=8, height=6)
ggplot(resNorm, aes(x=log2(trainClasses+1), y=validAvgLoss, color=dsname)) + geom_point() + geom_smooth(method="loess") + theme_minimal() + 
  xlab("Log2 N of Training Compounds") + ylab("Normalized Average Validation Loss by compound") + ggtitle("Validation Loss vs Number of Training compounds") + 
  geom_hline(yintercept=-1, color="black", linetype="dashed") + scale_color_discrete(name="cell line")
dev.off()

pdf(file=file.path(outdir, "fig4Sa_normTrainLoss.pdf"), width=8, height=6)
ggplot(resNorm, aes(x=log2(trainClasses+1), y=trainAvgLoss, color=dsname)) + geom_point() + geom_smooth(method="loess") + theme_minimal() + 
  xlab("Log2 N of Training Compounds") + ylab("Normalized Average Train Loss by compound") + ggtitle("Training Loss vs Number of Training compounds") + 
  geom_hline(yintercept=-1, color="black", linetype="dashed") + scale_color_discrete(name="cell line") + ylim(-10, 0)
dev.off()

pdf(file=file.path(outdir, "fig4b_deltaValidAURank.pdf"), width=8, height=6)
ggplot(resNorm, aes(x=log2(trainClasses+1), y=validAURank, color=dsname)) + geom_point() + geom_smooth(method="loess") + theme_minimal() + 
  xlab("Log2 N of Training Compounds") + ylab("Delta Average Validation AURank by compound") + ggtitle("Change in Validation AURank vs Number of Training compounds") + 
  geom_hline(yintercept=0, color="black", linetype="dashed") + scale_color_discrete(name="cell line")
dev.off()


# Cell Painting
cpdir <- "~/Work/bhk/analysis/metric_learning/2022_cellpaint/"
braydir <- file.path(cpdir, "bray/models")
lincsdir <- file.path(cpdir, "lincs/models")

brayntrain <- makeNTrainingPlots(ntrainingResPath= file.path(braydir, "brayntraining_epch=10.rds"), saveFiles=1, fname="bray", mytitle="bray", 
                                 outpath=file.path(braydir, "ntrainFigs"))

lincsntrainf <- list.files(lincsdir, pattern="ntraining")
lincsCells <- c("Batch1", "Batch2")

lincsntrain <- lapply(seq_along(lincsntrainf), FUN=function(x) 
  makeNTrainingPlots(ntrainingResPath = file.path(lincsdir, lincsntrainf[x]), saveFiles=1, fname=lincsCells[x], mytitle=lincsCells[x], 
                     outpath=file.path(lincsdir, "ntrainFigs")))

cpntrain <- c(lincsntrain, list(brayntrain))
cpNorm <- bind_rows(lapply(seq_along(cpntrain), FUN=function(x) cpntrain[[x]]$resNorm))

pdf(file=file.path(outdir, "fig4c_normValidLossCP.pdf"), width=8, height=6)
ggplot(cpNorm, aes(x=log2(trainClasses+1), y=validAvgLoss, color=dsname)) + geom_point() + geom_smooth(method="loess") + theme_minimal() + 
  xlab("Log2 N of Training Compounds") + ylab("Normalized Average Validation Loss by compound") + ggtitle("Cell Painting Validation Loss vs Number of Training compounds") + 
  geom_hline(yintercept=-1, color="black", linetype="dashed") + scale_color_discrete(name="dataset")
dev.off()

pdf(file=file.path(outdir, "fig4Sc_normTrainLoss.pdf"), width=8, height=6)
ggplot(cpNorm, aes(x=log2(trainClasses+1), y=trainAvgLoss, color=dsname)) + geom_point() + geom_smooth(method="loess") + theme_minimal() + 
  xlab("Log2 N of Training Compounds") + ylab("Normalized Average Train Loss by compound") + ggtitle("Cell Painting Training Loss vs Number of Training compounds") + 
  geom_hline(yintercept=-1, color="black", linetype="dashed") + scale_color_discrete(name="dataset") 
dev.off()

pdf(file=file.path(outdir, "fig4d_deltaValidAURank.pdf"), width=8, height=6)
ggplot(cpNorm, aes(x=log2(trainClasses+1), y=validAURank, color=dsname)) + geom_point() + geom_smooth(method="loess") + theme_minimal() + 
  xlab("Log2 N of Training Compounds") + ylab("Delta Average Validation AURank by compound") + ggtitle("Cell Painting Change in Validation AURank vs Number of Training compounds") + 
  geom_hline(yintercept=0, color="black", linetype="dashed") + scale_color_discrete(name="dataset")
dev.off()


#### Figure 5: Cell line specificity ####




#### Figure 6: Biological Interpretation


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


 
#############################
########### Committee Figures
#############################


# L1K PCL losses
setwd("~/Work/bhk/analysis/metric_learning/2022_L1K/modelRuns/PCLFigs")
mycells <- c("A375", "A549", "ASC", "HCC515", "HEK293", "HEPG2", "MCF7", "NPC", "PC3", "BT20", "HUVEC")
pclf <- list.files(pattern="PCL_pcl_res.rds")

auRanks <- data.frame()

for (acell in mycells){
  a <- readRDS(sprintf("%s_PCL_pcl_res.rds", acell))
  
  mlbrk <- 1-mean(unlist(sapply(a$MLRanks, FUN=function(x) sample(x, min(1000, length(x))))))
  cosbrk <- 1-mean(unlist(sapply(a$cosranks, FUN=function(x) sample(x, min(1000, length(x))))))
  
  auRanks <- rbind(auRanks, data.frame(trainset="cell", testset=acell, 
                              auMLrank=(1-mean(unlist(a$MLRanks))), auCosRank=(1-mean(unlist(a$cosranks))), 
                              auMLBalrank=mlbrk, auCosBalRank=cosbrk))
  
  b <- readRDS(sprintf("all_%s_PCL_pcl_res.rds", acell))
  mlbrk <- 1-mean(unlist(sapply(b$MLRanks, FUN=function(x) sample(x, min(1000, length(x))))))
  cosbrk <- 1-mean(unlist(sapply(b$cosranks, FUN=function(x) sample(x, min(1000, length(x))))))
  
  auRanks <- rbind(auRanks, data.frame(trainset="all", testset=acell, 
                                       auMLrank=(1-mean(unlist(b$MLRanks))), auCosRank=(1-mean(unlist(b$cosranks))), 
                                       auMLBalrank=mlbrk, auCosBalRank=cosbrk))
}



pdf(file.path("../../figs", "L1K_PCL_deltaRank_barplot.pdf"), width=10, height=6)
ggplot(auRanks, aes(x=testset, y=auMLrank - auCosRank, fill=as.factor(trainset))) + geom_bar(position="dodge", stat="identity") + theme_minimal() +
  xlab("Cell Line") + ylab("ML Rank - cos Rank") + ggtitle("Delta Mean Rank, PCL data") + scale_fill_discrete(name = "Training Set")
dev.off()
pdf(file.path("../../figs", "L1K_PCL_MLRank_barplot.pdf"), width=10, height=6)
ggplot(auRanks, aes(x=testset, y=auMLrank, fill=as.factor(trainset))) + geom_bar(position="dodge", stat="identity") + theme_minimal() +
  xlab("Cell Line") + ylab("ML Rank") + ggtitle("Mean ML Rank, PCL data") + scale_fill_discrete(name = "Training Set")
dev.off()

auRanksCell <- rbind(auRanks[which(auRanks$trainset == "cell"),], auRanks[which(auRanks$trainset == "cell"),])
auRanksCell$metric <- c(rep("cosine", 11), rep("ML",11))
auRanksCell$rank <- 0
auRanksCell$rank[auRanksCell$metric == "cosine"] <- auRanksCell$auCosRank[auRanksCell$metric == "cosine"]
auRanksCell$rank[auRanksCell$metric == "ML"] <- auRanksCell$auMLrank[auRanksCell$metric == "ML"]

pdf(file.path("../../figs", "L1K_PCL_rankVsMetric_barplot.pdf"), width=10, height=6)
ggplot(auRanksCell, aes(x=testset, y=rank, fill=as.factor(metric))) + geom_bar(position="dodge", stat="identity") + theme_minimal() +
  xlab("Cell Line") + ylab("Mean Rank") + ggtitle("Mean metric Rank, PCL data") + scale_fill_discrete(name = "Training Set") + coord_cartesian(ylim=c(0.4,1)) + 
  geom_hline(yintercept=0.5, linetype="dashed", color="black")
dev.off()


a <- readRDS("HCC515_PCL_pcl_res.rds")

pdf(file.path("../../figs", "L1K_PCL_Zscore_HCC515.pdf"), width=8, height=6)
plot(density((unlist(a$cosSims) - mean(a$cosbg))/sd(a$cosbg), bw=0.05), col="blue", lwd=2, xlim=c(-5, 12), xlab="PCL Z-score", main="HCC515 PCL Z-scores")
lines(density((unlist(a$MLSims) - mean(a$MLbg))/sd(a$MLbg), bw=0.05), col="red", lwd=2)
legend(x="topright", legend=c("Cosine", "Metric Learning"), lwd=c(4,4), col=c("blue", "red"))
dev.off()


hccpclstats <- data.frame(ml=1-sapply(a$MLRanks, mean), cosine=1-sapply(a$cosranks, mean), size=sapply(a$MLRanks, length), name=names(a$MLRanks))

pdf(file.path("../../figs", "L1K_PCL_Scatter_HCC515.pdf"), width=10, height=8)
ggplot(hccpclstats[hccpclstats$size > 50,], aes(x=ml, y=cosine)) + geom_point(aes(size=log10(size)), color="blue") + 
  xlim(c(0.35,1)) + ylim(c(0.35,1)) + theme_minimal() + geom_abline(slope=1, intercept=0, linetype="dashed", color="black") + 
  geom_text_repel(aes(label=name), size=2) + ggtitle("Distribution of MoA ranks, HCC515") + xlab("Metric Learning") + ylab("Cosine")
dev.off()

pclstatssub <- hccpclstats[hccpclstats$size > 50,]
pdf("L1K_PCL_LeadingPCLsHCC515Table.pdf", width=8, height=8)
gridExtra::grid.table(pclstatssub[order(pclstatssub$ml - pclstatssub$cosine, decreasing=TRUE)[1:10],c(4,1,2,3)], rows=NULL)
dev.off()


### Cell Count

pdf(file.path("../../figs", "L1K_CellLines.pdf"), width=12, height=7)
ggplot(celldf, aes(x=sigcount.Freq, y=pertcount, color=factor(cell_lineage))) + geom_point(size=2) + theme_minimal() + 
  xlab("Signatures") + ylab("Unique Compounds") + geom_text_repel(aes(label=cell_id)) + ggtitle("Cell lines in L1000") + 
  scale_color_discrete(name = "Lineage")
dev.off()


# FDR of ranks:
auFDR <- data.frame()

for (acell in mycells){
  a <- readRDS(sprintf("%s_PCL_pcl_res.rds", acell))
  
  auFDR <- rbind(auFDR, data.frame(trainset="cell", testset=acell, 
                                   fdr=0.05, rate = mean(p.adjust(unlist(a$cosranks), method="fdr") < 0.05), metric="cosine"))
  auFDR <- rbind(auFDR, data.frame(trainset="cell", testset=acell, 
                                   fdr=0.05, rate = mean(p.adjust(unlist(a$MLRanks), method="fdr") < 0.05), metric="ml"))
  auFDR <- rbind(auFDR, data.frame(trainset="cell", testset=acell, 
                                   fdr=0.1, rate = mean(p.adjust(unlist(a$cosranks), method="fdr") < 0.1), metric="cosine"))
  auFDR <- rbind(auFDR, data.frame(trainset="cell", testset=acell, 
                                   fdr=0.1, rate = mean(p.adjust(unlist(a$MLRanks), method="fdr") < 0.1), metric="ml"))
  auFDR <- rbind(auFDR, data.frame(trainset="cell", testset=acell, 
                                   fdr=0.25, rate = mean(p.adjust(unlist(a$cosranks), method="fdr") < 0.25), metric="cosine"))
  auFDR <- rbind(auFDR, data.frame(trainset="cell", testset=acell, 
                                   fdr=0.25, rate = mean(p.adjust(unlist(a$MLRanks), method="fdr") < 0.25), metric="ml"))
}

pdf(file.path("../../figs", "L1K_PCL_FDR_10.pdf"), width=8, height=5)
ggplot(auFDR[auFDR$fdr == 0.1,], aes(x=testset, y=rate, fill=metric)) + geom_bar(position="dodge", stat="identity") + theme_minimal() + 
  ylim(c(0, 0.6)) + xlab("Cell Line") + ylab("Rate") + ggtitle("L1K PCL FDR=0.1")
dev.off()
pdf(file.path("../../figs", "L1K_PCL_FDR_25.pdf"), width=8, height=5)
ggplot(auFDR[auFDR$fdr == 0.25,], aes(x=testset, y=rate, fill=metric)) + geom_bar(position="dodge", stat="identity") + theme_minimal() + 
  ylim(c(0, 0.6)) + xlab("Cell Line") + ylab("Rate") + ggtitle("L1K PCL FDR=0.25")
dev.off()



# L1K Cell Line Specificity
CLS <- readRDS("L1KCellLineSpec_nLines=6_epoch=10_npcds=2000.rds")
CLSres <- CLS
CLSres$deltaAvgLoss <- 0
for (ii in seq((dim(CLS)[1])/3)){
  CLSres$deltaAvgLoss[(3*ii - 2):(3*ii)] <- CLSres$testAvgLoss[(3*ii - 2):(3*ii)] - CLSres$testAvgLoss[3*ii - 2]
}
CLSres <- CLSres[which(seq(dim(CLSres)[1]) %% 3 != 1),]

pdf(file.path("../figs", "CellLineSpec_Delta.pdf"), width=8, height=6)
ggplot(CLSres, aes(x=as.factor(trainSet), y=deltaAvgLoss, fill=cell_id)) + geom_boxplot() + theme_minimal() + 
  geom_hline(yintercept=0, linetype="dashed", color="black") + xlab("Training Set") + ylab("Test Loss - Cosine Test Loss")
dev.off()


# L1K cross validation
mcf7 <- readRDS("L1Kxval_epch=20_folds=5_cell=MCF7.rds")

pdf(file=file.path(outdir, "mcf7_xval_loss.pdf"), width=8, height=6)
plot(x=-20, y=-20, xlim=c(0, 20), ylim=c(-4,-2), xlab="Epoch", ylab="T loss", main="MCF7 5-Fold Cross validation Loss")
sapply(names(mcf7)[2:6], FUN=function(x) points(mcf7[[x]]$res$mean_train_losses, pch=16, col=alpha("blue",0.6)))
sapply(names(mcf7)[2:6], FUN=function(x) points(mcf7[[x]]$res$mean_valid_losses, pch=16, col=alpha("red",0.6)))
legend(x="topright", legend=c("Training loss", "Validation loss"), col=c("blue", "red"), lwd=c(4,4))
dev.off()

asc <- readRDS("L1Kntraining_epch=20_cell=ASC.rds")
plot(asc$trainClasses, asc$validLoss, pch=16, col=alpha("red", 0.6), xlab="Number of Training compounds", 
     ylab="T loss", main="Validation Loss")


# L1K cross validation - batch figures
jobdir <- './modelRuns'
figdir <- file.path(jobdir, "xvalFigs")

f <- list.files(jobdir, pattern="xval")

for (myf in f){
  myfname <- strsplit(myf, ".", fixed=TRUE)[[1]][1]
  mytitle <- strsplit(myfname, "cell=", fixed=TRUE)[[1]][2]
  makeXValPlots(file.path(jobdir, myf), fname=myfname, mytitle=mytitle, outpath=figdir)
}

# L1K ntraining - batch figures
jobdir <- './modelRuns'
figdir <- file.path(jobdir, "nTrainingFigs")

fnt <- list.files(jobdir, pattern="ntraining")

for (myf in fnt[2]){
  myfname <- strsplit(myf, ".", fixed=TRUE)[[1]][1]
  mytitle <- strsplit(myfname, "cell=", fixed=TRUE)[[1]][2]
  makeNTrainingPlots(file.path(jobdir, myf), fname=myfname, mytitle=mytitle, outpath = figdir)
}


####
# Cell Painting xval figures
####

setwd("~/Work/bhk/analysis/metric_learning/2022_cellpaint/models/")
dir.create("xvalFigs")
makeXValPlots("brayxval_epch=10_smp=20_folds=3.rds", fname="brayxval_epch=10_smp=20_folds", mytitle="Bray ds", outpath="xvalFigs")
makeXValPlots("brayxval_epch=2_folds=3.rds", fname="brayxval_epch=2_folds=3", mytitle="Bray ds", outpath="xvalFigs")



# Cell Painting xval replicates
topdir <- "~/Work/bhk/data/cellpainting/bray-2017"
brayds <- loadBrayData(file.path(topdir, "bray_2017_combined_profiles.rds"))
setwd("~/Work/bhk/analysis/metric_learning/2022_cellpaint/models/")

model <- torch_load("brayxval_epch=10_smp=20_folds=3_model1.pt")
brayres <- readRDS("brayxval_epch=10_smp=20_folds=3.rds")

mycps <- brayres$myfolds[[1]]

subds <- brayds$ds[brayds$metads$Metadata_pert_id %in% mycps | brayds$metads$Metadata_broad_sample_type == "control",]
submeta <- brayds$metads[brayds$metads$Metadata_pert_id %in% mycps | brayds$metads$Metadata_broad_sample_type == "control",]



######### Cell Painting Figures, Shantanu Meeting
cpdir <- "~/Work/bhk/analysis/metric_learning/2022_cellpaint/"
outdir <- file.path(cpdir, "bray")

brayds <- loadBrayData(file.path(topdir, "bray_2017_combined_profiles.rds"))
model1 <- torch_load(file.path(cpdir, "models", "brayxval_epch=5_smp=50_folds=3_model1.pt"))

pclds <- readRDS("~/Work/code/metriclearning/data/pclds.rds")


# Plot replicate similarities, ranks
ix <- sample(dim(brayds$ds)[1], 2000)

bgcos <- innerProduct("cosine", brayds$ds[ix,], brayds$ds[ix,])
bgcos <- bgcos[upper.tri(bgcos)]
bgsim <- innerProduct(model1, brayds$ds[ix,], brayds$ds[ix,])
bgsim <- bgsim[upper.tri(bgsim)]

cosperts <- list()
simperts <- list()

for (mypert in unique(brayds$metads$Metadata_pert_id[brayds$metads$Metadata_pert_type == "trt"])[1:5000]){
  if (length(cosperts) %% 100 == 0){print(length(cosperts))}
  jx <- which(brayds$metads$Metadata_pert_id == mypert)
  
  pertcos <- innerProduct("cosine", brayds$ds[jx, ], brayds$ds[jx,])
  pertsim <- innerProduct(model1, brayds$ds[jx, ], brayds$ds[jx,])
  cosperts[[mypert]] <- pertcos[upper.tri(pertcos)]
  simperts[[mypert]] <- pertsim[upper.tri(pertsim)]
}

pdf(file.path(outdir, "bray_cosine_replisim.pdf"), width=10, height=8)
plot(density(bgcos, bw=0.01), col="purple", lwd=2, xlim=c(-1,1), xlab="Cosine Similarity", ylab="Density", main="Bray Cosine Replicate similarity", ylim=c(0, 1.5))
lines(density(unlist(cosperts), bw=0.01), col="forestgreen", lwd=2)
legend(x="topleft", legend=c("All signatures", "Replicates"), col=c("purple", "forestgreen"), lwd=c(4,4))
dev.off()

pdf(file.path(outdir, "bray_learnedMetric_replisim_epch=5_smp=50_folds=3_mdl1.pdf"), width=10, height=8)
plot(density(bgsim, bw=0.01), col="blue", lwd=2, xlim=c(-1,1), xlab="Learned Metric", ylab="Density", main="Bray Learned Metric replicate similarity")
lines(density(unlist(simperts), bw=0.01), col="red", lwd=2)
legend(x="topleft", legend=c("All signatures", "Replicates"), col=c("blue", "red"), lwd=c(4,4))
dev.off()

rankcos <- rankVectors(unlist(cosperts), bgcos)
ranksim <- rankVectors(unlist(simperts), bgsim)

pdf(file.path(outdir, "bray_ReplicateRanks.pdf"), width=10, height=8)
plot(ecdf(ranksim), col="red", lwd=3, xlab=c("Replicate Rank"), ylab="Cumulative Density", xlim=c(0,1), main="Replicate Ranks on Bray dataset")
lines(ecdf(rankcos), col="blue", lwd=3)
legend(x="bottomright", legend=c(sprintf("Learned Metric = %0.3f", 1-mean(ranksim)), sprintf("Cosine = %0.3f", 1-mean(rankcos))), col=c("red", "blue"), lwd=c(4,4))
dev.off()



# PCLs

pclcos <- list()
pclsim <- list()

for (ii in seq(dim(pclds$pcldf)[1])){
  mypcl <- pclds$pcldf$pclid[ii]
  pclperts <- pclds$pertids[[ii]]
  
  brayperts <- intersect(brayds$metads$Metadata_pert_id, pclperts)
  
  if (length(brayperts > 1)){
    print(sprintf("%d: %s", ii, mypcl))
    jx <- which(brayds$metads$Metadata_pert_id %in% brayperts)
    
    filt <- outer(brayds$metads$Metadata_pert_id[jx], brayds$metads$Metadata_pert_id[jx], '!=')
    filt <- filt[upper.tri(filt)]
    
    pcos <- innerProduct("cosine", brayds$ds[jx, ], brayds$ds[jx,])
    psim <- innerProduct(model1, brayds$ds[jx, ], brayds$ds[jx,])
    
    pclcos[[mypcl]] <- (pcos[upper.tri(pcos)])[filt]
    pclsim[[mypcl]] <- (psim[upper.tri(psim)])[filt]
  }
}

pdf(file.path(outdir, "bray_cosine_pclsim.pdf"), width=10, height=8)
plot(density(unlist(pclcos), bw=0.01), col="forestgreen", lwd=2, xlim=c(-1,1), xlab="Cosine Similarity", ylab="Density", main="Bray Cosine PCL Similarity", ylim=c(0,1.5))
lines(density(bgcos, bw=0.01), col="purple", lwd=2)
legend(x="topleft", legend=c("All signatures", "Same PCL"), col=c("purple", "forestgreen"), lwd=c(4,4))
dev.off()

pdf(file.path(outdir, "bray_learnedMetric_pclsim_epch=5_smp=50_folds=3_mdl1.pdf"), width=10, height=8)
plot(density(unlist(pclsim), bw=0.01), col="red", lwd=2, xlim=c(-1,1), xlab="Learned Metric", ylab="Density", main="Bray Learned Metric PCL similarity", ylim=c(0,4))
lines(density(bgsim, bw=0.01), col="blue", lwd=2)
legend(x="topleft", legend=c("All signatures", "Same PCL"), col=c("blue", "red"), lwd=c(4,4))
dev.off()

rankpclcos <- rankVectors(unlist(pclcos), bgcos)
rankpclsim <- rankVectors(unlist(pclsim), bgsim)

pdf(file.path(outdir, "bray_PCLRanks.pdf"), width=10, height=8)
plot(ecdf(rankpclsim), col="red", lwd=3, xlab=c("Rank"), ylab="Cumulative Density", xlim=c(0,1), main="Same PCL Ranks on Bray dataset")
lines(ecdf(rankpclcos), col="blue", lwd=3)
legend(x="bottomright", legend=c(sprintf("Learned Metric = %0.3f", 1-mean(rankpclsim)), sprintf("Cosine = %0.3f", 1-mean(rankpclcos))), col=c("red", "blue"), lwd=c(4,4))
dev.off()




# Eigenvalue distribution 
mycell <- "HEPG2"

ds <- parse_gctx(get_level5_ds(datapath), cid=siginfo$sig_id[siginfo$cell_id == mycell & siginfo$pert_type == "trt_cp"], rid = landmarks$pr_gene_id)

mymodel <- torch_load(file.path(l1kdir, "models", "L1Kmetric_epch=30_cell=HEPG2_model.pt"))
modelmat <- mymodel(torch_tensor(t(ds@mat), dtype=torch_float()))

pcHEPG2 <- prcomp(t(ds@mat), center=TRUE)
pcMLHEPG2 <- prcomp(as.matrix(modelmat), center=TRUE)

pctvar <- pcHEPG2$sdev^2/sum(pcHEPG2$sdev^2)
pctvarML <- pcMLHEPG2$sdev^2/sum(pcMLHEPG2$sdev^2)

pdf(file.path(outdir, "fig6_VarianceExplained_HEPG2.pdf"), width=8, height=6)
ggplot(data.frame(index=seq(978), native=pctvar, embedded=pctvarML), aes(x=native, y=embedded)) + geom_point() + scale_x_log10() + scale_y_log10() + 
  geom_abline(intercept=0, slope=1, col="blue", lty=2) + theme_minimal() + coord_cartesian(xlim=c(1e-4, 3e-1), ylim=c(1e-4, 3e-1)) + xlab("Native Pct Variance") +
  ylab("Embedded Pct Variance") + ggtitle("Variance Explained by Metric Learning in HEPG2")
dev.off()

pdf(file.path(outdir, "fig6_EigenvalueRatio_HEPG2.pdf"), width=8, height=6)
ggplot(data.frame(index=seq(978), native=pctvar, embedded=pctvarML, ratio=pctvarML/pctvar), aes(x=index, y=ratio)) + geom_point() + theme_minimal() + 
  xlab("Eigenvalue index") + ylab("Ratio of eigenvalues of embedding to native") + ggtitle("Ratio of eigenvalues of HEPG2 Embedding vs native")
dev.off()