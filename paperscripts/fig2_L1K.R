source("figinit.R")
source("paperFigFuncs.R")
source("../learningFunctions.R")



#### Figure 2: L1000 performance ####
#### Compute L1000 data: ####
# Fig 2a: Pick a representative cell line, say A375; compute replicate distribution

if (!file.exists(file.path(outdir, "../figspaper_res", "A375_SimDS.rds"))){
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
  
  saveRDS(list(mlsame=mlsame, mldiff=mldiff, cossame=cossame, cosdiff=cosdiff), file=file.path(outdir, "../figspaper_res", "A375_SimDS.rds"))
} else {
  print("reading A375_SimDS.rds")
  a375simds <- readRDS(file.path(outdir, "../figspaper_res", "A375_SimDS.rds"))
  attach(a375simds)
}


#### Fig 2a - cosine similarity for replicates vs all ####
# Pick a representative cell line, say A375

pdf(file.path(outdir, sprintf("fig2a_%s_pdf.pdf", mycell)), width=8, height=6)
plot(density(unlist(cosdiff), bw=0.01), col="forestgreen", lwd=1.5, lty=5, xlim=c(-0.5,1), 
     xlab="Similarity", ylab="Density", main=sprintf("%s Balanced Similarity", mycell))
lines(density(unlist(cossame), bw=0.01), col="orange", lwd=1.5, lty=5)
lines(density(unlist(mldiff), bw=0.01), col="blue", lwd=1.5)
lines(density(unlist(mlsame), bw=0.01), col="firebrick3", lwd=1.5)
lines(c(0, 0), c(0, 100), col="grey", lty=2)
legend(x="topright", legend=c("Metric learning replicates", "Metric Learning all pairs", 
                              "Cosine replicates", "Cosine all pairs"), lwd=2.5, col=c("firebrick3", "blue", "orange", "forestgreen"), 
       lty=c(1,1,5,5))
dev.off()

pdf(file.path(outdir, sprintf("fig2a_%s_cosine_pdf.pdf", mycell)), width=8, height=6)
plot(density(unlist(cosdiff), bw=0.01), col="forestgreen", lwd=2, lty=1, xlim=c(-0.5,1), 
     xlab="Similarity", ylab="Density", main=sprintf("%s Balanced Similarity", mycell))
lines(density(unlist(cossame), bw=0.01), col="orange", lwd=2, lty=1)
lines(c(0, 0), c(0, 100), col="grey", lty=2)
legend(x="topright", legend=c("Cosine replicates", "Cosine all pairs"), lwd=2.5, col=c("orange", "forestgreen"), 
       lty=c(1,1))
dev.off()



#### Fig 2b - one example of cross validated replicate recall curves ####

pdf(file.path(outdir, sprintf("fig2b_%s_xvalBalRankCDF.pdf", mycell)), width=8, height=6)
plot(c(-10), c(-10), xlim=c(0,1), ylim=c(0,1), xlab="Balanced Rank", ylab="Cumulative Fraction", 
     main=sprintf("%s Balanced Rank", mycell))
for (ii in seq(5)){
  lines(ecdf(unlist(balancedSample(xvalreps$cosRanks[[ii]], 100))), col="orange", lty=5, lwd=1.5)
}
for (ii in seq(5)){
  lines(ecdf(unlist(balancedSample(xvalreps$mlRanks[[ii]], 100))), col="firebrick3", lwd=1.5)
}
legend(x="bottomright", legend=c("Metric Learning", "Cosine"), col=c("firebrick3", "orange"), lwd=2.5,
       lty=c(1,5), bg="white")
dev.off()



#### Fig 2c - summary cross validated replicate recall (bar plots, auRank) ####

# Get cell names
mycells <- sapply(strsplit(sapply(strsplit(list.files(l1kdir, pattern="L1Kxval"), "cell="), 
                                  FUN=function(x) x[[2]]), ".rds"), FUN=function(x) x[[1]])

# !!!!!! Fix this once you have the all_xval models
mycells <- setdiff(mycells, c("all", "MCF10A", "SKL"))

if (!file.exists(file.path(outdir, "../figspaper_res", "L1K_xvalres.rds"))){
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
} else {
  xvalres <- readRDS(file.path(outdir, "../figspaper_res/L1K_xvalres.rds"))
}

#### mean balanced AUC ####
# Consider also plotting all 100 balanced AUCs rather than taking the mean for each fold
# The results are stored as ranks on [0,1], where 0 is the best rank. 
# To convert to AUC, take 1 - mean balanced rank.

if (!file.exists(file.path(outdir, "../figspaper_res/L1K_PertRes.rds"))){
  L1KPertBalAUCML <- 1 - sapply(mycells, FUN=function(x) 
    sapply(xvalres[[x]]$mlRanks, FUN=function(y) 
      mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(y)))))))
  
  L1KPertBalAUCCos <- 1 - sapply(mycells, FUN=function(x) 
    sapply(xvalres[[x]]$cosRanks, FUN=function(y) 
      mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(y)))))))
  
  ### mean compound AUC - avoids sampling
  L1KPertMeanAUCML <- 1 - sapply(mycells, FUN=function(x) 
    mean(sapply(xvalres[[x]]$mlRanks, FUN=function(y)
      mean(sapply(y, FUN=mean)))))
  
  L1KPertMeanAUCCos <- 1 - sapply(mycells, FUN=function(x) 
    mean(sapply(xvalres[[x]]$cosRanks, FUN=function(y)
      mean(sapply(y, FUN=mean)))))
  
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
  
  saveRDS(list(L1KPertBalAUCML=L1KPertBalAUCML, 
               L1KPertBalAUCCos=L1KPertBalAUCCos, 
               L1KPertMeanAUCML=L1KPertMeanAUCML,
               L1KPertMeanAUCCos=L1KPertMeanAUCCos,
               L1KPertBalFDRML=L1KPertBalFDRML, 
               L1KPertBalFDRCos=L1KPertBalFDRCos), 
          file=file.path(outdir, "../figspaper_res/L1K_PertRes.rds"))
} else {
  pertRes <- readRDS(file.path(outdir, "../figspaper_res/L1K_PertRes.rds"))
  attach(pertRes)
}

#ggplot(df, aes(x=Var2, y=value, fill=method)) + geom_boxplot() + facet_wrap(~Var2, scale="free")
pdf(file=file.path(outdir, "fig2c_L1KMeanBalAUC.pdf"), width=8, height=6)
ggplot(rbind(data.frame(melt(L1KPertBalAUCCos), method="Cosine"), data.frame(melt(L1KPertBalAUCML), method="Metric Learning")), aes(x=Var2, y=value, fill=method)) + 
  geom_boxplot() + geom_point(position=position_jitterdodge()) + 
  theme_minimal() + ylim(c(0.5, 1)) + geom_hline(yintercept=0.5, col="grey", linetype=2) + 
  xlab("Cell Line") + ylab("Mean Balanced AUC") + ggtitle("L1K Cross-validated Mean Balanced AUC") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
dev.off()


pdf(file=file.path(outdir, "fig2c_L1KMeanPertAUC.pdf"), width=8, height=6)
ggplot(data.frame(ML=L1KPertMeanAUCML, Cos=L1KPertMeanAUCCos, cell=names(L1KPertMeanAUCML)), 
       aes(x=Cos, y=ML, color=cell)) + geom_point() + xlim(c(0.6, 1)) + theme_minimal() + 
  ylim(c(0.6,1)) + geom_abline(intercept=0, slope=1, linetype=2, col="grey") + xlab("Cosine Mean AUC") +
  ylab("Metric Learning Mean AUC") + ggtitle("L1K Cross Validated Mean Compound AUC") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
dev.off()


pdf(file=file.path(outdir, "fig2d_L1Kpertfdr05.pdf"), width=8, height=6)
ggplot(rbind(L1KPertBalFDRCos, L1KPertBalFDRML), aes(x=cell, y=fdr05, fill=method)) + 
  geom_bar(stat="identity", position="dodge") + theme_minimal() + ylab("Balanced FDR < 0.05") + 
  ggtitle("L1K Replicate Balanced FDR < 0.05") + ylim(c(0, 0.6)) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
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

#### Fig 2e - summary PCL recall (auRank?) ####
pclds <- readRDS("data/pclds.rds")

pclcells <- sapply(strsplit(sapply(strsplit(list.files(l1kdir, pattern="L1Kmetric"), "cell="), 
                                   FUN=function(x) x[[2]]), ".rds"), FUN=function(x) x[[1]])
pclcells <- setdiff(pclcells, "all")

if (!file.exists(file.path(outdir, "../figspaper_res/L1K_PCLRes.rds"))){
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
  
  # SNR 
  L1KMoASNRML <- sapply(pclres, FUN=function(x) (sapply(x$mlPCLs$setSims, mean) - mean(x$mlPCLs$allSims))/sd(x$mlPCLs$allSims))
  L1KMoASNRCos <- sapply(pclres, FUN=function(x) (sapply(x$cosPCLs$setSims, mean) - mean(x$cosPCLs$allSims))/sd(x$cosPCLs$allSims))
  L1KMoASize <- sapply(pclres, FUN=function(x) sapply(x$mlRanks, length))
  
  # MoA FDR
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
  
  
  saveRDS(list(pclres=pclres, 
               L1KMOABalAUCML=L1KMOABalAUCML,
               L1KMOABalAUCCos=L1KMOABalAUCCos, 
               L1KMoASNRML=L1KMoASNRML, 
               L1KMoASNRCos=L1KMoASNRCos,
               L1KMoABalFDRML=L1KMoABalFDRML,
               L1KMoABalFDRCos=L1KMoABalFDRCos), 
          file=file.path(outdir, "../figspaper_res/L1K_PCLRes.rds"))
} else {
  L1K_PCLRes <- readRDS(file.path(outdir, "../figspaper_res/L1K_PCLRes.rds"))
  attach(L1K_PCLRes)
}


# Alt metric PCLres
# Generate with l1kAltMetrics.R; add arguments below. 

# runGenerateAltMetrics()
L1K_PCLResAlt <- readRDS(file.path(outdir, "../figspaper_res/L1K_PCLResAltMetrics.rds"))
attach(L1K_PCLResAlt)

# Fix the error ranking wtcs guys because of ties
for (mycell in names(pclresAlt)){
  pclresAlt[[mycell]]$wtcsRanks <- listify(pclresAlt[[mycell]]$wtcsPCLs$setSims, rankVectors(unlist(pclresAlt[[mycell]]$wtcsPCLs$setSims), pclresAlt[[mycell]]$wtcsPCLs$allSims))
}

altMoAfdr <- data.frame(rbind(cbind(fdr=0.01, sapply(pclresAlt, FUN=function(mycell) sapply(mycell[5:8], FUN=function(x) mean(p.adjust(unlist(balancedSample(x, k=1000)), method="fdr") < 0.01)))),
  cbind(fdr=0.05, sapply(pclresAlt, FUN=function(mycell) sapply(mycell[5:8], FUN=function(x) mean(p.adjust(unlist(balancedSample(x, k=1000)), method="fdr") < 0.05)))),
  cbind(fdr=0.10, sapply(pclresAlt, FUN=function(mycell) sapply(mycell[5:8], FUN=function(x) mean(p.adjust(unlist(balancedSample(x, k=1000)), method="fdr") < 0.10))))))
                   

#### Example Alt metric comparison - FDR CDF ####
pdf(file.path(outdir, "Sfig_L1KMoAAltMetricFDRcomparison.pdf"), width=8, height=6)
for (mycell in names(pclresAlt)) {
  print(mycell)
  plot(ecdf(sort(p.adjust(unlist(balancedSample(pclresAlt[[mycell]]$cosRanks, k=1000)), method="fdr"))), col="#F8766D", ylab="Fraction", 
       xlab="FDR", main=sprintf("%s Metric FDR CDF comparison", mycell), xlim=c(0,1), lwd=2)
  lines(ecdf(sort(p.adjust(unlist(balancedSample(pclresAlt[[mycell]]$wtcsRanks, k=1000)), method="fdr"))), col=cbbPalette[7], lwd=2)
  lines(ecdf(sort(p.adjust(unlist(balancedSample(pclresAlt[[mycell]]$pearsonRanks, k=1000)), method="fdr"))), col=cbbPalette[2], lwd=2)
  lines(ecdf(sort(p.adjust(unlist(balancedSample(pclresAlt[[mycell]]$spearRanks, k=1000)), method="fdr"))), col=cbbPalette[4], lwd=2)
  lines(ecdf(sort(p.adjust(unlist(balancedSample(pclres[[mycell]]$mlRanks, k=1000)), method="fdr"))), col="#00BFC4", lwd=2)
  legend(x="bottomright", legend=c("Metric Learning", "Cosine", "Pearson", "Spearman", "WTCS"), 
         col=c("#00BFC4", "#F8766D", cbbPalette[c(2,4,7)]), lwd=5)
}
dev.off()


combpcl <- rbind(L1KMoABalFDRML[,c(2, 4, 5)], L1KMoABalFDRCos[,c(2,4,5)])
combpcl <- rbind(combpcl, data.frame(fdr05=as.numeric(altMoAfdr[6,2:13]), cell=colnames(altMoAfdr)[2:13], method="wtcs"))
combpcl <- rbind(combpcl, data.frame(fdr05=as.numeric(altMoAfdr[7,2:13]), cell=colnames(altMoAfdr)[2:13], method="pearson"))
combpcl <- rbind(combpcl, data.frame(fdr05=as.numeric(altMoAfdr[8,2:13]), cell=colnames(altMoAfdr)[2:13], method="spearman"))

# Fix colors
pdf(file.path(outdir, "Sfig_L1KMoAAltMetricRecallComparison.pdf"), width=8, height=6)
ggplot(combpcl[!(combpcl$cell %in% c("HUVEC", "JURKAT", "NPC", "BT20")),], aes(x=cell, y=fdr05, fill=method)) + 
  geom_bar(stat="identity", position="dodge") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + xlab("Cell Line") + ylab("MoA Recall at FDR = 0.05") + 
  ggtitle("Comparison of metrics on L1000: MoA Recall at FDR < 0.05") + scale_fill_manual(values=c("#F8766D", "#00BFC4", cbbPalette[c(2,4,7)]))
dev.off()

#### Metric Learning MoA analysis ####
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


# PCL SNR  
# This is mostly formatting
pclsnrML <- sapply(seq(92), FUN=function(y) sapply(L1KMoASNRML, FUN=function(x) x[y]))
pclsnrCos <- sapply(seq(92), FUN=function(y) sapply(L1KMoASNRCos, FUN=function(x) x[y]))
pclsnrSize <- sapply(seq(92), FUN=function(y) sapply(L1KMoASize, FUN=function(x) x[y]))

pclsnrDF <- data.frame(pcl=pclds$pcldf$pclname, 
                       cellid=as.character(sapply(rownames(pclsnrML), FUN=function(x) rep(x,92))),
                       size=as.numeric(pclsnrSize),
                       ML=as.numeric(pclsnrML), Cos=as.numeric(pclsnrCos))

pdf(file.path(outdir, "fig4_L1KSNR_MoA_scatter.pdf"), width=8, height=6)
ggplot(pclsnrDF[pclsnrDF$size > 50,], aes(x=Cos, y=ML, color=cellid)) + geom_point(aes(size=log10(size)), alpha=0.8) + 
  theme_minimal() + geom_abline(slope=1, intercept=0, linetype="dashed", color="black") + 
  geom_text_repel(aes(label=pcl, size=2)) + xlab("Cosine SNR") + ylab("Metric learning SNR") + 
  ggtitle("L1K MoA SNR")
dev.off()
  

ggplot(hccpclstats[hccpclstats$size > 50,], aes(x=ml, y=cosine)) + geom_point(aes(size=log10(size)), color="blue") + 
  xlim(c(0.35,1)) + ylim(c(0.35,1)) + theme_minimal() + geom_abline(slope=1, intercept=0, linetype="dashed", color="black") + 
  geom_text_repel(aes(label=name), size=2) + ggtitle("Distribution of MoA ranks, HCC515") + xlab("Metric Learning") + ylab("Cosine")



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
       aes(x=ml-cos, color=cell, fill=cell)) + geom_density(alpha=0.6) + theme_minimal() + #+ geom_density(alpha=0.6, position="stack")
  geom_vline(xintercept=0, linetype=2, col="black") + xlab("MoA SNR: Metric Learning - Cosine") +
  ggtitle("Change in MoA Mean SNR with metric learning vs cosine")
dev.off()


# Fig 2f - FDR plots for replicates, PCLs
pdf(file=file.path(outdir, "fig2f_L1KMoAfdr10.pdf"), width=8, height=6)
ggplot(rbind(L1KMoABalFDRCos, L1KMoABalFDRML), aes(x=cell, y=fdr10, fill=method)) + 
  geom_bar(stat="identity", position="dodge") + theme_minimal() + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  ylab("Balanced FDR < 0.1") + ggtitle("L1K MoA Balanced FDR < 0.1") + ylim(c(0, 0.5)) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
dev.off()

pdf(file=file.path(outdir, "fig2f_L1KMoAfdr05.pdf"), width=8, height=6)
ggplot(rbind(L1KMoABalFDRCos, L1KMoABalFDRML), aes(x=cell, y=fdr05, fill=method)) + 
  geom_bar(stat="identity", position="dodge") + theme_minimal() + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  ylab("Balanced FDR < 0.05") + ggtitle("L1K MoA Balanced FDR < 0.05") + ylim(c(0, 0.5)) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
dev.off()






##### Supp Figs - L1000 #####

cpCounts <- sapply(pclcells, FUN=function(x) length(unique(siginfo$pert_iname[siginfo$cell_id == x & siginfo$pert_type == "trt_cp"])))
sigCounts <- sapply(pclcells, FUN=function(x) sum(siginfo$cell_id == x & siginfo$pert_type == "trt_cp"))

l1kCounts <- data.frame(cellLine=pclcells, compounds=cpCounts, signatures=sigCounts)

l1kCounts <- rbind(l1kCounts, data.frame(cellLine="all", 
                                         compounds=length(unique(siginfo$pert_iname[siginfo$pert_type == "trt_cp"])), 
                                         signatures=sum(siginfo$pert_type == "trt_cp")))


l1kCountsAll <- data.frame(cellLine = unique(siginfo$cell_id[siginfo$pert_type == "trt_cp"]), 
                           compounds = sapply(unique(siginfo$cell_id[siginfo$pert_type == "trt_cp"]), 
                                              FUN=function(x) length(unique(siginfo$pert_iname[siginfo$cell_id == x & siginfo$pert_type == "trt_cp"]))),
                           signatures = sapply(unique(siginfo$cell_id[siginfo$pert_type == "trt_cp"]), 
                                               FUN=function(x) sum(siginfo$cell_id == x & siginfo$pert_type == "trt_cp")))

l1kCountsTop <- l1kCountsAll[l1kCountsAll$signatures > 2500, ]
l1kCountsTop <- rbind(l1kCountsTop, data.frame(cellLine="other",
                                               compounds=0, 
                                               signatures=sum(siginfo$pert_type == "trt_cp") - sum(l1kCountsTop$signatures)))
l1kCountsTop$cellLine <- factor(l1kCountsTop$cellLine, levels=c(l1kCountsTop$cellLine[-length(l1kCountsTop$cellLine)], "other"))

pdf(file=file.path(outdir, "Sfig_L1KCompoundCountsByCell.pdf"), width=8, height=6)
ggplot(l1kCounts, aes(x=compounds, y=signatures, color=cellLine)) + geom_point(size=2) + theme_minimal() + 
  scale_x_continuous(trans="log10", limits=c(100, 35000), breaks=c(100, 300, 1000, 3000, 10000, 30000)) + scale_y_continuous(trans="log10", limits=c(1000, 1e6)) + 
  xlab("Unique Compounds") + ylab("Compound Signatures") + ggtitle("Number of Signatures and Unique Compounds in L1000 dataset (2020) by Cell Line") 
dev.off()

pdf(file=file.path(outdir, "Sfig_L1KCompoundCountsByCellAll.pdf"), width=8, height=6)
ggplot(l1kCounts, aes(x=compounds, y=signatures, color=cellLine)) + geom_point(data=l1kCountsAll, aes(x=compounds, y=signatures), size=1, color="grey") + 
  geom_point(size=3) + theme_minimal() + 
  scale_x_continuous(trans="log10", limits=c(10, 35000), breaks=c(10, 30, 100, 300, 1000, 3000, 10000, 30000)) + scale_y_continuous(trans="log10", limits=c(10, 1e6)) + 
  xlab("Unique Compounds") + ylab("Compound Signatures") + ggtitle("Number of Signatures and Unique Compounds in L1000 dataset (2020) by Cell Line")
dev.off()

pdf(file=file.path(outdir, "Sfig_L1KCellLineBarChart.pdf"), width=6, height=6)
ggplot(l1kCountsTop, aes(fill=cellLine, y=signatures, x=1)) + geom_bar(stat="identity") + scale_fill_manual(values=rep(cbPalette, 4)) +
  geom_text(aes())
dev.off()





