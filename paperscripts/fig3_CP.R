source("figinit.R")
source("../experimentFunctions.R")

#### Figure 3: Cell Painting performance ####

# Fig 3a - cosine similarity example
# Fig 3b - cross validated replicate recall curves
# Fig 3c - summary cross validated replicate recall 
# Fig 3d - PCL recall example
# Fig 3e - summary PCL recall (auRank)
# Fig 3f - FDR plots for replicates, PCLs

braydir <- file.path(cpdir, "bray/models")
pclds <- readRDS("data/pclds.rds")


# Load or generate results by running models:
if (!file.exists(file.path(outdir, "../figspaper_res/cpSummaryRes.rds"))){
  brayds <- loadBrayData(braypath)
  lincsds1 <- loadLincsData(file.path(lincspath, "2016_04_01_a549_48hr_batch1_dmso_spherized_profiles_with_input_normalized_by_dmso.rds"), byCell = 1, splitGrps = 0)
  lincsds2 <- loadLincsData(file.path(lincspath, "2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_dmso.rds"), byCell = 0, splitGrps = 0)
  
  Lds2_a549 <- splitCPds(lincsds2, which(lincsds2$metads$Metadata_cell_id == "A549"))
  Lds2_mcf7 <- splitCPds(lincsds2, which(lincsds2$metads$Metadata_cell_id == "MCF7"))
  Lds2_u2os <- splitCPds(lincsds2, which(lincsds2$metads$Metadata_cell_id == "U2OS"))
  
  lincsds2_hidose <- c()
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
  
  
  #### Compound Replicates
  # Take care to exclude those compounds that are not compound treatments. Note that the LINCS datasets still have substantial class imbalances.
  brayReps <- getCPReps(brayds$ds, mymodel=braymodel, pertlabels=brayds$metads$Metadata_pert_id, datalabel="Bray DS centered full")
  lincs1Reps <- getCPReps(lincsds1$ds, mymodel=lincs1model, pertlabels=lincsds1$metads$Metadata_pert_id, datalabel="LINCS A549 DMSO full model")
  lincs2Reps <- getCPReps(lincsds2$ds, mymodel=lincs2model, pertlabels=lincsds2$metads$Metadata_pert_id, datalabel="LINCS Batch2 DMSO full model")
  
  lincs2Reps_A549 <- getCPReps(Lds2_a549$ds, mymodel=lincs2model, pertlabels=Lds2_a549$metads$Metadata_pert_id, datalabel="LINCS Batch2 DMSO full model, A549")
  lincs2Reps_MCF7 <- getCPReps(Lds2_mcf7$ds, mymodel=lincs2model, pertlabels=Lds2_mcf7$metads$Metadata_pert_id, datalabel="LINCS Batch2 DMSO full model, MCF7")
  lincs2Reps_U2OS <- getCPReps(Lds2_u2os$ds, mymodel=lincs2model, pertlabels=Lds2_u2os$metads$Metadata_pert_id, datalabel="LINCS Batch2 DMSO full model, U2OS")
  
  
  #### Cross validation results
  brayXvalSets <- list(file.path(braydir, list.files(braydir, pattern="brayxval_epch=10_smp=20")), 
                       file.path(braydir, list.files(braydir, pattern="brayxval_epch=5_smp=50")))
  
  brayXvalLabs <- c("brayxval_10x20", "brayxval_5x50") #, ", "lincs2_10x5")
  
  xvalres <- list()
  for (ii in seq_along(brayXvalSets)){
    myxvalset <- brayXvalSets[[ii]]
    modlab <- brayXvalLabs[ii]
    xvalds <- readRDS(myxvalset[grep("rds", myxvalset)])
    xvalres[[ii]] <- getCPXvalReps(brayds$ds, brayds$metads, xvalds$myfolds, myxvalset[-length(myxvalset)], modlab)
  }
  
  lincs1XvalSets <- list(file.path(cpdir, "lincs/models", list.files(file.path(cpdir, "lincs/models"), pattern="lincsxval_a549_2016_dmso_epch=10_folds=5")))
  lincs1XvalLabs <- c("lincs1_10x5")
  
  for (ii in seq_along(lincs1XvalSets)){
    myxvalset <- lincs1XvalSets[[ii]]
    modlab <- lincs1XvalLabs[ii]
    xvalds <- readRDS(myxvalset[grep("rds", myxvalset)])
    xvalres[[length(xvalres)+1]] <- getCPXvalReps(lincsds1$ds, lincsds1$metads, xvalds$myfolds, myxvalset[-length(myxvalset)], modlab)
  }
  
  # Add LINCS2
  
  
  #### MoA Similarities
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

  saveRDS(file = file.path(outdir, "../figspaper_res/cpSummaryRes.rds"), 
          object=list(lincsMoADS=lincsMoADS, 
                      brayReps=brayReps,
                      lincs1Reps=lincs1Reps,
                      lincs2Reps=lincs2Reps,
                      lincs2Reps_A549=lincs2Reps_A549,
                      lincs2Reps_MCF7=lincs2Reps_MCF7,
                      lincs2Reps_U2OS=lincs2Reps_U2OS,
                      xvalres=xvalres, 
                      brayMoA=brayMoA,
                      lincs1MoA=lincs1MoA,
                      lincs2MoA=lincs2MoA,
                      lincs2MoAByCell=lincs2MoAByCell, 
                      lincs2MoAByCellHighDose=lincs2MoAByCellHighDose,
                      brayMoANative=brayMoANative,
                      lincs1MoANative=lincs1MoANative,
                      lincs2MoANative=lincs2MoANative))
} else {
  cpds <- readRDS(file.path(outdir, "../figspaper_res/cpSummaryRes.rds"))
  attach(cpds)
}



##### Results tables

# xval summary
# use bray5x50
brayAUCML <- sapply(xvalres[[2]]$mlRanks, FUN=function(x) 1 - mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(x))))))
brayAUCCos <- sapply(xvalres[[2]]$cosRanks, FUN=function(x) 1 - mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(x))))))

brayFDRML <- data.frame(fdr01=sapply(xvalres[[2]]$mlRanks, FUN=function(x) mean(sapply(seq(10), FUN=function(z) mean(p.adjust(unlist(balancedSample(x)), "fdr") < 0.01)))),
                        fdr05=sapply(xvalres[[2]]$mlRanks, FUN=function(x) mean(sapply(seq(10), FUN=function(z) mean(p.adjust(unlist(balancedSample(x)), "fdr") < 0.05)))),
                        fdr10=sapply(xvalres[[2]]$mlRanks, FUN=function(x) mean(sapply(seq(10), FUN=function(z) mean(p.adjust(unlist(balancedSample(x)), "fdr") < 0.10)))))
brayFDRCos <- data.frame(fdr01=sapply(xvalres[[2]]$cosRanks, FUN=function(x) mean(sapply(seq(10), FUN=function(z) mean(p.adjust(unlist(balancedSample(x)), "fdr") < 0.01)))),
                        fdr05=sapply(xvalres[[2]]$cosRanks, FUN=function(x) mean(sapply(seq(10), FUN=function(z) mean(p.adjust(unlist(balancedSample(x)), "fdr") < 0.05)))),
                        fdr10=sapply(xvalres[[2]]$cosRanks, FUN=function(x) mean(sapply(seq(10), FUN=function(z) mean(p.adjust(unlist(balancedSample(x)), "fdr") < 0.10)))))



##### Mechanism of Action Summary #####
moaCPdf <- data.frame()
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(brayMoA, dsname="CDRP MoA", balSampT = 100))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs1MoA, dsname="Lincs1 MoA", balSampT = 100))

moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoAByCell$A549, dsname="Lincs2 MoA A549"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoAByCell$MCF7, dsname="Lincs2 MoA MCF7"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoAByCell$U2OS, dsname="Lincs2 MoA U2OS"))

moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs1MoANative, dsname="Lincs1 Native MoA"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoABCNative$A549, dsname="Lincs2 MoA Native A549"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoABCNative$MCF7, dsname="Lincs2 MoA Native MCF7"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoABCNative$U2OS, dsname="Lincs2 MoA Native U2OS"))

moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoAByCellHighDose$A549, dsname="Lincs2 MoA High Dose A549"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoAByCellHighDose$MCF7, dsname="Lincs2 MoA High Dose MCF7"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoAByCellHighDose$U2OS, dsname="Lincs2 MoA High Dose U2OS"))



##### Fig 3a - Cell Painting Cosine

pdf(file=file.path(outdir, "fig3a_bray_simpdfs.pdf"), width=8, height=6)
plot(density(brayReps$inProdReps$diff, bw=0.01), col="blue", lwd=1.5, xlim=c(-1, 1), 
     xlab="Similarity", ylab="Density", main=sprintf("CDRP similarity density"))
lines(density(unlist(balancedSample(brayReps$inProdReps$same, 100)), bw=0.01), col="firebrick3", lwd=1.5)
lines(density(brayReps$cosReps$diff, bw=0.01), col="forestgreen", lwd=1.5, lty=5)
lines(density(unlist(balancedSample(brayReps$cosReps$same, 100)), bw=0.01), col="orange", lwd=1.5, lty=5)
lines(c(0, 0), c(0, 100), col="grey", lty=2)
legend(x="topright", legend=c("Metric learning replicates", "Metric Learning all pairs", 
                              "Cosine replicates", "Cosine all pairs"), lwd=2.5, col=c("firebrick3", "blue", "orange", "forestgreen"), 
       lty=c(1,1,5,5))
dev.off()


L1Perts <- intersect(names(lincs1Reps$inProdReps$same), unique(lincsds1$metads$Metadata_pert_id[lincsds1$metads$Metadata_pert_type == "trt"]))
L2Perts <- intersect(names(lincs2Reps$inProdReps$same), unique(lincsds2$metads$Metadata_pert_id[lincsds2$metads$Metadata_pert_type == "trt"]))

pdf(file=file.path(outdir, "fig3a2_lincs1_simpdfs.pdf"), width=8, height=6)
plot(density(lincs1Reps$cosReps$diff, bw=0.01), col="forestgreen", lwd=1.5, lty=5, xlab="Similarity", ylab="Density", main=sprintf("Cell Health A549 similarity density"))
lines(density(unlist(balancedSample(lincs1Reps$cosReps$same[L1Perts], 100)), bw=0.01), col="orange", lwd=1.5, lty=5)
lines(density(lincs1Reps$inProdReps$diff, bw=0.01), col="blue", lwd=1.5)
lines(density(unlist(balancedSample(lincs1Reps$inProdReps$same[L1Perts], 100)), bw=0.01), col="firebrick3", lwd=1.5)
legend(x="topright", legend=c("Metric learning replicates", "Metric Learning all pairs", 
                              "Cosine replicates", "Cosine all pairs"), lwd=2.5, col=c("firebrick3", "blue", "orange", "forestgreen"), 
       lty=c(1,1,5,5))
dev.off()

pdf(file=file.path(outdir, "fig3a3_lincs2_simpdfs.pdf"), width=8, height=6)
plot(density(lincs2Reps$cosReps$diff, bw=0.01), col="forestgreen", lwd=1.5, lty=5, xlab="Similarity", ylab="Density", main=sprintf("Cell Health Batch 2 similarity density"))
lines(density(unlist(balancedSample(lincs2Reps$cosReps$same[L2Perts], 100)), bw=0.01), col="orange", lwd=1.5, lty=5)
lines(density(lincs2Reps$inProdReps$diff, bw=0.01), col="blue", lwd=1.5)
lines(density(unlist(balancedSample(lincs2Reps$inProdReps$same[L2Perts], 100)), bw=0.01), col="firebrick3", lwd=1.5)
legend(x="topright", legend=c("Metric learning replicates", "Metric Learning all pairs", 
                              "Cosine replicates", "Cosine all pairs"), lwd=2.5, col=c("firebrick3", "blue", "orange", "forestgreen"), 
       lty=c(1,1,5,5))
dev.off()


dsList <- list(lincs2Reps_A549, lincs2Reps_MCF7, lincs2Reps_U2OS)
dsnames <- c("LINCS2 A549", "LINCS2 MCF7", "LINCS2 U2OS")

pdf(file=file.path(outdir, "fig3a4_lincs2_simpdfs_bycell.pdf"), width=8, height=6)
for (ii in seq_along(dsList)){
  myds <- dsList[[ii]]
  dsname <- dsnames[ii]
  print(dsname)
  plot(density(myds$cosReps$diff, bw=0.01), col="forestgreen", lwd=1.5, lty=5, xlab="Similarity", ylab="Density", 
       main=sprintf("Cell Health Batch 2 similarity density, %s", dsname))
  lines(density(unlist(balancedSample(myds$cosReps$same[L2Perts], 100)), bw=0.01), col="orange", lwd=1.5, lty=5)
  lines(density(myds$inProdReps$diff, bw=0.01), col="blue", lwd=1.5)
  lines(density(unlist(balancedSample(myds$inProdReps$same[L2Perts], 100)), bw=0.01), col="firebrick3", lwd=1.5)
  legend(x="topright", legend=c("Metric learning replicates", "Metric Learning all pairs", 
                                "Cosine replicates", "Cosine all pairs"), lwd=2.5, col=c("firebrick3", "blue", "orange", "forestgreen"), 
         lty=c(1,1,5,5))
}
dev.off()


# Fig 3b: Cross validation

fdrdf <- rbind(data.frame(dataset="CDRP Replicate fdr=0.05", method="ML", sp=brayFDRML$fdr05), 
               data.frame(dataset="CDRP Replicate fdr=0.05", method="Cos", sp=brayFDRCos$fdr05),
               data.frame(dataset=("CDRP MoA fdr=0.05"), method=moaCPdf$method[1:2], sp=moaCPdf[1:2, 4]), 
               data.frame(dataset=("CDRP MoA fdr=0.1"), method=moaCPdf$method[1:2], sp=moaCPdf[1:2, 5]), 
               data.frame(dataset=("CDRP MoA fdr=0.25"), method=moaCPdf$method[1:2], sp=moaCPdf[1:2, 6]))
fdrdf$dataset <- factor(fdrdf$dataset, levels=c("CDRP Replicate fdr=0.05", "CDRP MoA fdr=0.05", 
                                                "CDRP MoA fdr=0.1", "CDRP MoA fdr=0.25"))

aucdf <- rbind(data.frame(dataset="CDRP Replicate", method="ML", auROC=brayAUCML), 
               data.frame(dataset="CDRP Replicate", method="Cos", auROC=brayAUCCos),
               moaCPdf[1:2, c(1,2,3)])

pdf(file=file.path(outdir, "fig3_CPAUCsummary.pdf"), width=5, height=6)
ggplot(aucdf, aes(x=factor(dataset, levels=c("CDRP Replicate", "CDRP MoA")), y=auROC, fill=method)) + geom_bar(stat="identity", position="dodge") + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05)) + theme_minimal() + ylab("Mean Balanced AUC") + 
  xlab("") + ggtitle("Cell Painting CDRP AUC") + ylim(c(0,1))
dev.off()

pdf(file=file.path(outdir, "fig3_CPFDRsummary.pdf"), width=8, height=6)
ggplot(fdrdf, aes(x=dataset, y=sp, fill=method)) + geom_bar(stat="identity", position="dodge") + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05)) + theme_minimal() + ylab("Statistical Power") + 
  xlab("") + ggtitle("Cell Painting CDRP Statistical Power ") + ylim(c(0,1))
dev.off()

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


##### moaCPdf Plot #####
pdf(file.path(outdir, "fig3S_CPMoAAUROC_bar.pdf"), width=8, height=6)
ggplot(moaCPdf, aes(x=dataset, y=auROC, fill=method)) + geom_bar(stat="identity", position="dodge") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + xlab("Cell Painting dataset") + ylab("AUROC") + 
  ggtitle("Mechanism of Action AUROC for cell painting datasets") + coord_cartesian(ylim=c(0.4,1))
dev.off()

pdf(file.path(outdir, "fig3S_CPMoAfdr10_bar.pdf"), width=8, height=6)
ggplot(moaCPdf, aes(x=dataset, y=fdr10, fill=method)) + geom_bar(stat="identity", position="dodge") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + xlab("Cell Painting dataset") + ylab("Balanced FDR < 0.1") + 
  ggtitle("Mechanism of Action FDR < 0.1 for cell painting datasets")
dev.off()

#### MoA Scatter for Cell Painting Datasets, Mean Rank and Mean FDR < 0.1
# Refactor this into a function:
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
                        pcl=pclds$pcldf$pclname[1:length(brayMoA$mlPCLs$setSims)],
                        mlBrayvals = sapply(brayMoA$mlRanks, FUN=mean), 
                        cosBrayvals = sapply(brayMoA$cosRanks, FUN=mean), 
                        mlBrayFDR = sapply(brayMoA$mlRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.1)), 
                        cosBrayFDR = sapply(brayMoA$cosRanks, FUN=function(y) mean(p.adjust(y, "fdr") < 0.1)),
                        mlSNR=(sapply(brayMoA$mlPCLs$setSims, mean) - mean(brayMoA$mlPCLs$allSims))/sd(brayMoA$mlPCLs$allSims),
                        cosSNR=(sapply(brayMoA$cosPCLs$setSims, mean) - mean(brayMoA$cosPCLs$allSims))/sd(brayMoA$cosPCLs$allSims),
                        pclSize=sapply(brayMoA$mlPCLs$setSims, length))

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

pdf(file.path(outdir, "fig4_CPBray_MoASNR_scatter.pdf"), width=8, height=6)
ggplot(BrayMoADF[BrayMoADF$pclSize > 50, ], aes(x=cosSNR, y=mlSNR)) + geom_point(aes(size=log10(pclSize)), color="blue") + 
  theme_minimal() + geom_abline(slope=1, intercept=0, linetype="dashed", color="black") + 
  geom_text_repel(aes(label=pcl, size=2)) + xlab("Cosine SNR") + ylab("Metric learning SNR") + 
  ggtitle("CDRP MoA SNR") + xlim(c(0,5)) + ylim(c(0,5))
dev.off()

ggplot(pclsnrDF[pclsnrDF$size > 50,], aes(x=Cos, y=ML, color=cellid)) + geom_point(aes(size=log10(size)), alpha=0.8) + 
  theme_minimal() + geom_abline(slope=1, intercept=0, linetype="dashed", color="black") + 
  geom_text_repel(aes(label=pcl, size=1.5)) + xlab("Cosine SNR") + ylab("Metric learning SNR") + 
  ggtitle("L1K MoA SNR")



# Supp Fig: doses
pdf(file=file.path(outdir, "fig3S_CPDoses.pdf"), width=8, height=6)
plot(density(log10(brayds$metads$Metadata_mmoles_per_liter + 1e-3), bw=0.03), col=cbbPalette[2], lwd=3, xlab="Log10 Dose uM", ylab="Density", main="Cell Painting Doses, uM")
lines(density(log10(lincsds1$metads$Metadata_mmoles_per_liter + 1e-3), bw=0.03), col=cbbPalette[4], lwd=3)
lines(density(log10(lincsds2$metads$Metadata_mmoles_per_liter + 1e-3), bw=0.03), col=cbbPalette[5], lwd=3)
legend(x="topleft", legend=c("CDRP 2017", "Cell Health 1", "Cell Health 2"), col=cbbPalette[c(2, 4, 5)], lwd=3)
dev.off()