source("figinit.R")
source("../plotFunctions.R")

#### Figure 4: How much data do you need to learn a useful representation? ####
# L1000

if (!file.exists(file.path(outdir, "../figspaper_res/ntrainingRes.rds"))){
  
  ntraindir <- "~/Work/bhk/analysis/metric_learning/2022_L1K/modelRuns"
  myf <- list.files(ntraindir, pattern="ntraining")
  myf <- setdiff(myf, "L1Kntraining_epch=15_cell=HCC515.rds")
  fcells <- sapply(strsplit(myf, split="cell="), FUN=function(x) strsplit(x[2], ".rds")[[1]])
  
  ntrainL1K <- lapply(seq_along(myf), FUN=function(x) 
    makeNTrainingPlots(ntrainingResPath = file.path(ntraindir, myf[x]), saveFiles=1, fname=fcells[x], mytitle=fcells[x], outpath=file.path(ntraindir, "nTrainingFigs")))
  
  resNorm <- bind_rows(lapply(seq_along(myf), FUN=function(x) ntrainL1K[[x]]$resNorm))
  
  # Cell Painting
  cpdir <- "~/Work/bhk/analysis/metric_learning/2022_cellpaint/"
  braydir <- file.path(cpdir, "bray/models")
  lincsdir <- file.path(cpdir, "lincs/models")
  
  brayntrain <- makeNTrainingPlots(ntrainingResPath= file.path(braydir, "brayntraining_epch=10.rds"), saveFiles=0, fname="bray", mytitle="CDRP", 
                                   outpath=file.path(braydir, "ntrainFigs"))
  
  lincsntrainf <- list.files(lincsdir, pattern="ntraining")
  lincsCells <- c("Batch1", "Batch2")
  
  lincsntrain <- lapply(seq_along(lincsntrainf), FUN=function(x) 
    makeNTrainingPlots(ntrainingResPath = file.path(lincsdir, lincsntrainf[x]), saveFiles=1, fname=lincsCells[x], mytitle=lincsCells[x], 
                       outpath=file.path(lincsdir, "ntrainFigs")))
  
  cpntrain <- c(lincsntrain, list(brayntrain))
  cpNorm <- bind_rows(lapply(seq_along(cpntrain), FUN=function(x) cpntrain[[x]]$resNorm))
  
  combResNorm <- rbind(resNorm, brayntrain$resNorm)
  combResNorm$dsname[combResNorm$dsname == "all"] <- "L1000 all"
  
  saveRDS(list(resNorm=resNorm, 
               brayntrain=brayntrain,
               lincsntrain=lincsntrain,
               combResNorm=combResNorm),
          file.path(outdir, "../figspaper_res/ntrainingRes.rds"))
} else {
  ntrainRes <- readRDS(file.path(outdir, "../figspaper_res/ntrainingRes.rds"))
  attach(ntrainRes)
}

# L1000 Figures
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


# Bray only while I debug LINCS
pdf(file=file.path(outdir, "fig4_BrayNormValidLoss.pdf"), width=8, height=6)
ggplot(brayntrain$resNorm, aes(x=log2(trainClasses+1), y=validAvgLoss, color=dsname)) + geom_point() + geom_smooth(method="loess") + theme_minimal() + 
  xlab("Log2 N of Training Compounds") + ylab("Normalized Average Validation Loss by compound") + ggtitle("Cell Painting Validation Loss vs Number of Training compounds") + 
  geom_hline(yintercept=-1, color="black", linetype="dashed") + scale_color_discrete(name="dataset")
dev.off()

pdf(file=file.path(outdir, "fig4_BraydeltaValidAURank.pdf"), width=8, height=6)
ggplot(brayntrain$resNorm, aes(x=log2(trainClasses+1), y=validAURank, color=dsname)) + geom_point() + geom_smooth(method="loess") + theme_minimal() + 
  xlab("Log2 N of Training Compounds") + ylab("Delta Average Validation AURank by compound") + ggtitle("Cell Painting Change in Validation AURank vs Number of Training compounds") + 
  geom_hline(yintercept=0, color="black", linetype="dashed") + scale_color_discrete(name="dataset")
dev.off()


#### Combined figure
pdf(file=file.path(outdir, "fig4S_ntrainingCombined_ValidLoss.pdf"), width=7, height=7)
ggplot(combResNorm, aes(x=trainClasses+1, y=validAvgLoss, color=dsname)) + geom_point() + geom_smooth(method="loess") + theme_minimal() + 
  xlab("Number of Training Compounds + 1") + ylab("Normalized Validation Loss") + ggtitle("Validation Loss vs Number of Training compounds") + 
  geom_hline(yintercept=-1, color="black", linetype="dashed") + scale_color_discrete(name="dataset") + scale_x_log10() + 
  annotate('rect', xmin=1, xmax=4500, ymin=-2.5, ymax=-1, alpha=0.1, fill="blue") + 
  annotate('rect', xmin=1, xmax=4500, ymin=-1, ymax=0, alpha=0.1, fill="red") + 
  theme(legend.position="bottom")
dev.off()

pdf(file=file.path(outdir, "fig4_ntrainingCombined_auROC.pdf"), width=7, height=7)
ggplot(combResNorm, aes(x=trainClasses+1, y=validAURank, color=dsname)) + geom_point() + geom_smooth(method="loess") + theme_minimal() + 
  xlab("Number of Training Compounds + 1") + ylab("Delta Validation auROC") + ggtitle("Change in auROC vs Number of Training compounds") + 
  geom_hline(yintercept=0, color="black", linetype="dashed") + scale_color_discrete(name="dataset") + scale_x_log10() + 
  annotate('rect', xmin=1, xmax=4500, ymin=0, ymax=0.18, alpha=0.1, fill="blue") + 
  annotate('rect', xmin=1, xmax=4500, ymin=-0.18, ymax=0, alpha=0.1, fill="red") + 
  theme(legend.position="bottom")
dev.off()

