source("figinit.R")



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