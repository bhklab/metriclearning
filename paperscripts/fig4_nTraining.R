source("figinit.R")

#### Figure 4: How much data do you need to learn a useful representation? ####
# L1000
ntraindir <- "~/Work/bhk/analysis/metric_learning/2022_L1K/modelRuns"
f <- list.files(ntraindir, pattern="ntraining")
f <- setdiff(f, "L1Kntraining_epch=15_cell=HCC515.rds")
fcells <- sapply(strsplit(f, split="cell="), FUN=function(x) strsplit(x[2], ".rds")[[1]])

ntrainL1K <- lapply(seq_along(f), FUN=function(x) 
  makeNTrainingPlots(ntrainingResPath = file.path(ntraindir, f[x]), saveFiles=1, fname=fcells[x], mytitle=fcells[x], outpath=file.path(ntraindir, "nTrainingFigs")))

resNorm <- bind_rows(lapply(seq_along(f), FUN=function(x) ntrainL1K[[x]]$resNorm))

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

brayntrain <- makeNTrainingPlots(ntrainingResPath= file.path(braydir, "brayntraining_epch=10.rds"), saveFiles=0, fname="bray", mytitle="CDRP", 
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
combResNorm <- rbind(resNorm, brayntrain$resNorm)
combResNorm$dsname[combResNorm$dsname == "all"] <- "L1000 all"

pdf(file=file.path(outdir, "fig4S_ntrainingCombined_ValidLoss.pdf"), width=8, height=6)
ggplot(combResNorm, aes(x=trainClasses+1, y=validAvgLoss, color=dsname)) + geom_point() + geom_smooth(method="loess") + theme_minimal() + 
  xlab("Number of Training Compounds + 1") + ylab("Normalized Validation Loss") + ggtitle("Validation Loss vs Number of Training compounds") + 
  geom_hline(yintercept=-1, color="black", linetype="dashed") + scale_color_discrete(name="dataset") + scale_x_log10() + 
  annotate('rect', xmin=1, xmax=4500, ymin=-2.5, ymax=-1, alpha=0.1, fill="blue") + 
  annotate('rect', xmin=1, xmax=4500, ymin=-1, ymax=0, alpha=0.1, fill="red")
dev.off()

pdf(file=file.path(outdir, "fig4_ntrainingCombined_auROC.pdf"), width=8, height=6)
ggplot(combResNorm, aes(x=trainClasses+1, y=validAURank, color=dsname)) + geom_point() + geom_smooth(method="loess") + theme_minimal() + 
  xlab("Number of Training Compounds + 1") + ylab("Delta Validation auROC") + ggtitle("Change in auROC vs Number of Training compounds") + 
  geom_hline(yintercept=0, color="black", linetype="dashed") + scale_color_discrete(name="dataset") + scale_x_log10() + 
  annotate('rect', xmin=1, xmax=4500, ymin=0, ymax=0.18, alpha=0.1, fill="blue") + 
  annotate('rect', xmin=1, xmax=4500, ymin=-0.18, ymax=0, alpha=0.1, fill="red")
dev.off()






#### Cell Specificity data ####

# The x_PCL_pcl_res files look like they are calls to benchmarkL1KCellModel; check to confirm they are
# applications of the latest trained metrics. 

specdir <- file.path(topdir, "modelRuns/PCLFigs")
mycells <- c("A375", "A549", "ASC", "HCC515", "HEK293", "HEPG2", "MCF7", "NPC", "PC3", "BT20", "HUVEC")
pclf <- list.files(specdir, pattern="PCL_pcl_res.rds")

specEvalDF <- data.frame()

for (acell in mycells){
  print(acell)
  a <- readRDS(file.path(specdir, sprintf("%s_PCL_pcl_res.rds", acell)))
  
  Amlbrk <- 1-mean(unlist(balancedSample(a$MLRanks, 1000)))
  Acosbrk <- 1-mean(unlist(balancedSample(a$cosranks, 1000)))
  
  Amlfdr05 <- mean(p.adjust(unlist(balancedSample(a$MLRanks, k=1000)), method="fdr") < 0.05)
  Acosfdr05 <- mean(p.adjust(unlist(balancedSample(a$cosranks, k=1000)), method="fdr") < 0.05)
  
  specEvalDF <- rbind(specEvalDF, data.frame(trainset="ML cell", testset=acell, 
                                       auMLrank=(1-mean(unlist(a$MLRanks))), 
                                       auCosRank=(1-mean(unlist(a$cosranks))), 
                                       auMLBalrank=Amlbrk, auCosBalRank=Acosbrk, 
                                       mlFDR05=Amlfdr05, cosFDR05=Acosfdr05))
  
  b <- readRDS(file.path(specdir, sprintf("all_%s_PCL_pcl_res.rds", acell)))
  
  Bmlbrk <- 1-mean(unlist(sapply(b$MLRanks, FUN=function(x) sample(x, min(1000, length(x))))))
  Bcosbrk <- 1-mean(unlist(sapply(b$cosranks, FUN=function(x) sample(x, min(1000, length(x))))))
  
  Bmlfdr05 <- mean(p.adjust(unlist(balancedSample(b$MLRanks, k=1000)), method="fdr") < 0.05)
  Bcosfdr05 <- mean(p.adjust(unlist(balancedSample(b$cosranks, k=1000)), method="fdr") < 0.05)
  
  
  specEvalDF <- rbind(specEvalDF, data.frame(trainset="ML all", testset=acell, 
                                       auMLrank=(1-mean(unlist(b$MLRanks))), 
                                       auCosRank=(1-mean(unlist(b$cosranks))), 
                                       auMLBalrank=Bmlbrk, auCosBalRank=Bcosbrk, 
                                       mlFDR05=Bmlfdr05, cosFDR05=Bcosfdr05))
}


specPlotDF1 <- rbind(reshape2::melt(specEvalDF, id.vars=c("trainset", "testset"), measure.vars=c("mlFDR05")), 
                     reshape2::melt(specEvalDF[specEvalDF$trainset == "ML all",], id.vars=c("trainset", "testset"), measure.vars=c("cosFDR05"))) 
specPlotDF1$trainset[specPlotDF1$variable == "cosFDR05"] <- "cosine"

specPlotDF2 <- rbind(reshape2::melt(specEvalDF, id.vars=c("trainset", "testset"), measure.vars=c("auMLBalrank")), 
                     reshape2::melt(specEvalDF[specEvalDF$trainset == "ML all",], id.vars=c("trainset", "testset"), measure.vars=c("auCosBalRank"))) 
specPlotDF2$trainset[specPlotDF2$variable == "auCosBalRank"] <- "cosine"


pdf(file=file.path(outdir, "fig4_cellSpecificAll_AUROC.pdf"), width=8, height=6)
ggplot(specPlotDF2, aes(x=testset, y=value, fill=as.factor(trainset))) + geom_bar(position="dodge", stat="identity") + theme_minimal() + 
  xlab("Cell Line") + ylab("Balanced auROC") + ggtitle("Context Specificity: MoA auROC") + scale_fill_discrete(name = "Similarity Function") + coord_cartesian(ylim=c(0.5,0.8))
dev.off()

pdf(file=file.path(outdir, "fig4_cellSpecificAll_FDR05.pdf"), width=8, height=6)
ggplot(specPlotDF1, aes(x=testset, y=value, fill=as.factor(trainset))) + geom_bar(position="dodge", stat="identity") + theme_minimal() + 
  xlab("Cell Line") + ylab("MoA Recall at FDR = 0.05") + ggtitle("Context Specificity: MoA Recall at FDR = 0.05") + scale_fill_discrete(name = "Similarity Function")
dev.off()

pdf(file=file.path(outdir, "fig4_cellSpecificVsCos_AUROC.pdf"), width=8, height=6)
ggplot(specEvalDF, aes(x=testset, y=(auMLBalrank - auCosBalRank), fill=as.factor(trainset))) + geom_bar(position="dodge", stat="identity") + theme_minimal() + 
  xlab("Cell Line") + ylab("Delta Balanced auROC relative to cosine") + ggtitle("Context Specificity: Improvement in MoA auROC") + scale_fill_discrete(name = "Training Set")
dev.off()

pdf(file=file.path(outdir, "fig4_cellSpecificVsCos_FDR05.pdf"), width=8, height=6)
ggplot(specEvalDF, aes(x=testset, y=(mlFDR05 - cosFDR05), fill=as.factor(trainset))) + geom_bar(position="dodge", stat="identity") + theme_minimal() + 
  xlab("Cell Line") + ylab("Delta MoA Recall at FDR = 0.05 compared to cosine") + ggtitle("Context Specificity: Improvement in MoA Recall") + scale_fill_discrete(name = "Training Set")
dev.off()




pdf(file.path("../../figs", "L1K_PCL_deltaRank_barplot.pdf"), width=10, height=6)
ggplot(specEvalDF, aes(x=testset, y=auMLrank - auCosRank, fill=as.factor(trainset))) + geom_bar(position="dodge", stat="identity") + theme_minimal() +
  xlab("Cell Line") + ylab("ML Rank - cos Rank") + ggtitle("Delta Mean Rank, PCL data") + scale_fill_discrete(name = "Training Set")
dev.off()
pdf(file.path("../../figs", "L1K_PCL_MLRank_barplot.pdf"), width=10, height=6)
ggplot(specEvalDF, aes(x=testset, y=auMLrank, fill=as.factor(trainset))) + geom_bar(position="dodge", stat="identity") + theme_minimal() +
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