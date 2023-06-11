
#### Cell Specificity data ####

# The x_PCL_pcl_res files look like they are calls to benchmarkL1KCellModel; check to confirm they are
# applications of the latest trained metrics. 

specdir <- file.path(topdir, "modelRuns/PCLFigs")
mycells <- c("A375", "A549", "ASC", "HCC515", "HEK293", "HEPG2", "MCF7", "NPC", "PC3", "BT20", "HUVEC")
pclf <- list.files(specdir, pattern="PCL_pcl_res.rds")

#if (!file.exists(file.path(outdir, "../figspaper_res/cpSummaryRes.rds"))){
if (!file.exists(file.path(outdir, "../figspaper_res/cellSpecificRes.rds"))){
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
  
  saveRDS(list(specEvalDF=specEvalDF, specPlotDF1=specPlotDF1, specPlotDF2=specPlotDF2), 
          file.path(outdir, "../figspaper_res/cellSpecificRes.rds"))
  
} else {
  cellSpecRes <- readRDS(file.path(outdir, "../figspaper_res/cellSpecificRes.rds"))
  attach(cellSpecRes)
}


pdf(file=file.path(outdir, "fig4_cellSpecificAll_AUROC.pdf"), width=7, height=3.5)
ggplot(specPlotDF2, aes(x=testset, y=value, fill=as.factor(trainset))) + geom_bar(position="dodge", stat="identity") + 
  theme_minimal() + xlab("Cell Line") + ylab("Balanced auROC") + ggtitle("Context Specificity: MoA auROC") + 
  scale_fill_discrete(name = "Similarity Function") + coord_cartesian(ylim=c(0.5,0.8)) + 
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))
dev.off()

pdf(file=file.path(outdir, "fig4_cellSpecificAll_FDR05.pdf"), width=7, height=3.5)
ggplot(specPlotDF1, aes(x=testset, y=value, fill=as.factor(trainset))) + geom_bar(position="dodge", stat="identity") + 
  theme_minimal() + xlab("Cell Line") + ylab("MoA Recall at FDR = 0.05") + 
  ggtitle("Context Specificity: MoA Recall at FDR = 0.05") + scale_fill_discrete(name = "Similarity Function") +
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))
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