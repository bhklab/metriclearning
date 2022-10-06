


# test metricNTraining
datapath <- "~/Work/bhk/data/l1k/2020/"
outpath <- "~/Work/bhk/analysis/metric_learning/2022_L1K/figs"
cell_id <- "A375"

l1k_meta <- read_l1k_meta(datapath, version=2020)
mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$cell_id == cell_id & l1k_meta$siginfo$pert_type == "trt_cp",]

ds <- parse_gctx(get_level5_ds(datapath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)

resa375 <- metricNTraining(t(ds@mat), mysigs$pert_iname, epochs=2, reps=1)



# Cell Line Specificity

# See how many unique perts each cell line has
# Get those cell lines with at least 1000 unique compounds assayed:
upertcnt <- sapply(names(sort(table(l1kmeta$siginfo$cell_id[l1kmeta$siginfo$pert_type == "trt_cp"]), decreasing=TRUE)), 
                   FUN=function(x) length(unique(siginfo$pert_iname[siginfo$cell_id == x & siginfo$pert_type == "trt_cp"])))

sort(upertcnt[upertcnt >= 1000], decreasing=TRUE)

# Look up the cell lineage of these cell lines:
cellinfo[match(names(upertcnt[upertcnt >= 1000]), cellinfo$cell_id), c(1, 10, 14, 15, 16, 18)]

mycellids <- c("HEPG2", "HCC515", "NPC", "ASC", "HEK293", "YAPC")

grpprts <- unique(siginfo$pert_iname[siginfo$cell_id %in% mycellids & siginfo$pert_type == "trt_cp"])
trainprts <- sample(grpprts, round(length(grpprts)/2))
testprts <- setdiff(grpprts, trainprts)


# Test Cell Line Specificity

clsres_fix <- analyzeL1KCellLineSpecificity(datapath, datapath, cell_ids=mycellids, outpath=".", ncpds=4000, epochs=10, iter=3)
  
clsres_fix$trainSet <- factor(clsres_fix$trainSet, levels=c("cosine", "all", "cell"))
pdf(file.path(outpath, "CellSpec_TestLoss.pdf"), width=9, height=7)
ggplot(clsres_fix, aes(fill=cell_id, y=testAvgLoss, x=trainSet)) + geom_bar(position="dodge", stat="identity") + xlab("Training Set") + ylab("Test Set Loss") + 
  theme_minimal() + ggtitle("Cell-specific test loss")
dev.off()

pdf(file.path(outpath, "CellSpec_TestAURank.pdf"), width=9, height=7)
ggplot(clsres_fix, aes(fill=cell_id, y=testAURank-0.5, x=trainSet)) + geom_bar(position="dodge", stat="identity") + xlab("Training Set") + ylab("Test Set AURank") + 
  theme_minimal() + ggtitle("Cell-specific test loss") + scale_y_continuous(limits=c(0, 0.5), labels=seq(0.5, 1, 0.1))
dev.off()

pdf(file.path(outpath, "ReplicateCosineSimilarityDist.pdf"), width=9, height=7)
plot(density(mysim$diff, bw=0.01), col="blue", lwd=3, xlab="Similarity", ylab="Density", main="Cosine similarity for A375 compounds, Loss = -1.721", xlim=c(-1, 1))
lines(density(unlist(mysim$same), bw=0.01), col="red", lwd=3)
legend(x="topright", legend=c("Replicates", "All Signatures"), col=c("red", "blue"), lwd=c(3,3))
dev.off()

pdf(file.path(outpath, "ReplicateAURank.pdf"), width=9, height=7)
plot(ecdf(rankVectors(unlist(mysim$same), mysim$diff)), col="black", lwd=2, xlab="Rank relative to background", ylab="Cumulative fraction", main="A375 Cosine Rank, mean=0.772")
dev.off()





##### nTraining
asc <- readRDS("L1Kntraining_epch=20_cell=ASC.rds")

pdf(file.path(outpath, "asc_LossVsTrainingSize.pdf"), width=9, height=7)
plot(asc$trainClasses, -log(-asc$validLoss), pch=16, col=alpha("red", 0.6), ylim=c(-4, 1), xlab="Number of Training Compounds", ylab="-Log -Loss", 
     main="ASC Training and Validation Loss vs Number of Training Compounds")
points(asc$trainClasses, -log(-asc$trainLoss), pch=16, col=alpha("blue", 0.6))
legend(x="bottomright", legend=c("Training Loss", "Validation Loss"), pch=16, col=c("blue", "red"))
lines(x=c(-10, 250), y=c(-0.2713, -0.2713), col="black", lty=2)
dev.off()

a375 <- readRDS("L1Kntraining_epch=20_cell=A375.rds")
hcc <- readRDS("L1Kntraining_epch=20_cell=HCC515.rds")

pdf(file.path(outpath, "a375_LossVsTrainingSize.pdf"), width=9, height=7)
plot(log2(a375$trainClasses+1 ), a375$trainLoss, pch=16, col=alpha("blue",0.6), xlab="Log2 Number of Training Compounds + 1", ylab="T Loss", 
     main="A375 Training and Validation Loss vs Number of Training Compounds")
points(log2(a375$trainClasses+1), a375$validLoss, pch=16, col=alpha("red",0.6))
legend(x="bottomright", legend=c("Training Loss", "Validation Loss"), pch=16, col=c("blue", "red"))
lines(x=c(-10, 1000), y=c(mean(a375$validLoss[a375$trainClasses == 0]), mean(a375$validLoss[a375$trainClasses == 0])), col="black", lty=2)
dev.off()

pdf(file.path(outpath, "a375_AURankVsTrainingSize.pdf"), width=9, height=7)
plot(log2(a375$trainClasses + 1), a375$validAURank, pch=16, col=alpha("red", 0.6), xlab="Log2 Number of Training Compounds+1", ylab="AU Rank",
       main="A375 Validation AU Rank vs Number of Training Compounds")
lines(x=c(-10, 1000), y=c(mean(a375$validAURank[a375$trainClasses == 0]), mean(a375$validAURank[a375$trainClasses == 0])), col="black", lty=2)
dev.off()

pdf(file.path(outpath, "hcc515_LossVsTrainingSize.pdf"), width=9, height=7)
plot(log2(hcc$trainClasses+1 ), hcc$trainLoss, pch=16, col=alpha("blue",0.6), xlab="Log2 Number of Training Compounds + 1", ylab="T Loss", 
     main="HCC515 Training and Validation Loss vs Number of Training Compounds", ylim=c(min(c(hcc$trainLoss, hcc$validLoss)), max(c(hcc$trainLoss, hcc$validLoss))))
points(log2(hcc$trainClasses+1), hcc$validLoss, pch=16, col=alpha("red",0.6))
legend(x="bottomright", legend=c("Training Loss", "Validation Loss"), pch=16, col=c("blue", "red"))
lines(x=c(-10, 1000), y=c(mean(hcc$validLoss[hcc$trainClasses == 0]), mean(hcc$validLoss[hcc$trainClasses == 0])), col="black", lty=2)
dev.off()

pdf(file.path(outpath, "hcc515_AURankVsTrainingSize.pdf"), width=9, height=7)
plot(log2(hcc$trainClasses + 1), hcc$validAURank, pch=16, col=alpha("red", 0.6), xlab="Log2 Number of Training Compounds+1", ylab="AU Rank",
     main="HCC515 Validation AU Rank vs Number of Training Compounds", ylim=c(0.5,1))
lines(x=c(-10, 1000), y=c(mean(hcc$validAURank[hcc$trainClasses == 0]), mean(hcc$validAURank[hcc$trainClasses == 0])), col="black", lty=2)
dev.off()
