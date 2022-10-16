library(gridExtra)
library(ggplot2)


outdir <- "./figspaper"

#### Fig 1 workflow - matrix
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


#### Fig 6: Cosine Theory
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
l1kmeta <- read_l1k_meta(datapath, version=2020)
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




########### Committee Figures

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


