
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
