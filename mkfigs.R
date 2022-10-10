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


#### Fig 2: Cosine Theory
# Figure 2b: see bhk/code/cosine/cosineScript.R
x10 <- getCosineDist(mydim=20, N=1000, mkfig=0)
x100 <- getCosineDist(mydim=100, N=1000, mkfig=0)
x1000 <- getCosineDist(mydim=500, N=1000, mkfig=0)

# Fig 2a:
pdf(file.path(outdir, "fig2a_cosineVsDim_dists.pdf"), width=7, height=5)
plot(density(x1000[upper.tri(x1000)], bw=0.01), col="forestgreen", lwd=3, xlim=c(-1,1), 
     main="Distribution of cosine similarities for standard normals of dimension N")
lines(density(x10[upper.tri(x10)], bw=0.01), col="red", lwd=3)
lines(density(x100[upper.tri(x100)], bw=0.01), col="blue", lwd=3)
legend(x="topleft", legend=c("N=20", "N=100", "N=500"), lwd=c(4,4,4), col=c("red", "darkblue", "forestgreen"))
dev.off()

# Supp Fig: 
pdf(file.path(outdir, "fig2e_cosineFit_supp.pdf"), width=7, height=5)
x100 <- getCosineDist(mydim=50, N=1000, mkfig=1)
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
