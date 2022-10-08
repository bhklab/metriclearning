outdir <- "./figs"

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
