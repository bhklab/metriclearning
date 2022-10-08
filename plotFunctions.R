

#### Plot functions for nTrainingExperiments
makeNTrainingPlots <- function(ntrainingResPath, saveFiles=1, fname="myNTrainingRes", mytitle="MyData", outpath="."){
  res <- readRDS(ntrainingResPath)
  
  if (saveFiles){
    pdf(file.path(outpath, sprintf("%s_LossVsTrainingSize.pdf", fname)), width=9, height=7)
    plot(log2(res$trainClasses+1 ), res$trainLoss, pch=16, col=alpha("blue",0.6), xlab="Log2 Number of Training Compounds + 1", ylab="T Loss", 
         main=sprintf("%s Training and Validation Loss vs Number of Training Compounds", mytitle), 
         ylim=c(min(c(res$trainLoss, res$validLoss)), max(c(res$trainLoss, res$validLoss))))
    points(log2(res$trainClasses+1), res$validLoss, pch=16, col=alpha("red",0.6))
    legend(x="bottomright", legend=c("Training Loss", "Validation Loss"), pch=16, col=c("blue", "red"))
    lines(x=c(-10, 1000), y=c(mean(res$validLoss[res$trainClasses == 0]), mean(res$validLoss[res$trainClasses == 0])), col="black", lty=2)
    dev.off()
    
    pdf(file.path(outpath, "a375_AURankVsTrainingSize.pdf"), width=9, height=7)
    plot(log2(res$trainClasses + 1), res$validAURank, pch=16, col=alpha("red", 0.6), xlab="Log2 Number of Training Compounds+1", ylab="AU Rank",
         main="A375 Validation AU Rank vs Number of Training Compounds")
    lines(x=c(-10, 1000), y=c(mean(res$validAURank[res$trainClasses == 0]), mean(res$validAURank[res$trainClasses == 0])), col="black", lty=2)
    dev.off()
  } else {
    
    
  }
}