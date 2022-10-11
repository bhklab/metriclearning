

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
    
    pdf(file.path(outpath, sprintf("%s_AURankVsTrainingSize.pdf", fname)), width=9, height=7)
    plot(log2(res$trainClasses + 1), res$validAURank, pch=16, col=alpha("red", 0.6), xlab="Log2 Number of Training Compounds+1", ylab="AU Rank",
         main=sprintf("%s Validation AU Rank vs Number of Training Compounds", mytitle))
    lines(x=c(-10, 1000), y=c(mean(res$validAURank[res$trainClasses == 0]), mean(res$validAURank[res$trainClasses == 0])), col="black", lty=2)
    dev.off()
  } else {
    
    
  }
}


### Plot functions for xval experiments
makeXValPlots <- function(xvalResPath, saveFiles=1, fname="XValRes", mytitle="MyData", outpath="."){
  res <- readRDS(xvalResPath)
  
  minloss <- min(sapply(seq(2, length(res)), FUN=function(x) min(res[[x]]$mean_train_losses)))
  epochcount <- length(res$fold1$mean_train_losses)
  
  pdf(file.path(outpath, sprintf("%s_xvalLoss.pdf", fname)), width=8, height=6)
  plot(-10, -10, xlim=c(0, epochcount), ylim=c(-6, 0), #ylim=c(1.1*min(res$fold1$mean_train_losses), 0), 
       main=sprintf("%s Cross Validation Loss, nFolds=%d, epochs=%d", mytitle, length(res)-1, epochcount), xlab="Batch Epoch", ylab="Cosine T Loss")
  for (myfold in res[seq(2, length(res))]){
    trainloss <- c(myfold$baseline_train_loss, myfold$mean_train_losses)
    validloss <- c(myfold$baseline_valid_loss, myfold$mean_valid_losses)
    
    lines(seq(0, length(trainloss)-1), trainloss, col="forestgreen", lwd=2, type="b")
    lines(seq(0, length(validloss)-1), validloss, col="red", lwd=2, type="b")
  }
  legend(x="topright", legend=c("Training", "Test"), col=c("forestgreen", "red"), lwd=c(4,4))
  grid()
  dev.off()
  
  
  lossdf <- data.frame(fold=character, type=character, loss=numeric)
  for (ii in seq(2, length(res))-1){
    myfold <- res[[ii+1]]
    lossdf <- rbind(lossdf, data.frame(fold=ii, type="training", loss=myfold$train_losses))
    lossdf <- rbind(lossdf, data.frame(fold=ii, type="testing", loss=myfold$valid_losses))
  }
  
  pdf(file.path(outpath, sprintf("%s_xvalFinalLosses.pdf", fname)), width=8, height=6)
  print(ggplot(lossdf, aes(x=as.factor(fold), y=loss, fill=type)) + geom_boxplot() + theme_minimal() + xlab("Fold") + ylab("Loss") + 
    ggtitle(sprintf("%s Final Losses, nFolds=%d, epochs=%d", mytitle, length(res)-1, epochcount)))
  dev.off()
}

getXValBalRank <- function(xvalResPath, saveFiles=1, fname="XValRes", datapath=".",){
  
}



#### Plot functions for all datasets:

makeAllPlots <- function(allrespath, saveFiles=1, fname="allres", mytitle="MyData", outpath="."){
  
}