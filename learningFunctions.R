library(caret)

#' metricCrossValidate 
#' takes a matrix and class vector, runs n-fold cross validation with a learned inner product,
#' and reports the test and training loss with training.
#' @param mat1 Numeric, a matrix of size N x d, where rows are samples and columns are features
#' @param classes Factor (or character, numeric), where two elements with the same value are assumed to be similar
#' @param nFolds Numeric, number of folds to run (default = 5)
#' @param metric Character, denotes which type of metric to learn. One of {}.
#' @param epochs Numeric, the number of epochs of training to run (default = 100)

metricCrossValidate <- function(mat1, classes, nFolds=5, metric="", epochs=100, loss=mycos_t_loss){
# Revise to intelligently choose folds in the case of a wide range of class sizes
# But for now, do it naively
  
  mygroups <- table(classes)
  myfolds <- createFolds(names(mygroups), k=nFolds)
  myfolds <- sapply(myfolds, FUN=function(x) names(mygroups[x]))

  res_all <- list(myfolds=myfolds)
  mymodels <- list()
  
  for (testset in myfolds){
    testix <- which(classes %in% testset)
    trainix <- setdiff(seq_along(classes), testix)
    
    # res <- learnInnerProduct(mat1[trainix, ], classes[trainix], metric=metric, epochs=epochs)  
    
    gends_train <- genDataset(mat1[trainix,], classes[trainix], pca_first = FALSE, scale = FALSE, center = FALSE)
    train_dl <- dataloader(gends_train, batch_size = 1, shuffle = TRUE)
    
    gends_valid <- genDataset(mat1[testix,], classes[testix], pca_first=FALSE, scale=FALSE, center=FALSE)
    valid_dl <- dataloader(gends_valid, batch_size = 1, shuffle = TRUE)
    
    model <- OneLayerLinear(dim(mat1)[2], dim(mat1)[2])
    device <- if (cuda_is_available()) torch_device("cuda:0") else "cpu"
    model <- model$to(device = device)
    optimizer <- optim_adam(model$parameters, lr = 0.01)
    
    res <- train_function(model=model, 
                          train_dl=train_dl, 
                          valid_dl=valid_dl, 
                          myloss=loss, 
                          device=device, 
                          optimizer=optimizer, 
                          epochs = epochs)
    
    print("Computing baseline sim")
    trainGrpSim <- innerProductGroups("cosine", mat1[trainix,], classes[trainix], compact=1)
    validGrpSim <- innerProductGroups("cosine", mat1[testix,], classes[testix], compact=1)
    
    res$baseline_train_loss <- mean((mean(trainGrpSim$diff) - sapply(trainGrpSim$same, mean))/sd(trainGrpSim$diff))
    res$baseline_valid_loss <- mean((mean(validGrpSim$diff) - sapply(validGrpSim$same, mean))/sd(validGrpSim$diff))
  
    # Separate the final model from each fold to save the results
    res_all[sprintf("fold%d", length(res_all))] <- list(res[names(res) != "model"])
    mymodels <- c(mymodels, res$model)
  }
  
  return(list(res_all=res_all, mymodels=mymodels))
}



#' learnInnerProduct
#' This function takes an input matrix with class labels and learns an inner product metric 
#' maximizing similarity between elements of the same class label and dissimilarity of 
#' elements of different class labels. 
#' @param mat1 Numeric, a matrix of size N x d, where rows are samples and columns are features
#' @param classes Factor (or character, numeric), where two elements with the same value are assumed to be similar
#' @param metric Character, denotes which type of metric to learn. One of {}.
#' @param epochs Numeric, the number of epochs of training to run (default = 100)

learnInnerProduct <- function(mat1, classes, metric="", epochs=100, loss=mycos_t_loss){
  gends <- genDataset(mat1, classes, pca_first = FALSE, scale = FALSE, center = FALSE)
  
  train_dl <- dataloader(gends, batch_size = 1, shuffle = TRUE)
  
  model <- OneLayerLinear(dim(mat1)[2], dim(mat1)[2])
  
  device <- if (cuda_is_available()) torch_device("cuda:0") else "cpu"
  model <- model$to(device = device)
  
  optimizer <- optim_adam(model$parameters, lr = 0.01)
  
  res <- train_function(model=model, 
                        train_dl=train_dl, 
                        valid_dl=train_dl, 
                        myloss=loss, 
                        device=device, 
                        optimizer=optimizer, 
                        epochs = epochs, 
                        #save_pars = "w"
                        )
  
  # Return training loss?
  return(res)
}


#' innerProduct
#' This function computes the inner product similarity between two matrices given a trained model object
#' @param model A torch nn_module trained on a dataset using learnInnerProduct
#' @param mat1 Numeric, a matrix of size N x d1
#' @param mat2 Numeric, a matrix of size N x d2
#' 
#' @return matrix
#' @export
innerProduct <- function(model, mat1, mat2){
  
  if (is.character(model)){
    if (model == "cosine"){
      return(cosine(mat1, mat2))
    }
  }
  # Check for integrity of variables?
  m1 <- model(torch_tensor(mat1, dtype=torch_float()))
  m2 <- model(torch_tensor(mat2, dtype=torch_float()))
  
  return(cosine(as.matrix(m1), as.matrix(m2)))
}


#' cosine
#' 
#' This function computes the cosine distance between all pairs of columns of two input matrices
#' @param mat1 Numeric, a matrix of size d1 x N
#' @param mat2 Numeric, a matrix of size d2 x N
#'
#' @return matrix
#'
#' @export
cosine <- function(mat1, mat2){
  a <- mat1 %*% t(mat2)
  
  #mag1 <- sqrt(diag(mat1 %*% t(mat1)))
  #mag2 <- sqrt(diag(mat2 %*% t(mat2)))
  
  mag1 <- sqrt(sapply(seq_len(dim(mat1)[1]), FUN=function(j) sum(mat1[j,]^2)))
  mag2 <- sqrt(sapply(seq_len(dim(mat2)[1]), FUN=function(j) sum(mat2[j,]^2)))
  
  return(a/outer(mag1, mag2, "*"))
}


#' innerProductGroups
#' 
#' This function takes a model, a data matrix, and group labels as input and computes the similarity
#' between elements in the matrix that belong to the same class compared with the similarity between 
#' elements in the matrix for different classes. 
#' @param model Torch nn_module, ideally trained with learnInnerProduct
#' @param mat1 Numeric, a matrix of size N x d, where rows are samples and columns are features
#' @param classes Factor (or character, numeric), where two elements with the same value are assumed to be similar
#' @param compact Boolean (default 0); whether to use sampling to compute the interclass similarity. For large matrices 
#' of more than 10^4 samples, computing the interclass similarity can be very expensive. Compact samples 1e6 interclass 
#' similarities, but computes all intraclass similarities. 
#' 
#' @return list of same = list of intraclass similarities, where each group is one element in the list; diff = interclass
#' similarities
#' @export
innerProductGroups <- function(model, mat1, classes, compact=0){
  if (compact){
    # Use compact if the matrix is too large to effectively compute
    # Compact samples from the space of unlike similarities rather than computing the entire matrix
    
    # This could likely be accomplished more elegantly with dplyr
    samesim <- sapply(names(table(classes)[table(classes)>1]), FUN=function(x){
      ix <- which(classes == x)
      a <- innerProduct(model, mat1[ix, ], mat1[ix, ])
      a[upper.tri(a)]
    })
    
    jx <- sample(dim(mat1)[1], min(1000, dim(mat1)[1]))
    diffsim <- innerProduct(model, mat1[jx, ], mat1[jx, ])
    for (mygrp in names(table(classes)[table(classes) > 1])){
      ix <- which(jx %in% which(classes == mygrp))
      if (length(ix) > 1){
        diffsim[ix, ix] <- NA
      }
    }
    diffsim <- diffsim[upper.tri(diffsim)]
    diffsim <- diffsim[!is.na(diffsim)]
    return(list(same=samesim, diff=diffsim))
    
  } else {
    sims <- innerProduct(model, mat1, mat1)

    samesim <- sapply(names(table(classes)[table(classes) > 1]), FUN=function(x) {
      a <- sims[which(classes == x), which(classes == x)]
      a[upper.tri(a)]
      })
    
    diffsim <- sims
    for (mygroup in classes){
      ix <- which(classes == mygroup)
      diffsim[ix,ix] <- NA
    }
    
    diffsim <- diffsim[!is.na(diffsim)]
    return(list(same=samesim, diff=diffsim))
  }
}


#' innerProductPairwiseGroups
#' 
#' This function takes a model, a data matrix, group labels, and sets of group labels as input and 
#' computes the similarity between signatures (rows) in the matrix that belong to the same set of group
#' labels but have different group labels. The use case is that there are sets of groups that are similar (e.g. 
#' compounds with the same mechanism of action), and we seek the similarities of elements of the same set that
#' have different group labels. 

innerProductPairwiseGroups <- function(model, mat1, classes, sets, compact=0){
  # Compact samples a block from the full similarity matrix, which can be efficient if mat1 is very large.
  # Note that compact does not exclude same elements, i.e. it is null distribution that includes both 
  # similarities where the null hypothesis and where the alternate hypothesis is true. The presumption is that 
  # the rate of instances of the alternate hypothesis is sufficiently low that this isn't a problem. 
  if (compact){
    if (dim(mat1)[1] < 4000){
      sims <- innerProduct(model, mat1, mat1)
    } else {
      ix <- sample(dim(mat1)[1], 4000)
      sims <- innerProduct(model, mat1[ix,], mat1[ix,])
    }
    
    # Consider filtering same-groups and same-sets elements here
    
  } else {
    sims <- innerProduct(model, mat1, mat1)
    
    # Consider filtering same-groups and same-sets elements here
  }
  
  sims <- sims[upper.tri(sims)]
    
  setSims <- list()
  length(setSims) <- length(sets)
  
  for (ii in seq_along(sets)){
    myset <- sets[[ii]]
    # Check how many elements of classes are in my set. We require at least two:
    mygrps <- intersect(classes, myset)
    
    if (length(mygrps) > 1){
      jx <- which(classes %in% myset)
      
      # Only include elements whose group labels are not equal
      filt <- outer(classes[jx], classes[jx], '!=')
      filt <- filt[upper.tri(filt)]
      
      mysim <- innerProduct(model, mat1[jx,], mat1[jx,])
      setSims[[ii]] <- (mysim[upper.tri(mysim)])[filt]
    } 
  }
  
  if (!is.null(names(sets))){
    names(setSims) <- names(sets)
  }
  
  return(list(setSims=setSims, allSims=sims))
}


#' metricNTraining
#' 
#' This function takes a dataset, partitions it into a 20% test set and 80% train set, then learns a similarity
#' metric on subsets of the training set to evaluate generalization error on the test set as a function of 
#' training examples. 
#' @param mat1 Numeric, a matrix of size N x d, where rows are samples and columns are features
#' @param classes Factor (or character, numeric), where two elements with the same value are assumed to be similar
#' @param metric character (currently has no function)
#' @param epochs Numeric, default 10
#' @param nvals (optional) List of training sizes to test
#' @param validClassLabs (optional) list of class labels to explicitly define the test set
#' @param reps Numeric, number of times to run the sampling (default 1).
#' @param samplePts Numeric, number of training sizes to test (default 10).
#'
#' @return matrix
#'
#' @export
metricNTraining <- function(mat1, classes, metric="", epochs=10, nvals=c(), validClassLabs=c(), reps=1, samplePts=10){
  mygroups <- table(classes)[table(classes) > 2]
  bggroups <- table(classes)[table(classes) <= 2]

  res <- data.frame(trainClasses=numeric(), 
                    rep=numeric(), 
                    trainLoss=numeric(), 
                    trainAvgLoss=numeric(),
                    validLoss=numeric(), 
                    validAvgLoss=numeric(), 
                    validAURank=numeric())
  
  for (ii in seq_len(reps)){
    print(sprintf("ii = %d", ii))
  
    if (length(validClassLabs) == 0){
      # Designate 20% of the classes as a test set
      validClassLabs <- union(sample(names(mygroups), ceiling(length(mygroups)*0.2)), sample(names(bggroups), ceiling(length(bggroups)*0.2)))
    }  
    trainClassLabs <- setdiff(names(mygroups), validClassLabs)
    trainBgLabs <- setdiff(bggroups, validClassLabs)
    allTrainLabs <- union(trainClassLabs, trainBgLabs)
  
    if(length(nvals) == 0 & length(trainClassLabs) > 10){
      nmax <- min(3000, length(trainClassLabs))
      nvals <- round(exp(seq(log(10), log(nmax), (log(nmax) - log(10))/(samplePts-1)))) 
    }

    validmat <- mat1[which(classes %in% validClassLabs),]
    validclasses <- classes[which(classes %in% validClassLabs)]
  
    # Add baseline for cosine:
    print("Computing group sim")
    trainGrpSim <- innerProductGroups("cosine", mat1[which(classes %in% allTrainLabs),], classes[which(classes %in% allTrainLabs)], compact=1)
    validGrpSim <- innerProductGroups("cosine", validmat, validclasses, compact=1)
    
    print("Initial group sim computed")

    res <- rbind(res, data.frame(trainClasses=0, 
                                 rep=ii, 
                                 trainLoss=(mean(trainGrpSim$diff) - mean(unlist(trainGrpSim$same)))/sd(trainGrpSim$diff),
                                 trainAvgLoss=mean((mean(trainGrpSim$diff) - sapply(trainGrpSim$same, mean))/sd(trainGrpSim$diff)),
                                 validLoss=(mean(validGrpSim$diff) - mean(unlist(validGrpSim$same)))/sd(validGrpSim$diff),
                                 validAvgLoss=mean((mean(validGrpSim$diff) - sapply(validGrpSim$same, mean))/sd(validGrpSim$diff)),
                                 validAURank=1-mean(rankVectors(unlist(validGrpSim$same), validGrpSim$diff))))
    
    for (nclasses in nvals){
      print(sprintf("nclasses=%d", nclasses))
      myclasses <- union(sample(trainClassLabs, nclasses), trainBgLabs)
  
      trainmat <- mat1[which(classes %in% myclasses),]
      trainclasses <- classes[which(classes %in% myclasses)]

      trainmodel <- learnInnerProduct(trainmat, trainclasses, epochs=epochs)
      
      trainGrpSim <- innerProductGroups(trainmodel$model, trainmat, trainclasses, compact=1)
      validGrpSim <- innerProductGroups(trainmodel$model, validmat, validclasses, compact=1)
      
      res <- rbind(res, data.frame(trainClasses=nclasses, 
                                   rep=ii, 
                                   trainLoss=(mean(trainGrpSim$diff) - mean(unlist(trainGrpSim$same)))/sd(trainGrpSim$diff),
                                   trainAvgLoss=mean((mean(trainGrpSim$diff) - sapply(trainGrpSim$same, mean))/sd(trainGrpSim$diff)),
                                   validLoss=(mean(validGrpSim$diff) - mean(unlist(validGrpSim$same)))/sd(validGrpSim$diff),
                                   validAvgLoss=mean((mean(validGrpSim$diff) - sapply(validGrpSim$same, mean))/sd(validGrpSim$diff)),
                                   validAURank=1-mean(rankVectors(unlist(validGrpSim$same), validGrpSim$diff))))
    }
    
    validClassLabs <- c()
  }
  
  return(res)
}

#' rankVectors
#' 
#' rankVectors - given two vectors a and b, find for each element k of a: mean(b > k)
#' That is, find the percent rank within b of each element of a.  O(|a|+|b|) time.
#' This handles ties by putting all elements of b after all elements of a. 
#' @export
rankVectors <- function(a,b){
  if (length(b) == 0){
    return(numeric(length(a)))
  }
  
  # Given two input vectors a and b, returns for each element a_i of a: what fraction of elements of b are greater than or equal to a_i
  myvals <- data.frame(val=sort(a), ix=order(a))
  myvals <- rbind(myvals, data.frame(val=sort(b), ix=0))
  myvals <- cbind(myvals, randSeed=rnorm(dim(myvals)[1]))
  myvals <- myvals[order(myvals$val, myvals$randSeed), ]
  myvals$count <- cumsum(myvals$ix == 0)
  
  arank <- numeric(length(a))
  arank[myvals$ix[myvals$ix > 0]] <- 1 - myvals$count[myvals$ix > 0]/length(b)
  return(arank)
}

#' listify
#' listify - takes a vector and a list and maps the vector into the same structure as the list
#' equivalent to relisting and unlisted list.
listify <- function(mylist, myvec){
  mylengths <- sapply(mylist, length)
  mynames <- names(mylist)
  
  newlist <- vector("list", length(mylist))
  # Need to vectorize, but for the lists of interest, this is not costly
  ix <- 1
  for (ii in seq_along(mylist)){
    if (mylengths[ii] > 0){
      newlist[[ii]] <- myvec[ix:(ix + mylengths[ii]-1)]
      ix <- ix + mylengths[ii]
    } 
  }
  names(newlist) <- mynames
  return(newlist)
}


#' matLength
#' matLength - takes a matrix and computes the Euclidean length of each row
matLength <- function(mat1){
  return(sqrt(sapply(seq_len(dim(mat1)[1]), FUN=function(j) sum(mat1[j,]^2))))
}
