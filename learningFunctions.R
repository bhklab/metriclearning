#' metricCrossValidate 
#' takes a matrix and class vector, runs n-fold cross validation with a learned inner product,
#' and reports the test and training loss with training.
#' @param mat1 Numeric, a matrix of size N x d, where rows are samples and columns are features
#' @param classes Factor (or character, numeric), where two elements with the same value are assumed to be similar
#' @param nFolds Numeric, number of folds to run (default = 5)
#' @param metric Character, denotes which type of metric to learn. One of {}.
#' @param epochs Numeric, the number of epochs of training to run (default = 100)

metricCrossValidate <- function(mat1, classes, nFolds=5, metric="", epochs=100){
# Revise to intelligently choose folds in the case of a wide range of class sizes
# But for now, do it naively
  
  mygroups <- table(classes)
  myfolds <- createFolds(names(mygroups), k=nFolds)
  #myfolds0 <- createFolds(names(mygroups)[mygroups == 1], k=nFolds)
  
  res_all <- list(myfolds=myfolds)
  
  for (testset in myfolds){
    testix <- which(classes %in% testset)
    trainix <- setdiff(seq_along(classes), testix)
    
    res <- learnInnerProduct(mat1[trainix, ], classes[trainix], metric=metric, epochs=epochs)  
    
    
    res_all[sprintf("fold%d", length(res_all))] <- list(res)
  }
  
  return(res_all)
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
  return(list(model=model, res=res))
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
  
  if (model == "cosine"){
    return(cosine(mat1, mat2))
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
  
  mag1 <- sqrt(diag(mat1 %*% t(mat1)))
  mag2 <- sqrt(diag(mat2 %*% t(mat2)))
  
  return(a/outer(mag1, mag2, "*"))
}



innerProductGroups <- function(model, mat1, groupings, compact=0){
  if (compact){
    # Use compact if the matrix is too large to effectively compute
    # Compact samples from the space of unlike similarities rather than computing the entire matrix
    
    # This could likely be accomplished more elegantly with dplyr
    samesim <- sapply(names(table(groupings)[table(groupings)>1]), FUN=function(x){
      ix <- which(groupings == x)
      a <- innerProduct(model, mat1[ix, ], mat1[ix, ])
      a[upper.tri(a)]
    })
    
    jx <- sample(dim(mat1)[1], 1000)
    diffsim <- innerProduct(model, mat1[jx, ], mat1[jx, ])
    for (mygrp in names(table(groupings)[table(groupings) > 1])){
      ix <- which(jx %in% which(groupings == mygrp))
      if (length(ix) > 1){
        diffsim[ix, ix] <- NA
      }
    }
    diffsim <- diffsim[upper.tri(diffsim)]
    diffsim <- diffsim[!is.na(diffsim)]
    return(list(same=samesim, diff=diffsim))
    
  } else {
    sims <- innerProduct(model, mat1, mat1)

    samesim <- sapply(names(table(groupings)[table(groupings) > 1]), FUN=function(x) {
      a <- sims[which(groupings == x), which(groupings == x)]
      a[upper.tri(a)]
      })
    
    diffsim <- sims
    for (mygroup in groupings){
      ix <- which(groupings == mygroup)
      diffsim[ix,ix] <- NA
    }
    
    diffsim <- diffsim[!is.na(diffsim)]
    return(list(same=samesim, diff=diffsim))
  }
}