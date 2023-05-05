library(CMAPToolkit)

#' innerProductPairwiseGroups
#' 
#' This function takes a metric argument, a data matrix, group labels, and sets of group labels as input and 
#' computes the similarity between signatures (rows) in the matrix that belong to the same set of group
#' labels but have different group labels. The use case is that there are sets of groups that are similar (e.g. 
#' compounds with the same mechanism of action), and we seek the similarities of elements of the same set that
#' have different group labels. 

simPairwiseGroups <- function(metric, mat1, classes, sets, compact=1){
  # Compact samples a block from the full similarity matrix, which can be efficient if mat1 is very large.
  # Note that compact does not exclude same elements, i.e. it is null distribution that includes both 
  # similarities where the null hypothesis and where the alternate hypothesis is true. The presumption is that 
  # the rate of instances of the alternate hypothesis is sufficiently low that this isn't a problem. 
  if (compact){
    if (dim(mat1)[1] < 2000){
      sims <- calcSimBlock(ds1=t(mat1), ds2=t(mat1), metric=metric, parallel=TRUE)
    } else {
      ix <- sample(dim(mat1)[1], 2000)
      sims <- calcSimBlock(ds1=t(mat1[ix,]), ds2=t(mat1[ix,]), metric=metric, parallel=TRUE)
    }
    
    # Consider filtering same-groups and same-sets elements here
    
  } else {
    sims <- calcSimBlock(ds1=t(mat1), ds2=t(mat1), metric=metric, parallel=TRUE)
    
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
      
      mysim <- calcSimBlock(ds1=t(mat1[jx,]), ds2=t(mat1[jx,]), metric=metric, parallel=TRUE) 
      setSims[[ii]] <- (mysim[upper.tri(mysim)])[filt]
    } 
  }
  
  if (!is.null(names(sets))){
    names(setSims) <- names(sets)
  }
  
  return(list(setSims=setSims, allSims=sims))
}