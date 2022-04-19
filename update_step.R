
update_step <- function(mat, classes, wts, m, v, tparam, b1=0.9, b2=0.999){
  
  eta <- 0.01
  mt <- numeric(length(m))
  vt <- numeric(length(v))
  wt <- numeric(length(wts))
  
  wts[wts < 0] <- 0
  
  # First, compute the gradient
  a <- CMAPToolkit::cosine(sqrt(wts) * mat, sqrt(wts)*mat)
  
  for (ix in seq_along(wts)){
    delta <- 0.01
    wts2 <- wts
    wts2[ix] <- wts[ix] + delta
    
    a2 <- CMAPToolkit::cosine((sqrt(wts2) * mat), sqrt(wts2) * mat)
    
    # Evaluate loss
    mygrad <- (calc_loss(a2, classes) - calc_loss(a, classes))/delta
    
    mt[ix] <- b1*m[ix] + (1-b1)*mygrad
    vt[ix] <- b2*v[ix] + (1-b2)*mygrad^2
    
    mcorr <- mt[ix]/(1 - b1^tparam)
    vcorr <- vt[ix]/(1 - b2^tparam)
    
    wt[ix] <- wts[ix] + eta * mcorr / (sqrt(vcorr) + 1e-8)
  }
  
  return(list(w=wt, m=mt, v=vt))
}


update_step_par <- function(mat, classes, wts, m, v, tparam, b1=0.9, b2=0.999){
  
  eta <- 0.01
  mt <- numeric(length(m))
  vt <- numeric(length(v))
  wt <- numeric(length(wts))
  
  wts[wts < 0] <- 0
  
  # First, compute the gradient
  a <- CMAPToolkit::cosine(sqrt(wts) * mat, sqrt(wts)*mat)
  myloss <- calc_loss(a, classes)
  
  mygrad <- unlist(mclapply(seq(length(wts)), FUN=function(x) inner_update(mat, classes, wts, myloss, ix=x), mc.cores=detectCores()))
  
  mt <- b1*m + (1-b1)*mygrad
  vt <- b2*v + (1-b2)*mygrad^2
  
  mcorr <- mt/(1 - b1^tparam)
  vcorr <- vt/(1 - b2^tparam)
  
  wt <- wts + eta * mcorr / (sqrt(vcorr) + 1e-8)
  
  return(list(w=wt, m=mt, v=vt))
}

# This is clumsy, and probably calculating the gradient for each ix is all we need
inner_update <- function(mat, classes, wts, myloss, ix){
  delta <- 0.01
  wts2 <- wts
  wts2[ix] <- wts[ix] + delta
  
  a2 <- CMAPToolkit::cosine((sqrt(wts2) * mat), sqrt(wts2) * mat)
  
  mygrad <- (calc_loss(a2, classes) - myloss)/delta
  return(mygrad)
}


calc_loss <- function(simmat, classes){
  
  ix <- sapply(sort(unique(classes)), FUN=function(x) which(classes==x))
  diag(simmat) <- NA
  
  same_sims <- unlist(sapply(ix[2:length(ix)], FUN=function(x) simmat[x,x]))
  other_sims <- sum(simmat, na.rm=TRUE) - sum(same_sims, na.rm=TRUE)
  
  # This is gross, but I need all the similarities from different classes to compute the sd 
  enum_other_sims <- c(simmat[ix[[1]], ix[[1]]], 
                       unlist(sapply(seq_along(ix), FUN=function(x) 
                         sapply(setdiff(seq_along(ix), x), FUN=function(y) simmat[ix[[x]], ix[[y]]]))))
  
  possim <- sum(same_sims, na.rm=TRUE)/sum(!is.na(same_sims))
  negsim <- other_sims/(sum(!is.na(simmat)) - sum(!is.na(same_sims)))
  
  myloss <- (possim - negsim)/sd(enum_other_sims, na.rm=1)
  return(myloss)
}


calc_simmat <- function(ds1, ds2, rotmat, wts){
  if (is.numeric(ds1)){
    # Assume ds1 is a matrix
    xmat1 <- t(rotmat) %*% ds1
    xmat2 <- t(rotmat) %*% ds2
    return (CMAPToolkit::cosine(sqrt(wts) * xmat1, sqrt(wts) * xmat2))
  } else {
    xmat1 <- t(rotmat) %*% ds1@mat
    xmat2 <- t(rotmat) %*% ds2@mat
    return (CMAPToolkit::cosine(sqrt(wts) * xmat1, sqrt(wts) * xmat2))
  }
}


embedspace <- function(ds, rotmat, wts){
  if (is.numeric(ds)){
    return(sqrt(wts) * (t(rotmat) %*% ds))
  } else {
    return(sqrt(wts) * (t(rotmat) %*% ds@mat))
  }
}