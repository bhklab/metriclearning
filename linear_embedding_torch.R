library(cmapR)
library(CMAPToolkit)
library(ggplot2)
library(parallel)
library(coop)
library(torch)

l1k_dataset <- dataset(
  
  name = "L1K Sigs",
  
  initialize = function(datapath, metapath, cell_id="HEPG2", pca_first=FALSE, ndim=978) {
    
    l1k_meta <- read_l1k_meta(metapath, version=2020)
  	attach(l1k_meta)
  
  	self$mysigs <- siginfo[siginfo$cell_id == cell_id & siginfo$pert_type == "trt_cp", ]
  
  	ds <- parse_gctx(get_level5_ds(datapath), rid=l1k_meta$landmarks$pr_gene_id, cid=self$mysigs$sig_id)
    
    if(pca_first){
    	 pcspace <- prcomp(t(ds@mat), scale=FALSE, center=FALSE)
    	 self$sigs <- torch_tensor(pcspace$x)
    } else {
    	self$sigs <- torch_tensor(t(ds@mat))
    }

  	cpcount <- table(self$mysigs$pert_iname)
  	self$mycmpds <- names(cpcount)[which(cpcount >= 12)]
  	self$ndim <- ndim
  },
  
  .getitem = function(i) {
  	mycp <- self$mycmpds[[i]]
      
    ix <- which(self$mysigs$pert_iname == mycp)
    jx <- sample(setdiff(seq_len(dim(self$sigs)[1]), ix), 2*length(ix))

    mymat <- torch_cat(list(self$sigs[ix, 1:self$ndim], self$sigs[jx, 1:self$ndim]))
    myclass <- torch_tensor(c(rep(1, length(ix)), rep(0, length(jx))))


    list(x = mymat, y = myclass)
    
  },
  
  .length = function() {
    length(self$mycmpds)
  }
 
)


genDataset <- dataset(
  name = "Generic DS",
  
  initialize = function(mymat, identVec, pca_first, scale=FALSE, center=FALSE) {
    # mymat has signatures as rows, features as columns
    if(pca_first){
      pcspace <- prcomp(mymat, scale=scale, center=center)
      self$sigs <- torch_tensor(pcspace$x)
    } else {
      self$sigs <- torch_tensor(mymat)
    }
    
    cpcount <- table(identVec)
    self$mycmpds <- names(cpcount)[which(cpcount >= 3)]
    self$ndim <- dim(mymat)[2]
    self$identVec <- identVec
  },
  
  .getitem = function(i) {
    mycp <- self$mycmpds[[i]]
    
    ix <- which(self$identVec == mycp)
    jx <- sample(setdiff(seq_len(dim(self$sigs)[1]), ix), 2*length(ix))
    
    mymat <- torch_cat(list(self$sigs[ix, 1:self$ndim], self$sigs[jx, 1:self$ndim]))
    myclass <- torch_tensor(c(rep(1, length(ix)), rep(0, length(jx))))
    
    list(x = mymat, y = myclass)
  },
  
  .length = function() {
    length(self$mycmpds)
  }
)


torch_manual_seed(42)

OneLayerLinear <- nn_module(
  
  "linear_embedding",
  initialize = function(embedding_dim, nfeats) {
    self$fc1 <- nn_linear(in_features = nfeats, out_features = embedding_dim, bias=FALSE)
  },
  
  forward = function(x) {
    x %>% self$fc1() 
  }
)



DiagonalOnly <- nn_module(
  
  "Just scaling PCs",
  initialize = function(num_PCs) {
    self$w <- nn_parameter(torch_ones(num_PCs, requires_grad=TRUE))
  },
  
  forward = function(x) {
    torch_matmul(x,torch_diag(torch_sqrt(self$w))) 
  }
)




mycos_sim_loss <- function(x, y){

	sim_mat <- torch_triu(torch_cosine_similarity(x$reshape(c(-1, 1, dim(x)[2])), 
    								  x$reshape(c(1, -1, dim(x))), dim=3), diagonal=1)

	target <- torch_outer(y,y)

	torch_sum(target*(1-sim_mat) + (1-target)*(sim_mat))/((y$shape-1)*y$shape/2)
	
}



mycos_t_loss <- function(x, y){

  sim_mat <- torch_cosine_similarity(x$reshape(c(-1, 1, dim(x)[2])), 
                      x$reshape(c(1, -1, dim(x)[2])), dim=3)

  target <- torch_outer(y,y)

  target <- target + -(torch_max(y)+1)*torch_eye(target$shape[1]) ## set diagonal to negative numbers.

  same_mean <- torch_mean(torch_masked_select(sim_mat, target==1))
  diff_mean <- torch_mean(torch_masked_select(sim_mat, target==0))

  diff_std <- torch_std(torch_masked_select(sim_mat, target==0))

  -(same_mean-diff_mean)/(diff_std)
  
}






train_function <- function(model, epochs=5, train_dl, valid_dl, myloss, device, optimizer, save_pars=NA){
  par_list = list()
  
  mean_train_losses <- c()
  mean_valid_losses <- c()
  
  for (epoch in seq_len(epochs)) {

    model$train()
    #train_losses <- c()  

    coro::loop(for (b in train_dl) {
      optimizer$zero_grad()

      output <- model(b$x$to(device = device))

      loss <- myloss(output$reshape(c(-1, output$shape[3])),b$y$to(device = device)$reshape(c(-1)))
      if(is.nan(loss$item())){
        browser()
      }
      loss$backward()
      optimizer$step()
      #train_losses <- c(train_losses, loss$item())
      # loss <- nnf_binary_cross_entropy(output, b$y$to(dtype = torch_float(), device = device))
      # loss$backward()
      # optimizer$step()
      # train_losses <- c(train_losses, loss$item())
    })
    if(!is.na(save_pars)){
      # What?!  I have no idea what this does. 
      par_list[[epoch]] <- lapply(save_pars, function(par) return(as.array(model[[par]]$cpu())))
    }
    model$eval()
    train_losses <- c()
    valid_losses <- c()

    coro::loop(for (b in valid_dl) {
      output <- model(b$x$to(device = device))

      loss <- myloss(output$reshape(c(-1, output$shape[3])),b$y$to(device = device)$reshape(c(-1)))
      valid_losses <- c(valid_losses, loss$item())
    })
    coro::loop(for (b in train_dl){
      output <- model(b$x$to(device = device))
      
      loss <- myloss(output$reshape(c(-1, output$shape[3])), b$y$to(device = device)$reshape(c(-1)))
      train_losses <- c(train_losses, loss$item())
    })

    mean_train_losses <- c(mean_train_losses, mean(train_losses))
    mean_valid_losses <- c(mean_valid_losses, mean(valid_losses))
    cat(sprintf("Loss at epoch %d: training: %3f, validation: %3f\n", epoch, mean(train_losses), mean(valid_losses)))
  }
  return(list(model=model, mean_train_losses=mean_train_losses, mean_valid_losses=mean_valid_losses, 
              train_losses=train_losses, valid_losses=valid_losses,par_list=par_list))
}


