library(cmapR)
library(CMAPToolkit)
library(ggplot2)
library(parallel)
library(coop)

library(torch)

source("linear_embedding_torch.R")


torch_manual_seed(42)




hepg2_train <- l1k_dataset(metapath=metapath, datapath=datapath)
vcap_tune <- l1k_dataset(metapath=metapath, datapath=datapath, cell_id="VCAP")

train_dl <- dataloader(hepg2_train, batch_size = 1, shuffle = TRUE)
tune_dl <- dataloader(vcap_tune, batch_size = 1, shuffle = TRUE)



model <- OneLayerLinear(978)


device <- if (cuda_is_available()) torch_device("cuda:0") else "cpu"
model <- model$to(device = device)


optimizer <- optim_adam(model$parameters, lr = 0.01)


res <- train_function(model=model, 
					train_dl=train_dl, 
					valid_dl=valid_dl, 
					myloss=mycos_t_loss, 
					device=device, 
					optimizer=optimizer)


## pairwise_sim = torch_cosine_similarity(x1[:, None, :], x2[None, :, :], dim=-1)

## lets plot 


output <- model(hepg2_train$sigs$to(device=device))$to(device="cpu")

embedded_pca <-prcomp(output)

cols <- rainbow(hepg2_train$.length())



pairs(embedded_pca$x[,1:6], col=cols[as.numeric(factor(hepg2_train$mysigs$pert_iname))])


sims <- coop::cosine(t(as.array(output)))


same_cors <- numeric(0)
for(drug in hepg2_train$mycmpds){


	myx <- which(hepg2_train$mysigs$pert_iname==drug)
	same_cors <- c(same_cors,sims[myx,myx])


}



plot(density(sample(sims,10000)), col="blue", xlim=c(-1,1))
lines(density(same_cors),col='red')



output_valid <- model(vcap_tune$sigs)




sims_valid <- coop::cosine(t(as.array(output_valid)))



same_cors_valid <- numeric(0)
for(drug in vcap_tune$mycmpds){


	myx <- which(vcap_tune$mysigs$pert_iname==drug)
	same_cors_valid <- c(same_cors_valid,sims_valid[myx,myx])


}



plot(density(sample(sims,10000)),col='blue')
lines(density(same_cors_valid), col="red")



plot(density(same_cors), col="green")
lines(density(same_cors_valid),col='purple')







hepg2_train <- l1k_dataset(metapath=metapath, datapath=datapath, pca_first=TRUE)
vcap_tune <- l1k_dataset(metapath=metapath, datapath=datapath, cell_id="VCAP", pca_first=TRUE)

train_dl <- dataloader(hepg2_train, batch_size = 1, shuffle = TRUE)
tune_dl <- dataloader(vcap_tune, batch_size = 1, shuffle = TRUE)

model <- DiagonalOnly(978)



device <- if (cuda_is_available()) torch_device("cuda:0") else "cpu"
model <- model$to(device = device)


optimizer <- optim_adam(model$parameters, lr = 0.01)


res <- train_function(epochs=5,model=model, 
					  train_dl=train_dl, 
					  valid_dl=valid_dl, 
					  myloss=mycos_t_loss,
					   device=device, optimizer=optimizer, save_pars="w")




output <- model(hepg2_train$sigs$to(device=device))$to(device="cpu")

sims <- coop::cosine(t(as.array(output)))


same_cors <- numeric(0)
for(drug in hepg2_train$mycmpds){


	myx <- which(hepg2_train$mysigs$pert_iname==drug)
	same_cors <- c(same_cors,sims[myx,myx][upper.tri(sims[myx,myx])])


}



plot(density(sample(sims,10000)), col="blue", xlim=c(-1,1))
lines(density(same_cors),col='red')
