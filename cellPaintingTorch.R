library(cmapR)
library(CMAPToolkit)
library(ggplot2)
library(parallel)
library(torch)

topdir <- "~/Work/bhk/data/cellpainting/bray-2017"
outdir <- "~/Work/bhk/analysis/metric_learning/2022_cellpaint/"

combds <- readRDS(file.path(topdir, "bray_2017_combined_profiles.rds"))

# Split data 

metax <- grep("Meta", colnames(combds))
combMeta <- combds[, metax]
combData <- as.matrix(combds[, setdiff(seq(dim(combds)[2]), metax)])

pertCounts <- table(combMeta$Metadata_pert_id)
pertTrain <- sample(names(pertCounts), 24000)
pertTest <- setdiff(names(pertCounts), pertTrain)
ixTrain <- which(combMeta$Metadata_pert_id %in% pertTrain)
ixTest <- which(combMeta$Metadata_pert_id %in% pertTest)
  
bray_train <- genDataset(combData[ixTrain, 1:100], combMeta$Metadata_pert_id[ixTrain], 
                         pca_first = FALSE, scale = FALSE, center = FALSE)
bray_test <- genDataset(combData[ixTest, 1:100], combMeta$Metadata_pert_id[ixTest], 
                        pca_first = FALSE, scale = FALSE, center = FALSE)



train_dl <- dataloader(bray_train, batch_size = 1, shuffle = TRUE)
test_dl <- dataloader(bray_test, batch_size = 1, shuffle = TRUE)


model <- OneLayerLinear(100, 100)


device <- if (cuda_is_available()) torch_device("cuda:0") else "cpu"
model <- model$to(device = device)

optimizer <- optim_adam(model$parameters, lr = 0.01)

res <- train_function(model=model, 
                      train_dl=train_dl, 
                      valid_dl=test_dl, 
                      myloss=mycos_t_loss, 
                      device=device, 
                      optimizer=optimizer, 
                      epochs = 3, 
                      #save_pars = "w"
)