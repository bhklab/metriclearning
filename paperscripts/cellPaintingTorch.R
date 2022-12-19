library(cmapR)
library(CMAPToolkit)
library(ggplot2)
library(parallel)
library(torch)
library(pROC)

topdir <- "~/Work/bhk/data/cellpainting/bray-2017"
outdir <- "~/Work/bhk/analysis/metric_learning/2022_cellpaint/"

combds <- readRDS(file.path(topdir, "bray_2017_combined_profiles.rds"))

### Split data 
### Unnormalized data
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


### Normalized data
combMeta <- combds[, metax]
combData <- as.matrix(combds[, setdiff(seq(dim(combds)[2]), metax)])

combDataNorm <- scale(combData, center=TRUE, scale=TRUE)
# apparently some axes do not vary
combDataNorm[is.na(combDataNorm)] <- 0

norm_train <- genDataset(combDataNorm[ixTrain, 1:100], combMeta$Metadata_pert_id[ixTrain], 
                         pca_first = FALSE, scale = FALSE, center = FALSE)
norm_test <- genDataset(combDataNorm[ixTest, 1:100], combMeta$Metadata_pert_id[ixTest], 
                        pca_first = FALSE, scale = FALSE, center = FALSE)



train_norm <- dataloader(norm_train, batch_size = 1, shuffle = TRUE)
test_norm <- dataloader(norm_test, batch_size = 1, shuffle = TRUE)


model <- OneLayerLinear(100, 100)


device <- if (cuda_is_available()) torch_device("cuda:0") else "cpu"
model <- model$to(device = device)

optimizer <- optim_adam(model$parameters, lr = 0.01)

resnorm <- train_function(model=model, 
                      train_dl=train_norm, 
                      valid_dl=test_norm, 
                      myloss=mycos_t_loss, 
                      device=device, 
                      optimizer=optimizer, 
                      epochs = 10, 
                      #save_pars = "w"
)

mysimsnorm <- innerProductGroups(resnorm$model, combDataNorm[ixTest2, 1:100], combMeta$Metadata_pert_id[ixTest2])

ixTest2 <- which(combMeta$Metadata_pert_id %in% pertTest[1:500])
cos100 <- getCosineDist(mydim=100, N=1000, mkfig=0)

pdf(file=file.path(outdir, "bray", "bray_renormalized_simdensities_d=100.pdf"), width=10, height=8)
plot(density(cos100[upper.tri(cos100)], bw=0.01), col="black", lty=2, xlab="Modified cosine", ylab="Density", xlim=c(-1,1), main="Bray Learned Cosine Renormalized similarities, d=100 dims")
lines(density(unlist(mysimsnorm$same), bw=0.01), col="orange")
lines(density(mysimsnorm$diff, bw=0.01), col="darkgreen")
legend(x="topleft", legend=c("Same pert_id", "Different", "Null 100-dim"), col=c("orange", "darkgreen", "black"), lwd=c(2,2,2), lty=c(1,1,2))
dev.off()

# Unnormalized data for comparison
mysims <- innerProductGroups(res$model, combData[ixTest2, 1:100], combMeta$Metadata_pert_id[ixTest2])

pdf(file=file.path(outdir, "bray", "bray_nonorm_simdensities_d=100.pdf"), width=10, height=8)
plot(density(unlist(mysims$same), bw=0.01), col="red", lwd=2, xlab="Modified cosine", ylab="Density", xlim=c(-1,1), main="Bray Learned Cosine Unnormed similarities, d=100 dims")
lines(density(cos100[upper.tri(cos100)], bw=0.01), col="black", lty=2, lwd=2)
lines(density(unlist(mysims$diff), bw=0.01), col="blue", lwd=2)
legend(x="topleft", legend=c("Same pert_id", "Different", "Null 100-dim"), col=c("red", "blue", "black"), lwd=c(2,2,2), lty=c(1,1,2))
dev.off()


### Test cross validation refactoring

lip <- learnInnerProduct(combDataNorm[ixTrain2, 1:100], combMeta$Metadata_pert_id[ixTrain2], epochs=100)

ipgrps <- innerProductGroups(lip$model, combDataNorm[ixTrain2, 1:100], combMeta$Metadata_pert_id[ixTrain2])
ipgrps_comp <- innerProductGroups(lip$model, combDataNorm[ixTrain2, 1:100], combMeta$Metadata_pert_id[ixTrain2], compact = 1)

plot(density(ipgrps$diff, bw=0.01), col="blue", xlim=c(-1,1))
lines(density(unlist(ipgrps$same), bw=0.01), col="red")

plot(density(ipgrps_comp$diff, bw=0.01), col="blue", xlim=c(-1,1))
lines(density(unlist(ipgrps_comp$same), bw=0.01), col="red")

ipgrps_test <- innerProductGroups(lip$model, combDataNorm[ixTest2, 1:100], combMeta$Metadata_pert_id[ixTest2])
ipgrps_test2 <- innerProductGroups(lip$model, combDataNorm[ixTest, 1:100], combMeta$Metadata_pert_id[ixTest], compact=1)

cosgrps_test <- innerProductGroups(model="cosine", combDataNorm[ixTest2, 1:100], combMeta$Metadata_pert_id[ixTest2])


plot(density(ipgrps_test$diff, bw=0.01), col="blue", xlim=c(-1,1))
lines(density(unlist(ipgrps_test$same), bw=0.01), col="red")

plot(density(ipgrps_test2$diff, bw=0.01), col="blue", xlim=c(-1,1))
lines(density(unlist(ipgrps_test2$same), bw=0.01), col="red")



mxval <- metricCrossValidate(combDataNorm[ixTrain2, 1:100], combMeta$Metadata_pert_id[ixTrain2], nFolds=5, epochs=3) 




### Run on LINCS data
lincsdatadir <- "/Users/iansmith/Work/code/lincs-cell-painting/spherized_profiles/profiles"

ds <- loadLincsData("2016_04_01_a549_48hr_batch1_dmso_spherized_profiles_with_input_normalized_by_dmso.rds")

tres <- analyzeNormCPData(cppath="2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_whole_plate.rds",
                          outpath=".", method="xval", epochs=1, dsname="lincsB2Plate")