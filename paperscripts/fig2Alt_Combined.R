source("figinit.R")
source("paperFigFuncs.R")
source("../learningFunctions.R")


#### Figure 2: L1000 performance ####

l1kdir <- file.path(topdir, "modelRuns")

### Load the L1000 data:
mycells <- sapply(strsplit(sapply(strsplit(list.files(l1kdir, pattern="L1Kxval"), "cell="), 
                                  FUN=function(x) x[[2]]), ".rds"), FUN=function(x) x[[1]])

# !!!!!! Fix this once you have the all_xval models
mycells <- setdiff(mycells, c("all", "MCF10A", "SKL"))

if (!file.exists(file.path(outdir, "../figspaper_res", "L1K_xvalres.rds"))){
  xvalres <- list()
  
  for (ii in seq_along(mycells)){
    mycell <- mycells[ii]
    print(sprintf("%d: %s", ii, mycell))
    
    xvalds <- readRDS(file.path(l1kdir, list.files(l1kdir, pattern=sprintf("L1Kxval.*%s.rds", mycell))))
    f <- list.files(file.path(l1kdir, "models"), pattern=sprintf("L1Kxval.*%s", mycell))
    
    xvalres[[mycell]] <- getL1KXValReps(myfolds=xvalds$myfolds, 
                                        models=file.path(l1kdir, "models", f), 
                                        datapath=datapath,
                                        l1kmeta=l1kmeta, 
                                        cellid=mycell)
    
  }
  
  saveRDS(xvalres, file=file.path(outdir, "../figspaper_res/L1K_xvalres.rds"))
} else {
  xvalres <- readRDS(file.path(outdir, "../figspaper_res/L1K_xvalres.rds"))
}

if (!file.exists(file.path(outdir, "../figspaper_res/L1K_PertRes.rds"))){
  L1KPertBalAUCML <- 1 - sapply(mycells, FUN=function(x) 
    sapply(xvalres[[x]]$mlRanks, FUN=function(y) 
      mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(y)))))))
  
  L1KPertBalAUCCos <- 1 - sapply(mycells, FUN=function(x) 
    sapply(xvalres[[x]]$cosRanks, FUN=function(y) 
      mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(y)))))))
  
  ### mean compound AUC - avoids sampling
  L1KPertMeanAUCML <- 1 - sapply(mycells, FUN=function(x) 
    mean(sapply(xvalres[[x]]$mlRanks, FUN=function(y)
      mean(sapply(y, FUN=mean)))))
  
  L1KPertMeanAUCCos <- 1 - sapply(mycells, FUN=function(x) 
    mean(sapply(xvalres[[x]]$cosRanks, FUN=function(y)
      mean(sapply(y, FUN=mean)))))
  
  # FDR fractions - used balanced 
  L1KPertBalFDRML <- list(fdr01=sapply(mycells, FUN=function(x) 
    sapply(xvalres[[x]]$mlRanks, FUN=function(y) 
      mean(sapply(seq(100), FUN=function(z) mean(p.adjust(unlist(balancedSample(y)), "fdr") < 0.01))))), 
    fdr05=sapply(mycells, FUN=function(x) 
      sapply(xvalres[[x]]$mlRanks, FUN=function(y) 
        mean(sapply(seq(100), FUN=function(z) mean(p.adjust(unlist(balancedSample(y)), "fdr") < 0.05))))),
    fdr10=sapply(mycells, FUN=function(x) 
      sapply(xvalres[[x]]$mlRanks, FUN=function(y) 
        mean(sapply(seq(100), FUN=function(z) mean(p.adjust(unlist(balancedSample(y)), "fdr") < 0.10)))))
  )
  
  L1KPertBalFDRCos <- list(fdr01=sapply(mycells, FUN=function(x) 
    sapply(xvalres[[x]]$cosRanks, FUN=function(y) 
      mean(sapply(seq(100), FUN=function(z) mean(p.adjust(unlist(balancedSample(y)), "fdr") < 0.01))))), 
    fdr05=sapply(mycells, FUN=function(x) 
      sapply(xvalres[[x]]$cosRanks, FUN=function(y) 
        mean(sapply(seq(100), FUN=function(z) mean(p.adjust(unlist(balancedSample(y)), "fdr") < 0.05))))),
    fdr10=sapply(mycells, FUN=function(x) 
      sapply(xvalres[[x]]$cosRanks, FUN=function(y) 
        mean(sapply(seq(100), FUN=function(z) mean(p.adjust(unlist(balancedSample(y)), "fdr") < 0.10)))))
  )
  
  L1KPertBalFDRCos$method <- "cosine"
  L1KPertBalFDRML$method <- "metric learning"
  L1KPertBalFDRCos$cell <- mycells
  L1KPertBalFDRML$cell <- mycells
  
  saveRDS(list(L1KPertBalAUCML=L1KPertBalAUCML, 
               L1KPertBalAUCCos=L1KPertBalAUCCos, 
               L1KPertMeanAUCML=L1KPertMeanAUCML,
               L1KPertMeanAUCCos=L1KPertMeanAUCCos,
               L1KPertBalFDRML=L1KPertBalFDRML, 
               L1KPertBalFDRCos=L1KPertBalFDRCos), 
          file=file.path(outdir, "../figspaper_res/L1K_PertRes.rds"))
} else {
  pertRes <- readRDS(file.path(outdir, "../figspaper_res/L1K_PertRes.rds"))
  attach(pertRes)
}

if (!file.exists(file.path(outdir, "../figspaper_res/L1K_PCLRes.rds"))){
  pclds <- readRDS("data/pclds.rds")
  
  pclcells <- sapply(strsplit(sapply(strsplit(list.files(l1kdir, pattern="L1Kmetric"), "cell="), 
                                     FUN=function(x) x[[2]]), ".rds"), FUN=function(x) x[[1]])
  pclcells <- setdiff(pclcells, "all")
  
  pclres <- list()
  
  for (ii in seq_along(pclcells)){
    acell <- pclcells[ii]
    print(acell)
    modelds <- readRDS(file.path(l1kdir, list.files(l1kdir, pattern=sprintf("L1Kmetric.*%s.rds", acell))))
    f <- list.files(file.path(l1kdir, "models"), pattern=sprintf("L1Kmetric.*%s", acell))
    
    pclres[[acell]] <- getL1KMoA(datapath, l1kmeta, acell, mymodel=file.path(l1kdir, "models", f[1]), pclds=pclds)
  }
  
  L1KMOABalAUCML <- 1 - sapply(pclcells, FUN=function(x) 
    mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(pclres[[x]]$mlRanks, 1000))))))
  
  L1KMOABalAUCCos <- 1 - sapply(pclcells, FUN=function(x) 
    mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(pclres[[x]]$cosRanks, 1000))))))
  
  # SNR 
  L1KMoASNRML <- sapply(pclres, FUN=function(x) (sapply(x$mlPCLs$setSims, mean) - mean(x$mlPCLs$allSims))/sd(x$mlPCLs$allSims))
  L1KMoASNRCos <- sapply(pclres, FUN=function(x) (sapply(x$cosPCLs$setSims, mean) - mean(x$cosPCLs$allSims))/sd(x$cosPCLs$allSims))
  L1KMoASize <- sapply(pclres, FUN=function(x) sapply(x$mlRanks, length))
  
  # MoA FDR
  L1KMoABalFDRML <- data.frame(fdr01=sapply(pclres, FUN=function(x) 
    mean(sapply(seq(10), FUN=function(y) mean(p.adjust(unlist(balancedSample(x$mlRanks, k=1000)), "fdr") < 0.01)))), 
    fdr05=sapply(pclres, FUN=function(x) 
      mean(sapply(seq(10), FUN=function(y) mean(p.adjust(unlist(balancedSample(x$mlRanks, k=1000)), "fdr") < 0.05)))), 
    fdr10=sapply(pclres, FUN=function(x) 
      mean(sapply(seq(10), FUN=function(y) mean(p.adjust(unlist(balancedSample(x$mlRanks, k=1000)), "fdr") < 0.10)))) 
  )
  
  L1KMoABalFDRCos <- data.frame(fdr01=sapply(pclres, FUN=function(x) 
    mean(sapply(seq(10), FUN=function(y) mean(p.adjust(unlist(balancedSample(x$cosRanks, k=1000)), "fdr") < 0.01)))), 
    fdr05=sapply(pclres, FUN=function(x) 
      mean(sapply(seq(10), FUN=function(y) mean(p.adjust(unlist(balancedSample(x$cosRanks, k=1000)), "fdr") < 0.05)))), 
    fdr10=sapply(pclres, FUN=function(x) 
      mean(sapply(seq(10), FUN=function(y) mean(p.adjust(unlist(balancedSample(x$cosRanks, k=1000)), "fdr") < 0.10)))) 
  )
  
  L1KMoABalFDRML$cell <- row.names(L1KMoABalFDRML)
  L1KMoABalFDRML$method <- "metric learning"
  L1KMoABalFDRCos$cell <- row.names(L1KMoABalFDRCos)
  L1KMoABalFDRCos$method <- "cosine"
  
  
  saveRDS(list(pclres=pclres, 
               L1KMOABalAUCML=L1KMOABalAUCML,
               L1KMOABalAUCCos=L1KMOABalAUCCos, 
               L1KMoASNRML=L1KMoASNRML, 
               L1KMoASNRCos=L1KMoASNRCos,
               L1KMoABalFDRML=L1KMoABalFDRML,
               L1KMoABalFDRCos=L1KMoABalFDRCos), 
          file=file.path(outdir, "../figspaper_res/L1K_PCLRes.rds"))
} else {
  L1K_PCLRes <- readRDS(file.path(outdir, "../figspaper_res/L1K_PCLRes.rds"))
  attach(L1K_PCLRes)
}


### Load the Cell Painting data (CDRP):

braydir <- file.path(cpdir, "bray/models")
pclds <- readRDS("data/pclds.rds")


# Load or generate results by running models:
if (!file.exists(file.path(outdir, "../figspaper_res/cpSummaryRes.rds"))){
  brayds <- loadBrayData(braypath)
  lincsds1 <- loadLincsData(file.path(lincspath, "2016_04_01_a549_48hr_batch1_dmso_spherized_profiles_with_input_normalized_by_dmso.rds"), byCell = 1, splitGrps = 0)
  lincsds2 <- loadLincsData(file.path(lincspath, "2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_dmso.rds"), byCell = 0, splitGrps = 0)
  
  Lds2_a549 <- splitCPds(lincsds2, which(lincsds2$metads$Metadata_cell_id == "A549"))
  Lds2_mcf7 <- splitCPds(lincsds2, which(lincsds2$metads$Metadata_cell_id == "MCF7"))
  Lds2_u2os <- splitCPds(lincsds2, which(lincsds2$metads$Metadata_cell_id == "U2OS"))
  
  lincsds2_hidose <- c()
  lincsds2_hidose$ds <- lincsds2$ds[lincsds2$metads$Metadata_mmoles_per_liter >= 10, ]
  lincsds2_hidose$metads <- lincsds2$metads[lincsds2$metads$Metadata_mmoles_per_liter >= 10, ]
  
  l1m <- table(lincsds1$metads$Metadata_moa)
  l2m <- table(lincsds2$metads$Metadata_moa)
  
  # The first element is "", i.e. empty and unassigned, so exclude this.
  # It is very important that lincsds2 be read with byCell = 0 and splitGrps = 0, otherwise cell IDs and groupIDs get appended to pert identifiers. 
  # The by Cell and splitGrps flags are important for training the metric and defining batch sizes, not for evaluation. 
  lincsMoADS <- list(moaNames = c(names(l1m)[-1], names(l2m)[-1]), 
                     pertids = c(sapply(names(l1m)[-1], FUN=function(x) unique(lincsds1$metads$Metadata_pert_id[lincsds1$metads$Metadata_moa == x])), 
                                 sapply(names(l2m)[-1], FUN=function(x) unique(lincsds2$metads$Metadata_pert_id[lincsds2$metads$Metadata_moa == x]))), 
                     pertnames = c(sapply(names(l1m)[-1], FUN=function(x) unique(lincsds1$metads$Metadata_pert_iname[lincsds1$metads$Metadata_moa == x])), 
                                   sapply(names(l2m)[-1], FUN=function(x) unique(lincsds2$metads$Metadata_pert_iname[lincsds2$metads$Metadata_moa == x]))))
  
  
  braymodel <- torch::torch_load(file.path(braydir, "braymetric_epch=10_smp=50_model.pt"))
  lincs1model <- torch::torch_load(file.path(cpdir, "lincs/models/lincsmetric_a549_2016_dmso_epch=10_model.pt"))
  lincs2model <- torch::torch_load(file.path(cpdir, "lincs/models/lincsmetric_batch2_dmso_bycell_epch=10_model.pt"))
  
  
  #### Compound Replicates
  # Take care to exclude those compounds that are not compound treatments. Note that the LINCS datasets still have substantial class imbalances.
  brayReps <- getCPReps(brayds$ds, mymodel=braymodel, pertlabels=brayds$metads$Metadata_pert_id, datalabel="Bray DS centered full")
  lincs1Reps <- getCPReps(lincsds1$ds, mymodel=lincs1model, pertlabels=lincsds1$metads$Metadata_pert_id, datalabel="LINCS A549 DMSO full model")
  lincs2Reps <- getCPReps(lincsds2$ds, mymodel=lincs2model, pertlabels=lincsds2$metads$Metadata_pert_id, datalabel="LINCS Batch2 DMSO full model")
  
  lincs2Reps_A549 <- getCPReps(Lds2_a549$ds, mymodel=lincs2model, pertlabels=Lds2_a549$metads$Metadata_pert_id, datalabel="LINCS Batch2 DMSO full model, A549")
  lincs2Reps_MCF7 <- getCPReps(Lds2_mcf7$ds, mymodel=lincs2model, pertlabels=Lds2_mcf7$metads$Metadata_pert_id, datalabel="LINCS Batch2 DMSO full model, MCF7")
  lincs2Reps_U2OS <- getCPReps(Lds2_u2os$ds, mymodel=lincs2model, pertlabels=Lds2_u2os$metads$Metadata_pert_id, datalabel="LINCS Batch2 DMSO full model, U2OS")
  
  
  #### Cross validation results
  brayXvalSets <- list(file.path(braydir, list.files(braydir, pattern="brayxval_epch=10_smp=20")), 
                       file.path(braydir, list.files(braydir, pattern="brayxval_epch=5_smp=50")))
  
  brayXvalLabs <- c("brayxval_10x20", "brayxval_5x50") #, ", "lincs2_10x5")
  
  cpxvalres <- list()
  for (ii in seq_along(brayXvalSets)){
    myxvalset <- brayXvalSets[[ii]]
    modlab <- brayXvalLabs[ii]
    xvalds <- readRDS(myxvalset[grep("rds", myxvalset)])
    cpxvalres[[ii]] <- getCPXvalReps(brayds$ds, brayds$metads, xvalds$myfolds, myxvalset[-length(myxvalset)], modlab)
  }
  
  lincs1XvalSets <- list(file.path(cpdir, "lincs/models", list.files(file.path(cpdir, "lincs/models"), pattern="lincsxval_a549_2016_dmso_epch=10_folds=5")))
  lincs1XvalLabs <- c("lincs1_10x5")
  
  for (ii in seq_along(lincs1XvalSets)){
    myxvalset <- lincs1XvalSets[[ii]]
    modlab <- lincs1XvalLabs[ii]
    xvalds <- readRDS(myxvalset[grep("rds", myxvalset)])
    cpxvalres[[length(cpxvalres)+1]] <- getCPXvalReps(lincsds1$ds, lincsds1$metads, xvalds$myfolds, myxvalset[-length(myxvalset)], modlab)
  }
  
  # Add LINCS2
  
  
  #### MoA Similarities
  # These MoA results use the L1000 PCLs (pclds) as the definitive labels for what constitutes a MoA class
  brayMoA <- getCPMoA(brayds$ds, brayds$metads, mymodel=braymodel, pertlabels=brayds$metads$Metadata_pert_id, 
                      pclds=pclds, datalabel="Bray MoA, centered full")
  lincs1MoA <- getCPMoA(lincsds1$ds, lincsds1$metads, mymodel=lincs1model, pertlabels=lincsds1$metads$Metadata_pert_id, 
                        pclds=pclds, datalabel="LINCS A549 MoA DMSO full model")
  lincs2MoA <- getCPMoA(lincsds2$ds, lincsds2$metads, mymodel=lincs2model, pertlabels=lincsds2$metads$Metadata_pert_id, 
                        pclds=pclds, datalabel="LINCS Batch2 MoA DMSO full model")
  lincs2MoAByCell <- getCPMoA(lincsds2$ds, lincsds2$metads, mymodel=lincs2model, pertlabels=lincsds2$metads$Metadata_pert_id, 
                              pclds=pclds, datalabel="LINCS Batch2 MoA DMSO full model")
  
  lincs2MoAByCellHighDose <- getCPMoA(lincsds2_hidose$ds, lincsds2_hidose$metads, mymodel=lincs2model, pertlabels=lincsds2_hidose$metads$Metadata_pert_id, 
                                      pclds=pclds, datalabel="LINCS Batch2 MoA DMSO full model High Dose")
  
  # The 'Native' MoA results use mechanism of action annotations directly from the LINCS Cell Painting datasets.
  brayMoANative <- getCPMoA(brayds$ds, brayds$metads, mymodel=braymodel, pertlabels=brayds$metads$Metadata_pert_id, 
                            pclds=lincsMoADS, datalabel="Bray Native MoA, centered full")
  lincs1MoANative <- getCPMoA(lincsds1$ds, lincsds1$metads, mymodel=lincs1model, pertlabels=lincsds1$metads$Metadata_pert_id, 
                              pclds=lincsMoADS, datalabel="LINCS A549 Native MoAs DMSO full model")
  lincs2MoABCNative <- getCPMoA(lincsds2$ds, lincsds2$metads, mymodel=lincs2model, pertlabels=lincsds2$metads$Metadata_pert_id, 
                                pclds=lincsMoADS, datalabel="LINCS Batch2 Native MoAs DMSO full model by cell")
  
  saveRDS(file = file.path(outdir, "../figspaper_res/cpSummaryRes.rds"), 
          object=list(lincsMoADS=lincsMoADS, 
                      brayReps=brayReps,
                      lincs1Reps=lincs1Reps,
                      lincs2Reps=lincs2Reps,
                      lincs2Reps_A549=lincs2Reps_A549,
                      lincs2Reps_MCF7=lincs2Reps_MCF7,
                      lincs2Reps_U2OS=lincs2Reps_U2OS,
                      cpxvalres=cpxvalres, 
                      brayMoA=brayMoA,
                      lincs1MoA=lincs1MoA,
                      lincs2MoA=lincs2MoA,
                      lincs2MoAByCell=lincs2MoAByCell, 
                      lincs2MoAByCellHighDose=lincs2MoAByCellHighDose,
                      brayMoANative=brayMoANative,
                      lincs1MoANative=lincs1MoANative,
                      lincs2MoABCNative=lincs2MoABCNative))
} else {
  cpds <- readRDS(file.path(outdir, "../figspaper_res/cpSummaryRes.rds"))
  attach(cpds)
}



##### Results tables

# xval summary
# use bray5x50
brayAUCML <- sapply(cpxvalres[[2]]$mlRanks, FUN=function(x) 1 - mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(x))))))
brayAUCCos <- sapply(cpxvalres[[2]]$cosRanks, FUN=function(x) 1 - mean(sapply(seq(100), FUN=function(z) mean(unlist(balancedSample(x))))))

brayFDRML <- list(fdr01=sapply(cpxvalres[[2]]$mlRanks, FUN=function(x) mean(sapply(seq(10), FUN=function(z) mean(p.adjust(unlist(balancedSample(x)), "fdr") < 0.01)))),
                  fdr05=sapply(cpxvalres[[2]]$mlRanks, FUN=function(x) mean(sapply(seq(10), FUN=function(z) mean(p.adjust(unlist(balancedSample(x)), "fdr") < 0.05)))),
                  fdr10=sapply(cpxvalres[[2]]$mlRanks, FUN=function(x) mean(sapply(seq(10), FUN=function(z) mean(p.adjust(unlist(balancedSample(x)), "fdr") < 0.10)))))
brayFDRCos <- list(fdr01=sapply(cpxvalres[[2]]$cosRanks, FUN=function(x) mean(sapply(seq(10), FUN=function(z) mean(p.adjust(unlist(balancedSample(x)), "fdr") < 0.01)))),
                   fdr05=sapply(cpxvalres[[2]]$cosRanks, FUN=function(x) mean(sapply(seq(10), FUN=function(z) mean(p.adjust(unlist(balancedSample(x)), "fdr") < 0.05)))),
                   fdr10=sapply(cpxvalres[[2]]$cosRanks, FUN=function(x) mean(sapply(seq(10), FUN=function(z) mean(p.adjust(unlist(balancedSample(x)), "fdr") < 0.10)))))

##### Mechanism of Action Summary #####
moaCPdf <- data.frame()
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(brayMoA, dsname="CDRP MoA", balSampT = 100))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs1MoA, dsname="Lincs1 MoA", balSampT = 100))

moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoAByCell$A549, dsname="Lincs2 MoA A549"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoAByCell$MCF7, dsname="Lincs2 MoA MCF7"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoAByCell$U2OS, dsname="Lincs2 MoA U2OS"))

moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs1MoANative, dsname="Lincs1 Native MoA"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoABCNative$A549, dsname="Lincs2 MoA Native A549"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoABCNative$MCF7, dsname="Lincs2 MoA Native MCF7"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoABCNative$U2OS, dsname="Lincs2 MoA Native U2OS"))

moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoAByCellHighDose$A549, dsname="Lincs2 MoA High Dose A549"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoAByCellHighDose$MCF7, dsname="Lincs2 MoA High Dose MCF7"))
moaCPdf <- rbind(moaCPdf, summarizeCPMoA(lincs2MoAByCellHighDose$U2OS, dsname="Lincs2 MoA High Dose U2OS"))

# fdrdf <- rbind(data.frame(dataset="CDRP Replicate fdr=0.05", method="ML", sp=brayFDRML$fdr05), 
#                data.frame(dataset="CDRP Replicate fdr=0.05", method="Cos", sp=brayFDRCos$fdr05),
#                data.frame(dataset=("CDRP MoA fdr=0.05"), method=moaCPdf$method[1:2], sp=moaCPdf[1:2, 4]), 
#                data.frame(dataset=("CDRP MoA fdr=0.1"), method=moaCPdf$method[1:2], sp=moaCPdf[1:2, 5]), 
#                data.frame(dataset=("CDRP MoA fdr=0.25"), method=moaCPdf$method[1:2], sp=moaCPdf[1:2, 6]))
# fdrdf$dataset <- factor(fdrdf$dataset, levels=c("CDRP Replicate fdr=0.05", "CDRP MoA fdr=0.05", 
#                                                 "CDRP MoA fdr=0.1", "CDRP MoA fdr=0.25"))
# 
# aucdf <- rbind(data.frame(dataset="CDRP Replicate", method="ML", auROC=brayAUCML), 
#                data.frame(dataset="CDRP Replicate", method="Cos", auROC=brayAUCCos),
#                moaCPdf[1:2, c(1,2,3)])



# Make pert figures

aucdf <- rbind(cbind(reshape2::melt(L1KPertBalAUCCos), method="cosine"), 
               cbind(reshape2::melt(L1KPertBalAUCML), method="metric learning"), 
               data.frame(Var1=seq(3), Var2="Cell Painting CDRP", value=brayAUCCos, method="cosine"), 
               data.frame(Var1=seq(3), Var2="Cell Painting CDRP", value=brayAUCML, method="metric learning"))

fdrdf <- rbind(cbind(reshape2::melt(L1KPertBalFDRCos$fdr05), method="cosine"), 
               cbind(reshape2::melt(L1KPertBalFDRML$fdr05), method="metric learning"), 
               data.frame(Var1=seq(3), Var2="Cell Painting CDRP", value=brayFDRCos$fdr05, method="cosine"), 
               data.frame(Var1=seq(3), Var2="Cell Painting CDRP", value=brayFDRML$fdr05, method="metric learning"))

#M is for MoA
Maucdf <- rbind(data.frame(dataset=names(L1KMOABalAUCCos), auROC=as.numeric(L1KMOABalAUCCos), method="cosine"), 
                data.frame(dataset=names(L1KMOABalAUCML), auROC=as.numeric(L1KMOABalAUCML), method="metric learning"), 
                data.frame(dataset="Cell Painting CDRP", auROC=moaCPdf[1:2, 3], method=c("metric learning", "cosine")))

Mfdrdf <- rbind(data.frame(dataset=rownames(L1KMoABalFDRCos), value=as.numeric(L1KMoABalFDRCos$fdr05), method="cosine"), 
                data.frame(dataset=rownames(L1KMoABalFDRML), value=as.numeric(L1KMoABalFDRML$fdr05), method="metric learning"), 
                data.frame(dataset="Cell Painting CDRP", value=moaCPdf[1:2, 4], method=c("metric learning", "cosine")))
Mfdrdf <- rbind(Mfdrdf, 
                data.frame(dataset="CDRP FDR=0.25", value=moaCPdf[1:2, 6], method=c("metric learning", "cosine")))

Maucdf$dataset <- factor(Maucdf$dataset, levels = c(sort(setdiff(Maucdf$dataset, "Cell Painting CDRP")), "Cell Painting CDRP"))
Mfdrdf$dataset <- factor(Mfdrdf$dataset, levels = c(sort(setdiff(Mfdrdf$dataset, c("Cell Painting CDRP", "CDRP FDR=0.25"))), c("Cell Painting CDRP", "CDRP FDR=0.25")))

# Reorder axes, maybe split based on fill so ML on top, cosine on bottom
pdf(file=file.path(outdir, "Nfig2c_ReplicateAUROC.pdf"), width=6, height=5)
ggplot(aucdf, aes(x=value, y=Var2, fill=method, color=method)) + geom_violin() + theme_minimal() + xlim(c(0.5, 1)) + 
  scale_y_discrete(limits=rev) + xlab("AUROC") + ylab("") + ggtitle("Replicate Balanced AUROC")
dev.off()

pdf(file=file.path(outdir, "Nfig2d_ReplicateFDR.pdf"), width=6, height=5)
ggplot(fdrdf, aes(x=value, y=Var2, fill=method, color=method)) + geom_violin() + theme_minimal() + xlim(c(0, 1.1*max(fdrdf$value))) + 
  scale_y_discrete(limits=rev) + xlab("Recall at FDR = 0.05") + ylab("") + ggtitle("Replicate Balanced Recall at FDR = 0.05")
dev.off()

pdf(file=file.path(outdir, "Nfig3a_MoAAUROC.pdf"), width=6, height=5)
ggplot(Maucdf, aes(x=auROC, y=dataset, color=method)) + geom_point(size=2) + theme_minimal() + xlim(c(0.5, 0.8)) + 
  scale_y_discrete(limits=rev) + xlab("AUROC") + ylab("") + ggtitle("MoA Balanced AUROC")
dev.off()

pdf(file=file.path(outdir, "Nfig3b_MoAFDR.pdf"), width=6, height=5)
ggplot(Mfdrdf, aes(x=value, y=dataset, color=method)) + geom_point(size=2) + theme_minimal() +
  scale_y_discrete(limits=rev) + xlab("Recall at FDR = 0.05") + ylab("") + ggtitle("MoA Balanced Recall at FDR = 0.05")
dev.off()


