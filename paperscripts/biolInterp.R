library(fgsea)
library(CMAPToolkit)
library(cmapR)



biolInterp <- function(gopath, rotmat, landmarks, nPCs=300, minGoSize=30){
  
  geneRes <- list()
  
  goPws <- parse_gmt(file.path(gopath, "c5.go.v2023.1.Hs.entrez.gmt"))
  
  goLMs <- sapply(goPws, FUN=function(x) length(intersect(x$entry, landmarks$pr_gene_id)))
  
  goPwsLMs <- lapply(goPws, FUN=function(x) intersect(x$entry, landmarks$pr_gene_id))
  
  
  for (pcix in seq(min(nPCs, 978))){
    print(pcix)
    fgseaRes <- fgsea(goPwsLMs[goLMs >= minGoSize], rotmat[,pcix])
    geneRes[colnames(rotmat)[pcix]] <- list(fgseaRes[, 1:(dim(fgseaRes)[2]-1)])
  }
  
  return(geneRes)
}


runBioInterp <- function(outpath, eigenDataPath, l1kmetapath, gopath, nPCs=300){
  
  l1kmeta <- CMAPToolkit::read_l1k_meta(l1kmetapath, version=2020)
  attach(l1kmeta)
  
  eigendata <- readRDS(eigenDataPath)
  
  bioRes <- list()
  
  for (mycell in names(eigendata)){
    print(mycell)
    bioRes[[mycell]] <- biolInterp(gopath=gopath, rotmat=eigendata[[mycell]]$pcBase$rotation, landmarks=landmarks, nPCs=nPCs)
  }
  
  saveRDS(object=bioRes, file=file.path(outpath, "bioInterpRes.rds"))
}
