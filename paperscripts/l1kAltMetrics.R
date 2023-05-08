library(fgsea)
library(CMAPToolkit)


generateAltMetrics <- function(datapath, l1kmeta, outpath, pclds, pclcells) {
# Alt metric PCLres
  if (!file.exists(file.path(outpath, "L1K_PCLResAltMetrics.rds"))){
    pclresAlt <- list()
    
    for (ii in seq_along(pclcells)){
      
      acell <- pclcells[ii]
      print(acell)
      
      pclresAlt[[acell]] <- getL1KMoA(datapath, l1kmeta, acell, mymodel=c(), pclds=pclds, altMets=TRUE)
    }
    
    saveRDS(list(pclresAlt=pclresAlt), file=file.path(outpath, "L1K_PCLResAltMetrics.rds"))
    
    return(L1K_PCLResAlt)
  } else {
    L1K_PCLResAlt <- readRDS(file.path(outpath, "L1K_PCLResAltMetrics.rds"))
    return(L1K_PCLResAlt)
  }
}

runGenerateAltMetrics <- function(datapath, metapath, outpath, pcldspath, pclcells=c()){
  
  l1kmeta <- CMAPToolkit::read_l1k_meta(l1kpath, version=2020)

  # Specify which cells to run here
  if (length(pclcells) == 0){
    pclcells <- c("A375", "A549", "HA1E", "MCF10A", "MCF7", "PC3", "VCAP", "ASC", "BT20", "HCC515", "HEK293", "HEPG2", "HUVEC", "JURKAT", "NPC")
  }
  
  pclds <- readRDS(pcldspath)
  
  x <- generateAltMetrics(datapath=datapath, l1kmeta=l1kmeta, outpath=outpath, pclds=pclds, pclcells=pclcells)
  
  print("Executed successfully.")
}