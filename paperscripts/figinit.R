library(gridExtra)
library(ggplot2)
library(cmapR)
library(CMAPToolkit)
#library(dplyr)
library(torch)
library(reshape2)
library(hexbin)
library(RColorBrewer)

### This is a script to generate figures for the manuscript.

topdir <- "~/Work/bhk/analysis/metric_learning/2022_L1K/"
cpdir <- "~/Work/bhk/analysis/metric_learning/2022_cellpaint"
outdir <- file.path(topdir, "figspaper")

datapath <- "~/Work/bhk/data/l1k/2020/"
l1kmeta <- CMAPToolkit::read_l1k_meta(datapath, version=2020)
attach(l1kmeta)

braypath <- "~/Work/bhk/data/cellpainting/bray-2017/bray_2017_combined_profiles.rds"
lincspath <- "~/Work/bhk/data/cellpainting/lincs-cell-painting/lincs-cell-painting/spherized_profiles/profiles"
lincs1 <- file.path(lincspath, "2016_04_01_a549_48hr_batch1_dmso_spherized_profiles_with_input_normalized_by_dmso.rds")
lincs2 <- file.path(lincspath, "2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_dmso.rds")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7") #, "#F0E442"