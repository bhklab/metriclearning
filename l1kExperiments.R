


# test metricNTraining
datapath <- "~/Work/bhk/data/l1k/2020/"
cell_id <- "A375"

l1k_meta <- read_l1k_meta(datapath, version=2020)
mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$cell_id == cell_id & l1k_meta$siginfo$pert_type == "trt_cp",]

ds <- parse_gctx(get_level5_ds(datapath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)

resa375 <- metricNTraining(t(ds@mat), mysigs$pert_iname, epochs=2, reps=1)



# Cell Line Specificity

# See how many unique perts each cell line has
# Get those cell lines with at least 1000 unique compounds assayed:
upertcnt <- sapply(names(sort(table(l1kmeta$siginfo$cell_id[l1kmeta$siginfo$pert_type == "trt_cp"]), decreasing=TRUE)), 
                   FUN=function(x) length(unique(siginfo$pert_iname[siginfo$cell_id == x & siginfo$pert_type == "trt_cp"])))

sort(upertcnt[upertcnt >= 1000], decreasing=TRUE)

# Look up the cell lineage of these cell lines:
cellinfo[match(names(upertcnt[upertcnt >= 1000]), cellinfo$cell_id), c(1, 10, 14, 15, 16, 18)]

mycellids <- c("HEPG2", "HCC515", "NPC", "ASC", "HEK293", "YAPC")

grpprts <- unique(siginfo$pert_iname[siginfo$cell_id %in% mycellids & siginfo$pert_type == "trt_cp"])
trainprts <- sample(grpprts, round(length(grpprts)/2))
testprts <- setdiff(grpprts, trainprts)


# Test Cell Line Specificity

clsres_fix <- analyzeL1KCellLineSpecificity(datapath, datapath, cell_ids=mycellids, outpath=".", ncpds=2000, epochs=10, iter=1)
  