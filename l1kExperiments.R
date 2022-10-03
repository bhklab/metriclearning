


# test metricNTraining
datapath <- "~/Work/bhk/data/l1k/2020/"
cell_id <- "SKL"

l1k_meta <- read_l1k_meta(datapath, version=2020)
mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$cell_id == "SKL" & l1k_meta$siginfo$pert_type == "trt_cp",]

ds <- parse_gctx(get_level5_ds(datapath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)

res0 <- metricNTraining(t(ds@mat), mysigs$pert_iname, epochs=5, reps=1)
