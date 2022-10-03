


# test metricNTraining
datapath <- "~/Work/bhk/data/l1k/2020/"
cell_id <- "SKL"

l1k_meta <- read_l1k_meta(datapath, version=2020)
mysigs <- l1k_meta$siginfo[l1k_meta$siginfo$cell_id == "SKL" & l1k_meta$siginfo$pert_type == "trt_cp",]

ds <- parse_gctx(get_level5_ds(datapath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)

res <- metricNTraining(t(ds@mat), mysigs$pert_iname, epochs=10, reps=5)
#  ds <- parse_gctx(fname=get_level5_ds(dspath), rid=l1k_meta$landmarks$pr_gene_id, cid=mysigs$sig_id)

#metricNTraining <- function(mat1, classes, metric="", epochs=10, nvals=c(), validClassLabs=c(), reps=1){
  