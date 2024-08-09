options(warn=-1)
suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(muscat)
})
# read 
# https://stackoverflow.com/questions/68620638/how-to-execute-r-inside-snakemake

# get parameters
ns <- as.numeric(snakemake@params[['ns']])
ncps <- as.numeric(snakemake@params[['ncps']])
nc <- ns * ncps
out_path <- snakemake@params[['out_path']]

# muscat example data
data(example_sce)
ref_sce <- prepSim(example_sce, verbose=FALSE)

# simulate data
sim_sce <- simData(ref_sce,
                   p_dd = c(0,0,1,0,0,0), # only DE genes
                   nk = 1, # number of cell types
                   nc = 2 * nc, # number of cells
                   ns = 2 * ns, # number of subjects
                   force = TRUE
                   )

# mimic case-control
group_A <- paste0('sample', 1:ns, '.A')
group_B <- paste0('sample', (ns+1):(2*ns), '.B')
id_used <- colData(sim_sce)$sample_id %in% c(group_A, group_B)
save_sce <- sim_sce[,id_used]
colData(save_sce)$sample_id <- droplevels(colData(save_sce)$sample_id)

# sort according to sample id - for nebula
id_sort <- order(save_sce$sample_id)
save_sce <- save_sce[,id_sort]

# save data
print(ncol(save_sce))
saveRDS(save_sce, out_path)


