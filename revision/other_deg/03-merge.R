suppressPackageStartupMessages({
    library(Matrix)
    library(tidyverse)
    library(data.table)
    library(matrixStats)
    library(gridExtra)
    library(svglite)
})

# get parameters
fc <- snakemake@params[['fc']]
out_path <- snakemake@output[[1]]
in_path <- snakemake@input

# load data
p_list <- lapply(in_path, readRDS)
p_list <- lapply(p_list, grid.arrange)

# merge plot
p <- grid.arrange(grobs=p_list, ncol=3)
ggsave(out_path, plot=p, width=22, height=22)

