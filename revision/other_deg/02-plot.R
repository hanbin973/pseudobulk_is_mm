suppressPackageStartupMessages({
    library(Matrix)
    library(tidyverse)
    library(data.table)
    library(matrixStats)
    library(gridExtra)
    library(svglite)
})

# get parameters
ns <- as.numeric(snakemake@params[['ns']])
ncps <- as.numeric(snakemake@params[['ncps']])
nc <- ns * ncps
fc <- snakemake@params[['fc']]
out_path <- snakemake@output[[1]]
in_path <- snakemake@input
print(in_path)

# read file
df1 <- read.csv(in_path[[1]], row.names=1)
df2 <- read.csv(in_path[[2]], row.names=1)

# draw plot
title <- paste0("# subjects=",ns,", mean fold-change=",fc,", # cells per subject=",ncps)
p1 <- ggplot(df1, aes(x=nebula, y=edger), element_text(size=15)) +
    geom_point(alpha=0.5) +
    xlim(c(-3, 3)) +
    ylim(c(-3, 3)) +
    labs(
         x = "logFC - Negative-binomial mixed model",
         y = "logFC - edgeR"
         ) +
    geom_abline(color='red')

p2 <- ggplot(df2, aes(x=nebula, y=deseq2), element_text(size=15)) +
    geom_point(alpha=0.5) +
    xlim(c(-3, 3)) +
    ylim(c(-3, 3)) +
    labs(
         x = "logFC - Negative-binomial mixed model",
         y = "logFC - DESeq2"
         ) +
    geom_abline(color='red')

p3 <- grid.arrange(p1, p2, nrow=1, top=title)

# save image
saveRDS(p3, out_path)
#ggsave(out_path, plot = p3, width=8, height=4)

