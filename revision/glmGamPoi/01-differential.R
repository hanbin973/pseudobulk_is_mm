options(warn=-1)
suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(muscat)
    library(nebula)
    library(glmGamPoi)
    library(edgeR)
    library(DESeq2)
    library(Matrix)
})
# get parameters
ns <- as.numeric(snakemake@params[['ns']])
ncps <- as.numeric(snakemake@params[['ncps']])
nc <- ns * ncps
fc <- as.numeric(snakemake@params[['fc']])
in_path <- snakemake@params[['in_path']]

# adjust fc for output
fc_str <- snakemake@params[['fc']]

# read file
sce <- readRDS(in_path)
pb_sce <- glmGamPoi::pseudobulk(
    sce,
    group_by=vars(sample_id, group_id),
    n=n(),
    verbose=FALSE
    )

# === glmGamPoi vs nebula ===

# dataframe for storage
df <- data.frame(matrix(nrow=nrow(sce), ncol=0))
rownames(df) <- rownames(sce)
fit <- glmGamPoi::glm_gp(pb_sce, design = ~1 + group_id)
df[,'glmgampoi'] <- fit$Beta[,2]


# prepare nebula
pred <- model.matrix(
    ~ group_id,
    data=colData(sce)
    )

# run nebula
nebula_result <- nebula::nebula(
    counts(sce),
    colData(sce)$sample_id,
    offset=colMeans2(counts(sce)),
    pred=pred,
    cpc=0,
    mincp=0,
    verbose=FALSE
    )
df[nebula_result$summary$gene,'nebula'] <- nebula_result$summary[,2]

# save result
out_path <- paste0("degs/deg", "_n", ns, "_ncps", ncps, "_fc", fc_str, "_mglmgampoi.csv")
write.csv(df, out_path)
