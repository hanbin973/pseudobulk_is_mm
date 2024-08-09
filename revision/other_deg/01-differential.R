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

# === edgeR vs nebula ===

# dataframe for storage
df <- data.frame(matrix(nrow=nrow(sce), ncol=0))
rownames(df) <- rownames(sce)

# edgeR
model_matrix <- model.matrix(~ group_id, colData(pb_sce))
edgeR_data <- edgeR::DGEList(counts(pb_sce))
edgeR_data <- edgeR::estimateDisp(edgeR_data, model_matrix)
edgeR_data <- edgeR::calcNormFactors(edgeR_data)
edgeR_fit <- edgeR::glmFit(edgeR_data, design=model_matrix)
df[,'edger'] <- edgeR_fit$coefficients[,2]

# prepare nebula
pred <- model.matrix(
    ~ group_id,
    data=colData(sce)
    )

# deseq-subject offset to cell-offset
sp <- sparse.model.matrix(
    ~ 0 + sample_id,
    data=colData(sce)
    )
offset_subject <- as.vector(exp(edgeR_fit$offset)) / pb_sce$n
offset_cell <- (sp %*% offset_subject)[,1]

# run nebula
nebula_result <- nebula::nebula(
    counts(sce),
    colData(sce)$sample_id,
    offset=offset_cell,
    pred=pred,
    cpc=0,
    mincp=0,
    verbose=FALSE
    )
df[nebula_result$summary$gene,'nebula'] <- nebula_result$summary[,2]

# save result
out_path <- paste0("degs/deg", "_n", ns, "_ncps", ncps, "_fc", fc_str, "_medger.csv")
write.csv(df, out_path)




# === DESeq2 vs nebula ===

# dataframe for storage
df <- data.frame(matrix(nrow=nrow(sce), ncol=0))
rownames(df) <- rownames(sce)

# DESeq2
dds <- DESeq2::DESeqDataSetFromMatrix(
    counts(pb_sce),
    colData = colData(pb_sce),
    design = ~ group_id
    )
dds <- DESeq2::estimateSizeFactors(dds)
dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)
dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6)
# change log base from 2 to e
df[,'deseq2'] <- rowData(dds)$group_id_B_vs_A * log(2)

# prepare nebula
pred <- model.matrix(
    ~ group_id,
    data=colData(sce)
    )

# deseq-subject offset to cell-offset
sp <- sparse.model.matrix(
    ~ 0 + sample_id,
    data=colData(sce)
    )
offset_subject <- sizeFactors(dds) / pb_sce$n
offset_cell <- (sp %*% offset_subject)[,1]


# run nebula
nebula_result <- nebula::nebula(
    counts(sce),
    colData(sce)$sample_id,
    offset=offset_cell,
    pred=pred,
    cpc=0,
    mincp=0,
    verbose=FALSE
    )
df[nebula_result$summary$gene,'nebula'] <- nebula_result$summary[,2]

# save result
out_path <- paste0("degs/deg", "_n", ns, "_ncps", ncps, "_fc", fc_str, "_mdeseq2.csv")
write.csv(df, out_path)
