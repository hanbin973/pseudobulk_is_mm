options(warn=-1)
suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(nebula)
    library(glmGamPoi)
    library(Matrix)
})

# get parameters
ncps <- snakemake@params[['ncps']]
in_path <- snakemake@input
out_path <- snakemake@output[[1]]

# store
df <- data.frame(
                 matrix(nrow=0, ncol=3)
                 )
colnames(df) <- c('nc', 'method', 'time')

# benchmark
nit <- 10
it <- 1 
for (path in in_path){
    sce <- readRDS(path)
    pb_sce <- glmGamPoi::pseudobulk(
                                    sce,
                                    group_by=vars(sample_id, group_id),
                                    n=n(),
                                    verbose=FALSE
                                    )
    pred <- model.matrix(
                         ~group_id,
                         data=colData(sce)
                         )
    for (i in 1:nit){
        # glmGamPoi
        t1 <- Sys.time()
        fit <- glmGamPoi::glm_gp(pb_sce, design=~1+group_id)
        test <- glmGamPoi::test_de(fit, reduced_design=~1)
        t2 <- Sys.time()
        df[nrow(df) + 1,] <- c(ncps[it], 'glmgampoi', t2-t1)

        # nebula
        t1 <- Sys.time()
        fit <- nebula::nebula(
                              counts(sce),
                              colData(sce)$sample_id,
                              offset=colMeans2(counts(sce)),
                              pred=pred,
                              cpc=0,
                              mincp=0,
                              verbose=FALSE
                              )
        t2 <- Sys.time()
        df[nrow(df) + 1,] <- c(ncps[it], 'nebula', t2-t1)
    }
    it <- it+1
}

# save file
write.csv(df, out_path)

