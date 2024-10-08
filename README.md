Repository of "Pseudobulk with proper offsets has the same statistical properties as generalized linear mixed models in single-cell case-control studies" by [Lee and Han (2024)](https://doi.org/10.1093/bioinformatics/btae498).

## A simple example
The following code demonstrates the equivalence between offset pseudobulk and negative-binomial GLMM using the example from nebula repository [link](https://github.com/lhe17/nebula). 
```
# load the libraries
library(SingleCellExperiment)
library(glmGamPoi)
library(nebula)

data(sample_data)

# nebula
df <- model.matrix(~1+cc, data=sample_data$pred)
fit.nebula <- nebula::nebula(
    sample_data$count,
    sample_data$sid,
    pred=df, 
    offset=colMeans(sample_data$count)
)

# glmGamPoi
df <- data.frame(cc = sample_data$pred$cc, sid = sample_data$sid)
sce.obj <- SingleCellExperiment::SingleCellExperiment(
    sample_data$count, 
    colData=df
)
sce.pb <- glmGamPoi::pseudobulk(
    sce.obj,
    group_by=vars(sid, cc),
    verbose=FALSE
)
fit.glmgp <- glmGamPoi::glm_gp(sce.pb, design=~1+cc)

plot(fit.glmgp$Beta[,2],  fit.nebula$summary$logFC_cc)
lines(c(-0.2,0.2),c(-0.2,0.2), col='red')
```

![image](https://github.com/hanbin973/pseudobulk_is_mm/assets/17215340/5e1528a0-0015-478a-be24-a794e0874e89)

