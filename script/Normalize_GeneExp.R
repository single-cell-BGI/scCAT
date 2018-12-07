
args=commandArgs(TRUE)


library(Seurat)
library(SummarizedExperiment)

rna.counts <- read.table(args[1], header = T, stringsAsFactors = F)

rna.seurat <- CreateSeuratObject(raw.data = as.matrix(rna.counts), min.cells = 1, min.genes = 1000, project = "scCAT")
rna.seurat <- NormalizeData(object = rna.embryo, normalization.method = "LogNormalize", scale.factor = 1e4)

rna.seurat <- FindVariableGenes(object = rna.embryo, mean.function = ExpMean, dispersion.function = LogVMR, 
                                x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 1)

rna.seurat <- ScaleData(object = rna.embryo, vars.to.regress = c("nUMI"))

saveRDS(rna.seurat, args[2])