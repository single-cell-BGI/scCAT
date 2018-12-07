
args = commandArgs(TRUE)


library(chromVAR)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(motifmatchr)
library(MotifDb)
library(universalmotif)


if(args[1] =~ /RDS|rds|Rds/){
  atac.chromVAR = readRDS(args[1])
}else{
  atac.counts <- read.table(args[1], header = T, row.names = 1, stringsAsFactors = F)
  peak.bed <- strsplit(rownames(atac.counts), split = "_") 
  colnames(peak.bed) <- c("chr","start","end")
  peak.gr <- makeGRangesFromDataFrame(peak.bed, keep.extra.columns = T)
  names(peak.gr) <- rownames(atac.counts)
  peak.gr <- peak.gr[match(rownames(atac.counts),names(peak.gr))]
  atac.chromVAR <- SummarizedExperiment(assays = list(counts = as.matrix(atac.counts)),
                                rowRanges = peak.gr)
}

matrices.TF <- query(MotifDb, c('hsapiens', 'mmusculus'))
matrices.TF_v2 <- convert_motifs(matrices.TF, "TFBSTools-PFMatrix")
pfm.list <- do.call(PFMatrixList, matrices.TF_v2)


atac.chromVAR <- addGCBias(atac.chromVAR, genome = BSgenome.Hsapiens.UCSC.hg19)
atac.chromVAR.filter <- filterPeaks(atac.chromVAR, non_overlapping = F)

motif_ix_L <- matchMotifs(pfm.list, 
                          atac.chromVAR.filter, 
                          genome = BSgenome.Hsapiens.UCSC.hg19)

dev_L <- computeDeviations(object = atac.chromVAR.filter, 
                           annotations = motif_ix_L)

saveRDS(dev_L, args[2])