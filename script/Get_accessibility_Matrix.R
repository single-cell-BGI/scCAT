#!/usr/bin/Rscript

args = commandArgs(TRUE)

peak = args[1]
bed = args[2]
out = args[3]


library(SummarizedExperiment)
library(chromVAR)

peaks = read.table(peak,header=F,sep="\t")
peaks = peaks[,1:3]
name_row = paste0(peaks$V1,":",peaks$V2,"_",peaks$V3)
colnames(peaks) = c("chr","start","end")
peaks = makeGRangesFromDataFrame(peaks)

data = read.table(bed,header=F,sep="\t")
reads = as.vector(data$V2)
cellType=as.vector(data$V1)

fragment_counts <- getCounts(reads, peaks,paired =  FALSE,by_rg = FALSE,format = "bed",colData = DataFrame(celltype = c(cellType)))
read_count =as.matrix(assays(fragment_counts)$counts)
rownames(read_count)=name_row
write.table(read_count,file=paste0(out,"/PeakMatrix.xls"),quote=F,sep="\t")
saveRDS(fragment_counts,file=paste0(out,"/PeakMatrix.chromVAR.RDS"))
