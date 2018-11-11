# Usage: Rscript count_expression.r xxx.peak.gff xxx/bamfile/ outfile
# peak.gff was constructed by the reconstructed peak.
# bamfile fold stored the bam file which contain the unique reads.

args= commandArgs(T)

library(GenomicFeatures)
library(GenomicAlignments)

gtf=args[1]
bam=args[2]
geneLength=args[3]
outfile=args[4]

txdb=makeTxDbFromGFF(gtf)
genes=exonsBy(txdb,"gene")
ex=summarizeOverlaps(genes,bam,mode="Union",inter.feature=TRUE,fragment=TRUE,singleEnd=FALSE)
re=as.data.frame(assay(ex))

GeneLength = read.table(geneLength,col.names=c("GeneID","GeneLength"),sep="\t")
GeneLength=GeneLength[match(GeneLength$GeneID,rownames(re)),]
Read_normalize = as.numeric(as.vector(re[,1])) / as.numeric(as.vector(GeneLength[,2]))
re$TPM = Read_normalize / sum(Read_normalize) *1000000

write.table(re,file=outfile,quote=FALSE,sep="\t")
