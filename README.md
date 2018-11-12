### Introduction
___
The pipeline was used for processing and evaluation of the scCAT-seq datasets. Download our full datasets from here. 

### Dependency
___

* Software (bowtie-v1.2.1.1 , R, bedtools-v2.26.0, python-2.7, samtools-v1.5, SOAP2, perl, HISAT2, bamtools)
* Python packages (MACS2)
* R packages (GenomicAlignments, GenomicFeatures, ggplot2, SummarizedExperiment, reshape2, dplyr, BSgenome, Seurat, ChromVAR,motifmatchr,MotifDb)


### Quick Start
___

```sh
$ git clone --depth=1 xxx.git
$ cd xxx
$ chmod u+x bin/*
$ export PATH=$PATH:./bin/
$ ATACseq_singlecell.pl
	
	perl ATACseq_singlecell.pl -list fq.list -genome hg19.bowtie1_Index.fa -peak peaks.bed -promoter hg19_promoter.bed -blacklist blacklist_regions.bed  -tss hg19_TSS.bed -outdir ./Result

[Options]
  *-list:      list of cells and their fq file paths, format: name <Tab> read1 <Tab> read2, e.g. Cell_1 ATAC_Read_R1.fq.gz ATAC_Read_R2.fq.gz
  *-genome:    genome sequence of species, bowtie index.
   -peak:      reference peaks used for the generation of counting data matrix.
  *-promoter:  promoter regions of genes/transcripts.
  *-blacklist: blacklist regions of genome.
  *-tss:       transcriptional start sites of genes/transcripts.
  *-outdir:    directory for output results. 

 Parameters with * are required, the others are optional.


$ RNAseq_singlecell.pl
	
	perl RNAseq_singlecell.pl -samplelist sample.list -genome human.hisa2.fa -rRNA rRNA.soap2.fa -gff human_Gene_annotation.gff -flatGene refGene.bed -ExonsLength ExonsLength.txt -outdir ./outdir

[options]
	-samplelist:   list of cells and their fq file paths,format: name <Tab> read1 <Tab> read2, e.g. Cell_1 RNA_Read_R1.fq.gz RNA_Read_R2.fq.gz
	-genome:       the HISAT index of genome sequence.
	-rRNA:         the SOAP2 index of rRNA sequence.
	-gff:          GFF file of the genome.
	-bedGene:      bed format of the genes.
	-outdir:       the folder store the results.
	-ExonsLength:  exons length of Genes, format: GeneID <Tab> length.

```

### An example for scCAT-seq data analyses.
___

#### Step1: Convert the scRNAseq fastq files to gene expression data matrix.
The analyses of scRNA-seq part of scCAT-seq data have following steps:

1. Remove the Ribosome reads from fastq via SOAP2 and custom script.
2. Clean fastqs were aligned to human genome via HISAT2.
3. Remove the low mapping quality (MAPQ <30) alignments and unproperly paired reads via samtools.
4. Remove the duplicated fragments via samtools.
5. Count the reads of each gene in each single cell via R.

Here, you can use the perl script to generate gene expresion matrix directly. The sample list contains 3 columns: cell name, the paths of read1 fastq files, and the paths of read2 fastq files. 

```sh
$ RNAseq_singlecell.pl -samplelist sample.list \
                            -genome hg19.hisa2.fa \
                            -rRNA rRNA.soap2.fa \
                            -gff gencode.v19.annotation.gtf \
                            -ExonsLength GeneExonsLength.txt \
                            -bedGene hg19_RefSeq.bed \
                            -outdir ./scRNAseq
```

#### Step2: Convert the scATAC-seq fastq files to deduplecate bam files.
The analyses of scATAC-seq part of scCAT-seq data have following steps:

1. Alignment by bowtie1 
2. Remove the low mapping quality (MAPQ <30) alignments and unproperly paired reads via samtools.
3. Remove duplicate fragments via picard.
4. Remove the fragments which overlap with blacklist regions.
5. Count the usable fragments of each cell within the specific regions (including: promoters and peaks) via bedtools.

Note that if you already have the reference peak file, you can generate the chromatin accessibility count matrix by setting **-peak** parameter in the perl script. The sample name should be the same as that in scRNAseq analysis. 

```sh
$ ATACseq_singlecell.pl -list sample.list \ 
                             -genome hg19.bowtie1.fa \
                             -promoter hg19_promoter.bed \
                             -blacklist hg19_blacklist.bed \
                             -tss hg19_TSS.bed \
                             -outdir ./scATACseq
```


#### Step3: filter out the low quality single cells
Filter the cells according to the threshold below:

* Number of usable fragment of each cell > 10,000
* Ratio of reads in promoter > 15%
* Detected gene number (TMP > 1 ) > 4000 

```R
$ R 

ATAC_Summary <- read.table("./scATACseq/scATACseq_Summary.xls", header=T,
                          col.names=c("Cell", "Raw.Fragments", "Mapped.Fragments", 
                                      "chrM.Fragments", "Usable.Fragments", "Fragments.in.Promoter",
                                       "Perc.Fragments.in.Promoter"))
                          
RNA_Summary <- read.table("./scRNAseq/scRNA_Summary.xls", header=T,
                          col.names=c("Cell", "Raw.Fragments", "Mapped.Ratio", 
                                      "Detected.Gene.Number"))
                          
Summary_df <- merge(ATAC_Summary, RNA_Summary,by="Cell")

col <- rep("gray",length(Summary_df$Cell))

col[which(Summary_df$Usable.Fragments > 1000 & Summary_df$Perc.Fragments.in.Promoter > 15 & Summary_df$Detected.Gene.Number > 4000)] = "black"

usableCell <- Summary_df$Cell[which(col == "black")]

write.table(Summary, file = "./scCAT.data.summary.txt", sep="\t", quote=F, row.names=F)
write.table(as.data.frame(usableCell),file="./scCAT_UsableCell.txt",quote=F)

### you also can filter the matrix of gene expresssion.
Gene_matrix <- read.table("./scRNAseq/scRNAseq_Gene_Count.xls",header=T,row.names=1)
Gene_matrix_usable <- Gene_matrix[,colnames(Gene_matrix) %in% usableCell]
write.table(Gene_matrix_Usable,file="./scRNAseq/scRNA_GeneTPM.UsableCells.xls",quote=F,sep="\t")


### Plot the data quality
plot( log10(Summary_df$Usable.Fragments), Summary_df$Frat.Fragments.in.Promoter, col = col, 
     xlab = "Number of fragments (log10)", ylab = "signal-to-noise ratio", pch = 20)
ablines( v = 4, col = "red", lty = 2)
ablines( h = 15, col = "red", lty = 2)
legend("topright",legend= c("Pass","Falss"), col=c("black","gray"), pch = 20, nrow = 1)

```

#### Step4: If you don't have a reference peak, this step will guide you to construct a reference peak of scATAC-seq and then generate the count matrix.

```sh 
### select the high quality cells
for i in `cat ./scCAT_UsableCell.txt`
  do 
      echo -e "./scATACseq/Result/Alighment/$i/$i.nodup.bam" 
  done
  > ./scATACseq/List/Usable_Cell_bam.list

  
for i in `cat ./scCAT_UsableCell.txt`
  do 
      echo -e "$i\t./scATACseq/Result/Alighment/$i/$i.UsableFragments.bed" 
  done
  > ./scATACseq/List/Usable_Cell_Fragment.list  
  
### merge bam files
bamtools merge -list ./Usable_Cell_bam.list -o ./Usable_Cells.bam

### Peak Calling via MACS2
macs2 callpeak -t ./Usable_Cells.bam \
               -f BAMPE -n all_Cells \ 
               -g hs -B -q 10e-5 \
               --nomodel --nolambda --extsize 250 \
               --outdir ./PeakCalling
               
### count reads in each peak of each single cell
Get_accessibility_Matrix.R  ./PeakCalling/all_Cells_peaks.narrowPeak \
                            ./scATACseq/List/Usable_Cell_Fragnment.list \
                            ./scATACseq/Result/Matrix

```

#### step5: Calculated TF deviation of scATACseq data of Embryo scCAT-seq.
Here, Calculate the TF deviation of scATACseq data of the embryo scCAT-seq via chromVAR

```sh
$ mkdir Embryo_Analysis

$ Rscript TF_deviation.R ./scATACseq/Result/Matrix/UsableCells_PeakMatrix.xls ./Embryo_Analysis/Embryo_TF_Deviation.RDS 

```

#### step6: Normalization of scRNAseq data of embryo scCAT-seq
Here, Normalize the gene expression of scRNAseq data of embryo scCAT-seq via seurat

```sh
$ Rscript Normalize_GeneExp.R ./scRNAseq/scRNA_GeneTPM.UsableCells.xls ./Embryo_Analysis/Embryo_rna_Normalized.RDS

```

#### step7: Plot the accessibility and expression of candidate TFs.

```R
###loading packages
library(Seurat)
library(chromVAR)
library(GenomicRanges)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(motifmatchr)
library(MotifDb)
library(universalmotif)
library(ggplot2)
library(pheatmap)
library(viridis)
library(ggrepel)
library(gridExtra)

### loading data
atac.dev.motifDB <- readRDS("data/scCAT/Embryo_TF_Deviation.RDS")

rownames(atac.dev.motifDB) <- rowData(atac.dev.motifDB)$name
atac.dev.motifDB <- atac.dev.motifDB[!is.na(rownames(atac.dev.motifDB)),]
motif.list <- unique(rownames(atac.dev.motifDB))
motif.list.ix <- sapply(motif.list, function(motif) which(rownames(atac.dev.motifDB) == motif)[1])
atac.dev.motifDB <- atac.dev.motifDB[motif.list.ix,]

TF.variability <- computeVariability(atac.dev.motifDB)

variability.atac.counts.filter <- TF.variability
n = 100

tf.lable <- c("PITX2","HOXD8","MEF2B","TBX2","GRHL1","POU5F1","GATA5","HNF1A","YY1","ONECUT1","KLF1","OTX2","NANOG","GSC","Sox4")

### defined the function.
plotTFvariation <- function(variability.atac.counts.filter,n = 20,tf.lable){
  res_df <- cbind(variability.atac.counts.filter, 
                  rank = rank(-1 * variability.atac.counts.filter$variability, 
                              ties.method = "random"), 
                  annotation = variability.atac.counts.filter$name)
  top_df <- res_df[res_df$name %in% tf.lable, ]
  
  p <- ggplot(res_df, aes_string(x = "rank", 
                      y = "variability", 
                      min = "bootstrap_lower_bound", 
                      max = "bootstrap_upper_bound", 
                      label = "annotation")) + 
          geom_errorbar(colour = "grey") + 
          geom_point(colour = "purple",size = 1) + 
          xlab("Sorted TFs") + ylab("Variability") + 
          theme_bw() + 
          scale_y_continuous(expand= c(0, 0),
                             limits = c(0, max(res_df$bootstrap_upper_bound, na.rm = TRUE) * 1.05)) +
          geom_text_repel(data = top_df, size = 3, col = "Black") # require ggrepel package
  return(p)
}


par(mfrow=c(1,3))

plotTFvariation(TF.variability,50,tf.lable)

# Plot TSNE results
set.seed(1)
tsne_results <- deviationsTsne(atac.dev.motifDB, threshold = 1.2, perplexity = 10, shiny = FALSE)

plotDeviationsTsne(atac.dev.motifDB, tsne_results,sample_column = "stage", shiny = FALSE)

plotDeviationsTsne(atac.dev.motifDB, tsne_results, annotation = "POU5F1", shiny = FALSE)

```

#### step8: Integrative analysis of transcription factor accessibility and RNA expression in the morula and blastocyst

```R 
# loading Normalized RNA data by seurat
rna.seurat <- readRDS("data/scCAT/Embryo_rna_Normalized.RDS") 


colData(atac.dev.motifDB)$group <- "group1"
colData(atac.dev.motifDB)$group[grep("scCAT_Blastocyst",colnames(atac.dev.motifDB))] <- "group2"
colData(atac.dev.motifDB)$group[grep("scCAT_Blastocyst_504|scCAT_Blastocyst_522|scCAT_Blastocyst_539",colnames(atac.dev.motifDB))] <- "group1"


plot_deviation_heatmap <- function(genelist) {
  tf <- genelist
  tf <- tf[tf %in% rownames(atac.dev.motifDB)]
  dat <- assays(atac.dev.motifDB)$deviations[match(tf,rownames(atac.dev.motifDB)),]
  dat <- dat[,colnames(dat)[c(grep("group1",colData(atac.dev.motifDB)$group),
                              grep("group2",colData(atac.dev.motifDB)$group))]]
  dat <- apply(dat, 1, function(tf) {
    tf.box <- boxplot(tf)
    out <- names(tf.box$out)
    tf[out] <- ifelse(tf[out] > 0, tf.box$stats[5,],tf.box$stats[1,])
    return(tf)
  })  
  label = data.frame(group = c(rep("group1",32),rep("group2",40)),
                     row.names = colnames(t(dat)))
  pheatmap(t(dat),
          scale = "row",
          cluster_cols = F,
          cluster_rows = F,
         show_rownames = T, 
         show_colnames = F,
         clustering_method = "ward.D2",
         fontsize_row = 5,
         color = viridis(100),
         border_color = NA,
         annotation_col = label)
}

plot_gene_heatmap <- function(genelist) {
  dat <- rna.seurat@scale.data[,c(grep("group1",colData(atac.dev.motifDB)$group),grep("group2",colData(atac.dev.motifDB)$group))]
  gene_id <- hsapiens.genes$gene_id[match(genelist,hsapiens.genes$gene_name,nomatch = 0)]
  dat <- dat[match(gene_id,rownames(dat),nomatch = 0),]
  rownames(dat) <- hsapiens.genes$gene_name[match(rownames(dat),hsapiens.genes$gene_id,nomatch = 0)]
  dat <- apply(dat, 1, function(tf) {
    tf.box <- boxplot(tf)
    out <- names(tf.box$out)
    tf[out] <- ifelse(tf[out] > 0, tf.box$stats[5,],tf.box$stats[1,])
    return(tf)
  })  
  label = data.frame(group = c(rep("group1",32),rep("group2",40)),
                     row.names = colnames(t(dat)))
  pheatmap(t(dat),
          scale = "row",
          cluster_cols = F,
          cluster_rows = F,
         show_rownames = T, 
         show_colnames = F,
         clustering_method = "ward.D2",
         fontsize_row = 5,
         color = magma(100),
         border_color = NA,
         annotation_col = label)
}

# provide TF list, here is an example list

genelist <- c("NANOG","POU5F1","POU2F1","POU3F1","POU3F2","POU3F3","KLF4","KLF7","KLF14","KLF16","SOX2","SOX5","SOX6","SOX9","MEF2A","MEF2B","MEF2C","MEF2D","OTX1","OTX2","TBX1","TBX3","TBX4","TBX20","HOXD4","HOXD8","HOXD9","HOXD11","GATA2","GATA3","GATA4","GATA6","GRHL1","GRHL2") # An example

plot_deviation_heatmap(genelist)

plot_gene_heatmap(genelist)

```






