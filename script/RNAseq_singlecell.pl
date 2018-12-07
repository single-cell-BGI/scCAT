#!/usr/bin/perl 
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);


my($samplelist,$outdir,$genome,$gff,$bedGene,$rRNA,$ExonsLength);

GetOptions(
  "samplelist:s" => \$samplelist,
  "outdir:s" => \$outdir,
  "genome:s" => \$genome,
  "rRNA:s"  => \$rRNA,
  "gff:s" => \$gff,
  "bedGene:s" => \$bedGene,
  "ExonsLength:s" => \$ExonsLength, 
);

if(!defined $samplelist || !defined $outdir || !defined $genome || !defined $gff || !defined $bedGene || !defined $rRNA || !defined $ExonsLength){
   print "\n\tperl $0 -samplelist sample.list -genome human.fa -rRNA rRNA.fa -gff human_Gene_annotation.gff -flatGene refGene.bed -ExonsLength ExonsLength.txt -outdir ./outdir

[options]
\t-list:       list of Cells with fq files,format: Cell_name <Tab> read1.fq.gz <Tab> read2.fq.gz.
\t-genome:     the HISAT index of genome sequence.
\t-rRNA:       the SOAP2 index of rRNA sequence.
\t-gff:        GFF file of the genome.
\t-bedGene:    bed format of the genes.
\t-outdir:     the folder store the results.
\t-ExonsLength:  Exons length of Genes, format: GeneID <Tab> length.\n\n";

  exit(1);
}

mkdir "$outdir" unless(-e "$outdir");
mkdir "$outdir/1_rmrRNA" unless(-e "$outdir/1_rmrRNA");
mkdir "$outdir/2_Alignment" unless(-e "$outdir/2_Alignment");
mkdir "$outdir/3_GeneExp" unless(-e "$outdir/3_GeneExp");
mkdir "$outdir/Shell" unless(-e "$outdir/Shell");

open MAIN,">$outdir/Shell/main.sh" || die "$!";
open PARA,">$outdir/Shell/Part1.sh" || die "$!";
open PARB,">$outdir/Shell/Part2.sh" || die "$!";

print MAIN "sh $outdir/Shell/Part1.sh\n";
print MAIN "sh $outdir/Shell/Part2.sh\n";
close MAIN;

open LIST1,">$outdir/GeneCount.list" || die "$!";


open IN, "$samplelist" || die "$!";
while(<IN>){
   chomp;
   my ($name,$fq1,$fq2)=split(/\s+/,$_);  
   system("mkdir $outdir/1_rmrRNA/$name") unless (-e "$outdir/1_rmrRNA/$name");
   
   open OUT,">$outdir/Shell/$name.sh" || die "$!";
   print PARA "sh $outdir/Shell/$name.sh\n";

   print OUT "soap_mm_gz -a $fq1 -b $fq2 -D $rRNA -m 0 -x 1000 -s 28 -l 31 -v 5 -r 1 -p 3 -o $outdir/1_rmrRNA/$name/$name\_rRNA.PESoap.gz -2 $outdir/1_rmrRNA/$name/$name\_rRNA.SESoap.gz\n";
   print OUT "perl $Bin/rRNAFilter.pl -fq $fq1,$fq2 -soap $outdir/1_rmrRNA/$name/$name\_rRNA.PESoap.gz,$outdir/1_rmrRNA/$name/$name\_rRNA.SESoap.gz -output $outdir/1_rmrRNA/$name/$name\_Clean\n";

   system("mkdir -p $outdir/2_Alignment/$name") unless (-e "$outdir/2_Alignment/$name");

   print OUT "hisat2 -p 8 --phred33 --sensitive --no-discordant --no-mixed -I 1 -X 1000 -x $genome -1 $outdir/1_rmrRNA/$name/$name\_Clean_1.fq.gz -2 $outdir/1_rmrRNA/$name/$name\_Clean_2.fq.gz 2>$outdir/2_Alignment/$name/$name.Hisat2Genome.MapReadsStat.txt | samtools view -b -S - | samtools sort -@ 5 -o $outdir/2_Alignment/$name/$name.genome.bam -\n";
   
   system("mkdir $outdir/3_GeneExp/$name") unless (-e "$outdir/3_GeneExp/$name");
   
   print OUT "samtools view -h -q 30 -f 2 $outdir/2_Alignment/$name/$name.genome.bam |samtools view -bS - | samtools rmdup -S - $outdir/2_Alignment/$name/$name.genome.nodup.bam\n";
   print OUT "samtools index $outdir/2_Alignment/$name/$name.genome.nodup.bam\n";
   print OUT "geneBody_coverage.py -r $bedGene -i $outdir/2_Alignment/$name/$name.genome.nodup.bam -o $outdir/2_Alignment/$name/$name.genome.nodup\n";
   print OUT "Rscript $Bin/count_GeneReadNumber.R $gff $outdir/2_Alignment/$name/$name.genome.nodup.bam $ExonsLength $outdir/3_GeneExp/$name/$name.GeneCounts.txt\n";
   print LIST1 "$name\t$outdir/3_GeneExp/$name/$name.GeneCounts.txt\n";

}close OUT;close LIST1;close PARA;

my ($pre) = `basename $outdir`;
chomp($pre);

print PARB "perl $Bin/RNAseq_Summary.pl -i $outdir -prefix $outdir/$pre.Summary\n";
print PARB "perl $Bin/Get_GeneExpMatrix.pl -SampleList $outdir/GeneCount.list -prefix $outdir/$pre\n";

close PARB;

system("nohup sh $outdir/Shell/main.sh >$outdir/Shell/main.sh.out 2>&1 &");

