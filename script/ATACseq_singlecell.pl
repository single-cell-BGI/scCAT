#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);

my ($list,$outdir,$promoter,$genome,$blacklist,$peak,$tss);

GetOptions(
  "list:s" => \$list,
  "outdir:s" => \$outdir,
  "peak:s" => \$peak,
  "genome:s" => \$genome,
  "promoter:s" => \$promoter,
  "blacklist:s" => \$blacklist,
  "tss:s" => \$tss,
);

 
if(!defined $list || !defined $outdir || !defined $promoter || !defined $genome || !defined $blacklist || !defined $tss){
   print "\n\tperl $0 -list sample.list -genome hg19.bowtie1.fa -promoter -outdir out\n\n";
   print "[Options]
\t*-list:      a list of sample, format: name <Tab> read1.fq.gz <Tab> read2.fq.gz.
\t*-genome:    genome sequence of species, bowtie1 index.
\t -peak:      standard peak for count coverage in sinle cells, bed formation. 
\t*-promoter:  promoter region of genes/transcripts.
\t*-blacklist: blacklist regions of genome.
\t*-tss:       transcript start site of genes/transcripts.
\t*-outdir:    output fold.

The parameter with * is required\n\n";
   exit(1);
}

mkdir "$outdir" unless(-e "$outdir");
system("mkdir -p $outdir/Result/Alignment") unless( -e "$outdir/Result/Alignment");
system("mkdir -p $outdir/Result/TSS_Enrichment")unless (-e "$outdir/Result/TSS_Enrichment");
system("mkdir -p $outdir/Result/Matrix") unless(-e "$outdir/Result/Matrix");
mkdir "$outdir/Shell";
mkdir "$outdir/List";
system("mkdir -p $outdir/Result/Alignment/") unless (-e "$outdir/Result/Alignment/");
system("mkdir -p $outdir/Result/TSS_Enrichment/") unless (-e "$outdir/Result/TSS_Enrichment/");

open SH,">$outdir/main.sh" || die "$!";
open SHA,">$outdir/Shell/FirstPart.sh" || die "$!";
open SHB,">$outdir/Shell/SecondPart.sh" || die "$!";

open LIST,">$outdir/List/PeakCount.list" || die "$!";
open LIST2,">$outdir/List/UsableFragnment.list" || die "$!";

print SH "sh $outdir/Shell/FirstPart.sh\nsh $outdir/Shell/SecondPart.sh\n";

open IN,"$list" || die "$!";
while(<IN>){
     chomp;
     my ($sample,$fq1,$fq2) = split(/\s+/,$_); 
     mkdir "$outdir/Result/Alignment/$sample" unless(-e "$outdir/Result/Alignment/$sample");
     mkdir "$outdir/Result/TSS_Enrichment/$sample" unless(-e "$outdir/Result/TSS_Enrichment/$sample"); 
	  
     open OUT ,">$outdir/Shell/$sample\_ATACline_run.sh" || die "$!"; 
     print SHA "sh $outdir/Shell/$sample\_ATACline_run.sh\n";
     print OUT "gunzip -dc $fq1 > $outdir/Result/Alignment/$sample/$sample\_R1.fq\n";
     print OUT "gunzip -dc $fq2 > $outdir/Result/Alignment/$sample/$sample\_R2.fq\n";
     print OUT "bowtie -X 2000 -m 1 -l 25 -p 15 --chunkmbs 256 $genome -1 $outdir/Result/Alignment/$sample/$sample\_R1.fq -2 $outdir/Result/Alignment/$sample/$sample\_R2.fq  2>$outdir/Result/Alignment/$sample/$sample.bowtie.log -S | samtools view -bS -q 30 -f 2 - | samtools sort -@ 5 -o $outdir/Result/Alignment/$sample/$sample.q30.srt.bam -\n";
     print OUT "rm $outdir/Result/Alignment/$sample/$sample\_R1.fq\n";
     print OUT "rm $outdir/Result/Alignment/$sample/$sample\_R2.fq\n";
     print OUT "java -Xmx2g -jar $Bin/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true INPUT=$outdir/Result/Alignment/$sample/$sample.q30.srt.bam OUTPUT=$outdir/Result/Alignment/$sample/$sample.nodup.bam M=$outdir/Result/Alignment/$sample/$sample.sort.rmdupMetrics.xls\n";
     print OUT "samtools index $outdir/Result/Alignment/$sample/$sample.nodup.bam\n";
     print OUT "bamCoverage -b $outdir/Result/Alignment/$sample/$sample.nodup.bam -o $outdir/Result/Alignment/$sample/$sample.nodup.bw\n";
     print OUT "computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 -R $tss -S $outdir/Result/Alignment/$sample/$sample.nodup.bw --skipZeros -o $outdir/Result/TSS_Enrichment/$sample/$sample.nodup.TSS_Matrix.gz\n";
     print OUT "plotProfile -m $outdir/Result/TSS_Enrichment/$sample/$sample.nodup.TSS_Matrix.gz -out $outdir/Result/TSS_Enrichment/$sample/$sample.nodup.TSSenrichment.pdf\n";
     print OUT "samtools sort -@ 5 -n $outdir/Result/Alignment/$sample/$sample.nodup.bam -o $outdir/Result/Alignment/$sample/$sample.nodup.name.bam\n";
     print OUT "bedtools bamtobed -i $outdir/Result/Alignment/$sample/$sample.nodup.name.bam -bedpe | grep -v \"chrM\" | awk '{print \$1\"\\t\"\$2\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10}' >$outdir/Result/Alignment/$sample/$sample.nodup.Fragments.bed\n";
     print OUT "rm $outdir/Result/Alignment/$sample/$sample.nodup.name.bam\n";
     print OUT "bedtools intersect -a $outdir/Result/Alignment/$sample/$sample.nodup.Fragments.bed -b $blacklist -v  | sort -k 1,1 -k 2,2n > $outdir/Result/Alignment/$sample/$sample.UsableFragments.bed\n";
     print LIST2 "$sample\t$outdir/Result/Alignment/$sample/$sample.UsableFragments.bed\n";
     close OUT;
}close IN;close LIST;close LIST2;
    
    my $prefix = `basename $outdir`;
    print SHB "perl $Bin/ATACseq_Summary.pl $outdir $outdir $promoter\n";
    if(defined $peak){
       print SHB "Rscript ReadCount.R $peak $outdir/List/UsableFragnment.list $outdir/Result/Matrix\n";
    }
    close SHA;close SHB;


print SH "sh $outdir/Shell/FirstPart.sh\nsh $outdir/Shell/SecondPart.sh\n";
close SH;

system("nohup sh $outdir/main.sh >$outdir/Shell/main.sh.out 2>&1 &")
