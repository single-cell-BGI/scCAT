#!/usr/bin/perl -w 
use strict;

unless(@ARGV==3){
  print "\nperl $0 indir outdir promoter\n\n";
  exit(1);
}

my $indir = shift;
my $outdir =shift;
my $promoter = shift;

my $prefix = (split(/\//,$indir))[-1];
open OUT,">$outdir/$prefix.Summary.xls" || die "$!";
my @samples = `ls $indir/Result/Alignment |grep -Ev "\\.|\\.xls|\\.list"`;
print OUT "Samples\tTotal Fragments\tchrM Fragments\tUsable Fragments\tFragments within promoter\tPercentage of Flagments within promoter (%)\n";
foreach my $sample(@samples){
  chomp($sample);
  my $a = `ls -l $indir/Result/Alignment/$sample/$sample.bowtie.log`;
  my $b = (split(/\s+/,$a))[4];
  next if($b == "0");
  my $all = `grep "reads processed" $indir/Result/Alignment/$sample/$sample.bowtie.log`;
  my $all2 = (split(/\s+/,$all))[3];chomp($all2);
  my $chrM2 = `samtools view -F 4 $indir/Result/Alignment/$sample/$sample.q30.srt.bam | grep "chrM" |wc -l`;
  my $chrM = $chrM2/2;
  my $usable = `wc -l $indir/Result/Alignment/$sample/$sample.UsableFragments.bed | cut -f1 -d " "`;
  my $promoter_reads = `bedtools intersect -a $promoter -b $indir/Result/Alignment/$sample/$sample.UsableFragments.bed -wo |wc -l`;
  if($usable == 0){$usable=1};
  my $promoter_reads_rate = $promoter_reads/$usable*100;
  chomp($chrM);;chomp($usable);chomp($promoter_reads);

  print OUT "$sample\t$all2\t$chrM\t$usable\t$promoter_reads\t$promoter_reads_rate\n";
}close OUT;
