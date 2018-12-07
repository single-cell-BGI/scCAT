#!/usr/bin/perl 
use warnings; 
use strict;
use Getopt::Long;

my($SampleList,$prefix);

GetOptions(
   "SampleList:s" => \$SampleList,
   "prefix:s" => \$prefix,
);

if(!defined $SampleList || !defined $prefix ){
   print "\n\tperl $0 -SampleList sample.list -Length ExonsLength.txt -prefix outfile_prefix

[Options]              
     *-SampleList:    list of Gene Expression of Cells, format: Cell_name <Tab> ReadCount_of_Genes.txt.
     *-prefix:        prefix of output, e.g. E:/RNAseq/Test .\n\n";
  exit(1);
}

my(%Count,%TPM);


my (@samples,$name,@genes);

open IN, "$SampleList" || die "$!";
while(<IN>){
   chomp;
   my @a = split(/\s+/,$_);
   if(!-e $a[1]){print "$a[0] is NOT exists.\n";next};

   push @samples,$a[0];
   my $num=0; 
   open SAM, "$a[1]" || die "$!";
   while(<SAM>){
     chomp;
     next if($.==1);
     next if($_=~ /^\s*$/);

     my @tempA=split(/\s+/,$_);

     $Count{$tempA[0]}{$a[0]}=$tempA[1];
     $TPM{$tempA[0]}{$a[0]} = $tempA[2];
     $num+=$tempA[1];    
   }close SAM;
}
close IN;


open OUT ,">$prefix.GeneReadCounts.xls" || die "$!";
open OUT2,">$prefix.GeneTPM.xls" || die "$!";

my $header = join("\t",@samples);
print OUT "GeneID\t$header\n";
print OUT2 "GeneID\t$header\n";

foreach my $gene (sort keys %Count){
     print OUT "$gene";
     print OUT2 "$gene";
     foreach my $sample (@samples){
        if(!defined $Count{$gene}{$sample}){
           print OUT "\t0";
           print OUT2 "\t0";
        }else{
           print OUT "\t$Count{$gene}{$sample}";
           print OUT2 "\t$TPM{$gene}{$sample}";
        }
     }print OUT "\n";print OUT2 "\n";
}close OUT;close OUT2;
