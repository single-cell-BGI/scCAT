#!/usr/bin/perl 
use warnings;
use strict;
use Getopt::Long;

my ($in_dir,$out_prefix);

GetOptions(
  "i:s" => \$in_dir,
  "prefix:s" => \$out_prefix,
);

if(!defined $in_dir || !defined $out_prefix){
   print "\n\tperl $0 -i in_dir -prefix file_prefix\n\n";
   exit(1);
}

my %summery;

my @samples = `ls $in_dir/1_rmrRNA `;

open OUT,">$out_prefix.xls" || die "$!";

###0_remove rRNA result

print OUT "\tRaw reads\tGenome mapped rate(%)\tGene number(TPM >=1)\n";

foreach my $sample (@samples){
   chomp($sample);
   if(!-e "$in_dir/1_rmrRNA/$sample/$sample\_Clean.rrna.stat"){
   $summery{$sample} .= "0";
   }else{
      my $a = `cat $in_dir/1_rmrRNA/$sample/$sample\_Clean.rrna.stat |grep total`;
      my $raw_reads = `cat $in_dir/1_rmrRNA/$sample/$sample\_Clean.rrna.stat |grep total |cut -f2`;
      chomp($raw_reads);
   $summery{$sample} .= "$raw_reads";
   }
 
  my ($rate_genome);
  if(!-e "$in_dir/2_Alignment/$sample/$sample.Hisat2Genome.MapReadsStat.txt"){
     print "$sample not complete at genome alignment\n";
     $rate_genome=0;
     $summery{$sample} .= "\t$rate_genome";
  }else{
     open IN,"<$in_dir/2_Alignment/$sample/$sample.Hisat2Genome.MapReadsStat.txt" || die "$!";
      while(<IN>){
        chomp;
        if($_ =~ /overall alignment rate/){
           my @a = split(/\s+/,$_);
           $a[0]=~ s/%//;
           $rate_genome = $a[0];
        }else{$rate_genome=0;}
      }close IN;

    $summery{$sample} .= "\t$rate_genome";
 
  }

  my ($gene_number);  
  if(!-e "$in_dir/3_GeneExp/$sample/$sample.GeneCounts.txt"){
     print "$sample not complete at Count Gene ReadNumber\n";
     $gene_number=0;
     $summery{$sample} .= "\t$gene_number";
  }else{
     $gene_number = `awk '\$3>=1 && NR>1 {print }' $in_dir/3_GeneExp/$sample/$sample.GeneCounts.txt |wc -l`;
     chomp($gene_number);
     $summery{$sample} .= "\t$gene_number";
  }
  print OUT "$sample\t$summery{$sample}\n";
}close OUT;
