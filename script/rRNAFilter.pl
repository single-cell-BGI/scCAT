#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename qw(dirname);
use Cwd 'abs_path';
use PerlIO::gzip;
use File::Path;
use FindBin;
use lib("$FindBin::RealBin/PerlLib");

my ($fq, $soap, $output);

GetOptions(
   "fq:s" => \$fq, 
   "soap:s" => \$soap, 
   "output:s" => \$output, 
);

if (!defined $fq || !defined $soap || !defined $output) {
  print"
description: filter rRNA
usage: perl $0 [options]
options:
	*-fq:        fq files, separated by comma \",\"
	*-soap:      soap result files, separated by comma \",\"
	*-output:    prefix of output file.
e.g.
	perl $0 -fq a_1.fq,a_2.fq -soap soap.pe,soap.se -output new.a
	perl $0 -fq a.fq -soap soap.pe -output new.a\n";
	exit (1);
}

my @fqs = split /,+/, $fq;
my @soaps = split /,+/, $soap;

my $exit = 0;
for ((@fqs, @soaps)) {
	if (!-f $_) {
		print STDERR "file $_ not exists\n";
		$exit = 1;
	}
}

exit 1 if ($exit == 1);

my $outdir = dirname($output);
mkdir $outdir if (!-d $outdir);

my %rRNAs = ();
&getRRNA(\@soaps, \%rRNAs);
&filterReads(\@fqs, \%rRNAs, $output);
&showLog("done");

exit 0;

sub showLog {
	my ($info) = @_;
	my @times = localtime; # sec, min, hour, day, month, year
	print STDERR sprintf("[%d-%02d-%02d %02d:%02d:%02d] %s\n", $times[5] + 1900, $times[4] + 1, $times[3], $times[2], $times[1], $times[0], $info);
}

sub getRRNA {
	my ($soaps, $rRNAs) = @_;
	for my $soap (@{$soaps}) {
		&showLog("read soap file: $soap");
		#open SOAP, "< $soap" || die $!;
		my $fileformat=`file $soap `;
		if ($fileformat =~ /gzip compressed data/) {
			open SOAP, " gzip -dc $soap |" or die "can't open the soap file  $soap";
		}elsif ($fileformat =~ /bzip2 compressed data/){
			open SOAP, " bzip2 -cd $soap |" or die "can't open the soap file  $soap";
		}else{
			open SOAP, $soap or die "can't open the soap file  $soap";
		}
		while (<SOAP>) {
			next if (/^\s*$/);
			chomp;
			my @tabs = split /\t/, $_;
			$tabs[0] =~ s/\/.*//;
			$rRNAs->{$tabs[0]} = 1;
		}
		close SOAP;
	}
}

sub filterReads {
	my ($fqs, $rRNAs, $output) = @_;
	my ($i, $total, $delete);
	for my $fq (@{$fqs}) {
		&showLog("process reads file: $fq");
		$total = 0;
		$delete = 0;
		my $suffix = ".fq.gz";
		$i++;
		if (@{$fqs} > 0) {
			$suffix = "_$i.fq.gz";
		}
		my $fileformat=`file $fq `;
		if ($fileformat =~ /gzip compressed data/) {
                        open IN, " gzip -dc $fq |" or die "can't open the soap file  $fq";
                }elsif ($fileformat =~ /bzip2 compressed data/){
                        open IN, " bzip2 -cd $fq |" or die "can't open the soap file  $fq";
                }elsif($fq =~ /.gz$/){ # for link file...
			open IN, " gzip -dc $fq |" or die "can't open the soap file  $fq";
		}else{
                        open IN, $fq or die "can't open the soap file  $fq";
                }
		open OUT, ">:gzip", "$output$suffix" || die $!;
        while (<IN>) {
			$total++;
			my $name = $_;
			my $seq = <IN>;
			<IN>;
			my $quality = <IN>;
			chomp(my $temp = $name);
			if($temp =~ m/\//){ # hiseq2000->old read id
				$temp =~ s/^\@(.*)\/.*/$1/; 
			}elsif($temp =~ m/\s+(.+)/){ # hiseq4000->new read id
				$temp =~ s/^\@(.*)\s+.*/$1/; 
			}else{
				$temp =~ s/^\@(.*)$/$1/; # CG platform read id
			}
			if (exists $rRNAs->{$temp}) {
				$delete++;
			} else {
				print OUT "$name$seq+\n$quality";
			}
		}
		close OUT;
		close IN;
	}
	&showLog("output stat file: $output.rrna.stat");
	my $percent = sprintf("%.2f%%", $delete / $total * 100);
	open STAT, "> $output.rrna.stat" || die $!;
	print STAT << "EOF";
total\t$total
delete\t$delete\t$percent
EOF
	close STAT;
}
