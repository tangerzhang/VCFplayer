#!/usr/bin/perl -w

use Getopt::Std;
getopts "v:s:e:m:o:";


if ((!defined $opt_v)|| (!defined $opt_o) ) {
    die "****************************************************************
    Usage: perl $0 -v snp.vcf -o snp.fasta 
      -h : help and usage.
      -v : snp.vcf, input
      -o : snp.fasta, output
      -s:  sample included
      -e : sample excluded
      -m : missing rate, default 0.4
****************************************************************\n";
}else{
  print "************************************************************************\n";
  print "Version 1.1\n";
  print "Copyright to Tanger\n";
  print "RUNNING...\n";
  print "************************************************************************\n";        
     }

my $vcf    = $opt_v;
my $outfa  = $opt_o;
my $missing_rate = (defined $opt_m)?$opt_m:0.4;
my $sam_in = (defined $opt_s)?$opt_s:"NONE";
my $sam_ex = (defined $opt_e)?$opt_e:"NONE";

my %indb;
my %exdb;
if($sam_ex ne "NONE"){
	open(IN, $sam_ex) or die"";
	while(<IN>){
		chomp;
		my $sample = (split/\s+/,$_)[0];
		$exdb{$sample}++;
		}
	close IN;
	}
	
if($sam_in ne "NONE"){
	open(IN, $sam_in) or die"";
	while(<IN>){
		chomp;
		my $sample = (split/\s+/,$_)[0];
		$indb{$sample}++;
	  }
	close IN;
	}

my %infordb;
my %namedb;
open(IN, "grep -v '##' $vcf |") or die"";
my $head = <IN>;
my @headb = split(/\s+/,$head);
foreach my $i(9..$#headb){
	$namedb{$i} = $headb[$i];
	}
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $chrn = $data[0];
	my $posi = $data[1];
	my $refN = $data[3];
	my $altN = $data[4];
	next if((length $refN)>1 or (length $altN)>1);
	my $mr   = 0;
	my $count = 0;
	foreach my $i(9..$#data){
		my $sample = $namedb{$i};
		next if($sam_ex ne "NONE" and exists($exdb{$sample}));
		next if($sam_in ne "NONE" and !exists($indb{$sample}));
		$count++;
		$mr++ if($data[$i] eq "./.");
		}
	$mr = $mr/$count; $mr = sprintf("%.2f",$mr);
	next if($mr>$missing_rate);
	foreach my $i(9..$#data){
		my $sample = $namedb{$i};
		my $geno   = $data[$i];
		my $baseN;
		if($geno=~/0\/0/){
			$baseN = $refN;
		}elsif($geno=~/0\/1/){
			my $n = int rand 1;
			$baseN = $refN if($n==0);
			$baseN = $altN if($n==1);
		}elsif($geno=~/1\/1/){
			$baseN = $altN;
		}else{
			$baseN = '-';
			}
		next if($sam_ex ne "NONE" and exists($exdb{$sample}));
		next if($sam_in ne "NONE" and !exists($indb{$sample}));
		$infordb{$i}->{'sample'} = $sample;
	  $infordb{$i}->{'seq'}   .= $baseN;
		}
	}
close IN;

open(OUT, "> $outfa") or die"";
foreach my $i (sort {$a<=>$b} keys %infordb){
	my $sample = $infordb{$i}->{'sample'};
	my $seq    = $infordb{$i}->{'seq'};
	print OUT ">$sample\n$seq\n";
	}
close OUT;
