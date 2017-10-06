#!/usr/bin/perl -w

use Getopt::Std;
getopts "v:s:o:m:";


if ((!defined $opt_v)|| (!defined $opt_s) || (!defined $opt_o) ) {
    die "************************************************************************
    Usage: perl $0 -v file.vcf -s sample.included.txt -o out.vcf
      -h : help and usage.
      -v : file.vcf
      -s : included sample name
      -o : output vcf
      -m : missing rate (default 0.4)
************************************************************************\n";
}

my $sample_file = $opt_s;
my $vcf         = $opt_v;
my $outvcf      = $opt_o;
my $missing_rate = (defined $opt_m)?$opt_m:0.4;

my %samdb;
open(INN, $sample_file) or die"";
while(<INN>){
	chomp;
	my @sdb = split(/\s+/,$_);
	map {$samdb{$_}++} @sdb;
	}
close INN;

my %condb;
open(OUT, ">$outvcf") or die"";
open(IN, $vcf) or die"";
while(<IN>){
	chomp;
	if(/##/){
		print OUT "$_\n";
	 }elsif(/#CHROM/){
		my @tmpdb = split(/\s+/,$_);
	  map {$condb{$_}++ if(exists($samdb{$tmpdb[$_]})) } (9..$#tmpdb);
		foreach my $i(0..$#tmpdb){
			print OUT "$tmpdb[$i]	" if($i<9);
		  $condb{$i}++ if(exists($samdb{$tmpdb[$i]}));
		  print OUT "$tmpdb[$i]	" if(exists($samdb{$tmpdb[$i]}));
			}
		print OUT "\n";
	}else{
		my @data = split(/\s+/,$_);
		my $mr   = & missing_rate(@data);
		next if($mr>$missing_rate);
		foreach my $i(0..$#data){
			print OUT "$data[$i]	" if($i<9);
			print OUT "$data[$i]	" if(exists($condb{$i}));
			}
		print OUT "\n";
		}
	}
close IN;


my $cutoff_AD    = 1; ###Used for calculating missing rate
sub missing_rate{
	my @data = @_;
	my $num_of_samples = @data - 8;
	my $num_of_missing_data = 0;
	foreach my $i(9..$#data){
		$num_of_missing_data++ if($data[$i] eq "./.");
#		my @tmpdb = split(/:/,$data[$i]);
#		next if(@tmpdb<2);
#		my $ad = $tmpdb[1];
#		$num_of_missing_data++ if($ad eq ".");
#		next if($ad eq ".");
#		my ($ref_d,$alt_d) = split(/,/,$ad);
#		my $combine_ad     = $ref_d + $alt_d;
#		$num_of_missing_data++ if($combine_ad<$cutoff_AD);
		}
	my $mr_tmp = sprintf("%.2f",$num_of_missing_data/$num_of_samples) if($num_of_samples!=0);
	return $mr_tmp;
	}
