#!/usr/bin/perl -w
###Clean SNPs which are called by GATK and samtools/vcftools together

use Getopt::Std;
getopts "s:g:m:a:b:r:d:";


if ((!defined $opt_s)|| (!defined $opt_g) ) {
    die "************************************************************************
    Usage: perl $0 -s vcf.samtools -g vcf.gatk -r repeatMask.gff3
      -h : help and usage.
      -s : snps/indels vcf file called by samtools/bcftools
      -g : snps/indels vcf file called by GATK
      -r : repeatMasker.gff3
      -m : cutoff missing rate (default 0.4)
      -a : cutoff min depth (default 5)
      -b : cutoff max depth (default 100)
      -d : SNP distance (default 5 bp)
************************************************************************\n";
}else{
  print "************************************************************************\n";
  print "Version 1.1\n";
  print "Copyright to Tanger\n";
  print "RUNNING...\n";
  print "************************************************************************\n";
      }


my $vcf_sam      = $opt_s;
my $vcf_gatk     = $opt_g;
my $repeatLib    = (defined $opt_r)?$opt_r:"NONE";
die "Repeat library should be gff3 format\n" if(!($repeatLib=~/gff/) and $repeatLib ne "NONE");
my $missing_rate = (defined $opt_m)?$opt_m:0.4;
my $depth_min    = (defined $opt_a)?$opt_a:5;
my $depth_max    = (defined $opt_b)?$opt_b:100;
my $SNP_dis      = (defined $opt_d)?$opt_d:5;
my $cutoff_AD    = 1; ###Used for calculating missing rate

my %infordb;
###read GATK result
my %gatkdb;

print "1. reading GATK vcf file ...\n";
open(IN, $vcf_gatk) or die"";
while(<IN>){
	chomp;
	next if(/#/);
	my @data = split(/\s+/,$_);
	my $chrn = $data[0];
	my $posi = $data[1];
	next if((length $data[3])!=1 or (length $data[4])!=1);
	my $key  = $chrn.",".$posi;
	my $ad   = (split/:/,$data[9])[1];
	next if($data[3] =~ /N/);
	next if($data[4] =~ /N/);
  my $GT  = $data[3]."_".$data[4];	
	if(/(FS=\S+);/){
		$strand_bias = $1;
		$strand_bias =~ s/FS=//g;
		$strand_bias =~ s/;.*//g;
		}
	$infordb{$key} += 1;
	next if(!($data[8]=~/AD/));
	my $mr = & missing_rate (@data);	
	my $depth = 0;
	if($data[7]=~/DP=(\d+)/){
		$depth = $1;
	}else{
		$depth = 0;
		}
	$gatkdb{$key}->{'ref'}  = $data[3];
	$gatkdb{$key}->{'alt'}  = $data[4];
	$gatkdb{$key}->{'SB'}   = $strand_bias;
	$gatkdb{$key}->{'MissingRate'} = $mr;
	$gatkdb{$key}->{'depth'} = $depth;
	$gatkdb{$key}->{'GT'}   = $GT;
	}
close IN;

print "2. reading samtools vcf file ...\n";
###read samtools results
my %samdb;
open(IN, $vcf_sam) or die"";
while(<IN>){
	chomp;
	next if(/#/);
	my @data = split(/\s+/,$_);
	$data[4] =~ s/,\<\*\>//g;
	next if((length $data[3])!=1 or (length $data[4])!=1);
	my $chrn = $data[0];
	my $posi = $data[1];
	next if($data[3] =~ /N/);
	next if($data[4] =~ /N/);	
	my $key  = $chrn.",".$posi;
	$samdb{$key}->{'ref'}  = $data[3];
	$samdb{$key}->{'alt'}  = $data[4];
	if(/PV4=(.*)\s+/){
		$strand_bias = (split/,/,$1)[0];
		}
	$samdb{$key}->{'SB'}   = $strand_bias;
  my $GT    = $data[3]."_".$data[4];
  $samdb{$key}->{'GT'}   = $GT;
	$infordb{$key} += 1;
	}
close IN;

my %repdb;
if($repeatLib eq "NONE"){
	$repdb{'none'} = "none";
}else{
  print "3. reading repeat library ...\n";
  open(IN, $repeatLib) or die"";
  while(<IN>){
  	chomp;
  	next if(/#/);
  	my ($ctg,$a,$b) = (split/\s+/,$_)[0,3,4];
  	foreach my $i($a..$b){
  		my $key = $ctg.",".$i;
  		$repdb{$key}++;
  		}
  	}
  close IN;	
	}


print "3. Start filter step ...\n";
print "a. overlap GATK/SAMtools based on SNP location\n";
print "b. remove SNP in repeat regions\n";
print "c. SNP concordance in GATK and SAMtools vcf files\n";
print "d. remove SNPs with strand bias>5 in GATK and samtools\n";
print "e. remove SNPs with depth<=$depth_min and depth>=$depth_max \n";
print "f. remove SNPs with missing rate>$missing_rate\n";
print "g. remove SNPs with <$SNP_dis bp distance\n";


open(BUG, "> debug.txt") or die"";
my $last_posi  = 0;
my $posi       = 0;
open(OUT, "> filterSNP.vcf") or die "";
open(IN, $vcf_gatk) or die"";
while(<IN>){
	chomp;
	if(/#/){
		print OUT "$_\n";
		next;
		}
	$last_posi = $posi;
	my @data  = split(/\s+/,$_);
	my $refN  = $data[3];
	my $altN  = $data[4];	
	next if((length $data[3])!=1 or (length $data[4])!=1);
	my $key   = $data[0].",".$data[1];
	my $lg    = $data[0];
	$posi  = $data[1];
	print BUG "remove $key as SNP distance is smaller than $SNP_dis\n" if((($posi-$last_posi+1)<=$SNP_dis) and ($posi-$last_posi)>0);
	next if((($posi-$last_posi+1)<=$SNP_dis) and ($posi-$last_posi)>0);
	my $gatk_sb = $gatkdb{$key}->{'SB'};
	my $gatk_mr = $gatkdb{$key}->{'MissingRate'};
	my $gatk_dp = $gatkdb{$key}->{'depth'};
	my $gatk_gt = $refN."_".$altN;
	my $sam_sb  = $samdb{$key}->{'SB'};
	my $sam_gt  = $samdb{$key}->{'GT'};
	print BUG "remove $key as it cannot recalled by both gatk and samtools\n" if($infordb{$key}<2);
	next if($infordb{$key}<2);
	print BUG "remove $key as it is in repeat regions\n" if(exists($repdb{$key}));
	next if(exists($repdb{$key}));
	print BUG "remove $key as the SNP in samtools and gatk is not coordinate\n" if($gatk_gt ne $sam_gt);
	next if($gatk_gt ne $sam_gt);
	print BUG "remove $key as strand bias is not defined\n" if(!defined($gatk_sb));
	next if(!defined($gatk_sb));
	next if(!($gatk_sb=~/\d+/) or !($sam_sb=~/\d+/));
  print BUG "remove $key as the strand bias larger than 5 in gatk or samtools vcf\n" if(($gatk_sb>5) or ($sam_sb>5));
	next if(($gatk_sb>5) or ($sam_sb>5));
	print BUG "remove $key as depth is smaller than $depth_min or larger than $depth_max\n" if($gatk_dp<=$depth_min or $gatk_dp>=$depth_max);
	next if($gatk_dp<=$depth_min or $gatk_dp>=$depth_max);
  print BUG "remove $key as missing rate is larger than $missing_rate\n" if($gatk_mr>$missing_rate);
	next if($gatk_mr>$missing_rate);
	print OUT "$_\n";
	}
close IN;
close OUT;


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
