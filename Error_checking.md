# Some comparisons between the RADseq and HiSeq individuals indentified several concerns.

First, there appears to be some heterozygous chrX sites in males that snuck through the RADseq pipeline for some reason.  I wrote this script to identify how many there are in each male and also to quantify how many there are at eash position (Counts_chrX_het_sites_in_males.pl):

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;

# this script counts male chrX het sites in a tab delimited file.
# TO run type this
# Counts_chrX_het_sites_in_males.pl inputfile.tab 1111100110000111100011100110010100000000 3_6 outputfile.tab 
# where 1111100110000111100011100110010100000000 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) or is not (0) female ,and 3_6 refers to (i) the column that contains the 
# outgroup nucleotide (3 in this case for rhesus) and (ii) the column number of the first individual in the ingroup 
# (6 in this case which would be a file with human and baboon outgroup seqs already in it)

# Counts_chrX_het_sites_in_males.pl chrX_HiSeq_RADseq_combined.tab 1111100110000111100011100110010100000000 3_6 positions_of_male_chrX_hets.tab 


my $inputfile = $ARGV[0];
my @rhesus_and_begin = split("_",$ARGV[2]);
my $outputfile = $ARGV[3];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

my @temp;

# first define the males (0s)
my @sexes=split('',$ARGV[1]);
my $counter;
my @males;
my @bases;
my %hetsites;

foreach(@sexes){
	if($_ == 0){
		push(@males,$counter);
	}
	$counter+=1;
}

# so @males is an non-continuous array to indicate which column we have male genotypes
# we need to add $rhesus_and_begin[1]-1 to this column to get the actual column in the file

my @namez;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if(($temp[0] ne '#CHROM')&&($temp[0] eq 'chrX')){
		$counter=0;
		foreach(@males){
			@bases=split('/',$temp[$_+$rhesus_and_begin[1]-1]);
			if($bases[0] ne $bases[1]){
				$hetsites{$temp[0].'_'.$temp[1]}[$counter]=1;
			}
			$counter+=1;
		}
		# so $hetsites{chr_pos}[] refers to each male sequentially
	}
	elsif($temp[0] eq '#CHROM'){
		$counter=0;
		foreach(@males){
			push(@namez,$temp[$_+$rhesus_and_begin[1]-1]);
		}
		print "@namez\n";
		# so names is an array of length $#males		
	}
}
close DATAINPUT;

my @number_of_hets;

# initialize the array
$counter=0;
foreach(@namez) {
	$number_of_hets[$counter]=0;
	$counter+=1;
}

# so the @number_of_hets should have each male sequentially

# now cycle through each key to add up the het sites by individual
foreach my $chr_and_position (keys %hetsites ) {
	foreach my $individual (@ {$hetsites{$chr_and_position}} ) {
		if($hetsites{$chr_and_position}[$individual] == 1){
			$number_of_hets[$individual]+=1;
		}
	}
}

my %number_of_heterozygous_males_per_chrX_site;

# now cycle through each key of the array to add up the het sites by site
foreach my $chr_and_position (keys %hetsites ) {
	foreach my $individual (@ {$hetsites{$chr_and_position}} ) {
		if($hetsites{$chr_and_position}[$individual] == 1){
			$number_of_heterozygous_males_per_chrX_site{$chr_and_position}+=1;
		}
	}
}

$counter=0;
print "Here is a summary of the number of heterozygous sites in males\n";
foreach(@number_of_hets){
	print $namez[$counter],"\t",$number_of_hets[$_],"\n";
	$counter+=1;
}

my @temp2;

print "Here is a summary of the sites with one or more heterozygous male\n";
foreach my $chr_and_position (keys %number_of_heterozygous_males_per_chrX_site ) {
	if(defined($number_of_heterozygous_males_per_chrX_site{$chr_and_position})){
		@temp2=split('_',$chr_and_position);
		print OUTFILE $temp2[0],"\t",$temp2[1],"\t",$number_of_heterozygous_males_per_chrX_site{$chr_and_position},"\n";
	}
}

close OUTFILE;

```


