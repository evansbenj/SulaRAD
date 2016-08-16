# Some comparisons between the RADseq and HiSeq individuals indentified several concerns.

First, there appears to be some heterozygous chrX sites in males that snuck through the RADseq pipeline for some reason.  I wrote this script to identify how many there are in each male and also to quantify how many there are at eash position (Counts_chrX_het_sites_in_males.pl):

```perl
#!/usr/bin/env perl
use strict;
use warnings;

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

foreach(@sexes){
	if($_ == 0){
		push(@males,$counter);
	}
	$counter+=1;
}

# so @males is an non-continuous array to indicate which column we have male genotypes
# we need to add $rhesus_and_begin[1]-1 to this column to get the actual column in the file

my @namez;
my @number_of_hets;
my %number_of_heterozygous_males_per_chrX_site;
my @number_of_het_sites_per_male;
my $number_of_sites_with_at_least_one_het_male=0;

# initialize @number_of_het_sites_per_male
$counter=0;
foreach(@males){
	$number_of_het_sites_per_male[$counter]=0;
	$counter+=1;
}	

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if(($temp[0] ne '#CHROM')&&($temp[0] eq 'chrX')){
		$counter=0;
		foreach(@males){
			@bases=split('/',$temp[$_+$rhesus_and_begin[1]-1]);
			if($bases[0] ne $bases[1]){
				$number_of_heterozygous_males_per_chrX_site{$temp[0].'_'.$temp[1]}+=1;
				$number_of_het_sites_per_male[$counter]+=1;
			}
			$counter+=1;
		}
		if(defined($number_of_heterozygous_males_per_chrX_site{$temp[0].'_'.$temp[1]})){
			$number_of_sites_with_at_least_one_het_male+=1;
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

print "The number_of_sites_with_at_least_one_het_male is ",$number_of_sites_with_at_least_one_het_male,"\n";

foreach my $chr_and_position (keys %number_of_heterozygous_males_per_chrX_site ) {
	print $chr_and_position,"\t",$number_of_heterozygous_males_per_chrX_site{$chr_and_position},"\n";
}	

$counter=0;
foreach (@number_of_het_sites_per_male) {
	if(defined($number_of_het_sites_per_male[$counter])){
		print $namez[$counter],"\t",$number_of_het_sites_per_male[$counter],"\n";
	}
	else{
		print $namez[$counter],"\t0\n";
	}
	$counter+=1;	
}


```

Unfortunately there are alot (441 to be exact) and multiple in every of the 22 male individuals (range 28-119).


Another concern is that, in the HiSeq data there are lots of heterozygous sites on the chrX in males.  I wrote this script to identify them and to compare the RADseq nad HiSeq genotypes (Compares_RADseq_and_Hiseq.pl)"

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;

# this script compares RAdseq and HiSeq genotypes in a file.
# TO run type this
# Compares_RADseq_and_Hiseq.pl chrX_HiSeq_RADseq_combined.tab

my $inputfile = $ARGV[0];
my $outputfile3 = "nem_diffs.txt";
my $outputfile4 = "nigra_diffs.txt";
my $outputfile5 = "tonk_diffs.txt";
my $outputfile6 = "X_hiseq_het_sites_in_PM664.txt";
my $outputfile7 = "X_hiseq_het_sites_in_PF660.txt";
my $outputfile8 = "X_hiseq_het_sites_in_PM592.txt";

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}


unless (open(OUTFILE3, ">$outputfile3"))  {
	print "I can\'t write to $outputfile3\n";
	exit;
}
print "Creating output file: $outputfile3\n";

unless (open(OUTFILE4, ">$outputfile4"))  {
	print "I can\'t write to $outputfile4\n";
	exit;
}
print "Creating output file: $outputfile4\n";

unless (open(OUTFILE5, ">$outputfile5"))  {
	print "I can\'t write to $outputfile5\n";
	exit;
}
print "Creating output file: $outputfile5\n";

unless (open(OUTFILE6, ">$outputfile6"))  {
	print "I can\'t write to $outputfile6\n";
	exit;
}
print "Creating output file: $outputfile6\n";

unless (open(OUTFILE7, ">$outputfile7"))  {
	print "I can\'t write to $outputfile7\n";
	exit;
}
print "Creating output file: $outputfile7\n";

unless (open(OUTFILE8, ">$outputfile8"))  {
	print "I can\'t write to $outputfile8\n";
	exit;
}
print "Creating output file: $outputfile8\n";

my @temp;
my @temp2;
my @temp3;
my @temp4;
my $nem_diffs;
my $nigra_diffs;
my $tonk_diffs;
my $nem_RAD_HiSeq_count;
my $nigra_RAD_HiSeq_count;
my $tonk_RAD_HiSeq_count;

my $nem_het_count=0;
my $nigra_het_count=0;
my $tonk_het_count=0;
my $nem_total_count=0;
my $nigra_total_count=0;
my $tonk_total_count=0;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne '#CHROM'){
		if(($temp[22] ne $temp[45])&&($temp[22] ne './.')&&($temp[45] ne './.')){ # comparing PM664 nem
			print OUTFILE3 $temp[1],"\t",$temp[22],"\t",$temp[45],"\n";
			$nem_diffs+=1;
		}
		if(($temp[27] ne $temp[46])&&($temp[27] ne './.')&&($temp[46] ne './.')){ # comparing PF660 nigra
			print OUTFILE4 $temp[1],"\t",$temp[27],"\t",$temp[46],"\n";
			$nigra_diffs+=1;
		}
		if(($temp[43] ne $temp[47])&&($temp[43] ne './.')&&($temp[47] ne './.')){ # comparing PM592 tonk
			print OUTFILE5 $temp[1],"\t",$temp[43],"\t",$temp[47],"\n";
			$tonk_diffs+=1;
		}
		if(($temp[22] ne './.')&&($temp[45] ne './.')){
			$nem_RAD_HiSeq_count+=1;
		}
		if(($temp[27] ne './.')&&($temp[46] ne './.')){
			$nigra_RAD_HiSeq_count+=1;
		}
		if(($temp[43] ne './.')&&($temp[47] ne './.')){
			$tonk_RAD_HiSeq_count+=1;
		}

		# now check for het sites on male X chrs
		@temp2=split('/',$temp[45]);
		@temp3=split('/',$temp[46]);
		@temp4=split('/',$temp[47]);
		if($temp2[0] ne $temp2[1]){
			print OUTFILE6 $temp[1],"\t",$temp[45],"\t",$temp[47],"\n";
			$nem_het_count+=1;
		}
		if($temp3[0] ne $temp3[1]){
			print OUTFILE7 $temp[1],"\t",$temp[45],"\t",$temp[47],"\n";
			$nigra_het_count+=1;
		}
		if($temp4[0] ne $temp4[1]){
			print OUTFILE8 $temp[1],"\t",$temp[45],"\t",$temp[47],"\n";
			$tonk_het_count+=1;
		}
		if($temp[45] ne './.'){
			$nem_total_count+=1;
		}
		if($temp[46] ne './.'){
			$nigra_total_count+=1;
		}
		if($temp[47] ne './.'){
			$tonk_total_count+=1;
		}
	}
	else{
		print OUTFILE3 "POS\t$temp[21]\t$temp[44]\n";
		print OUTFILE4 "POS\t$temp[26]\t$temp[45]\n";
		print OUTFILE5 "POS\t$temp[42]\t$temp[46]\n";
		print OUTFILE6 "POS\t$temp[44]\t$temp[46]\n";
	}
}

print "The number of neme hets is ",$nem_het_count," out of ",$nem_total_count," sites, and the proportion is ",$nem_het_count/$nem_total_count,"\n";
print "The number of nigr hets is ",$nigra_het_count," out of ",$nigra_total_count," sites, and the proportion is ",$nigra_het_count/$nigra_total_count,"\n";
print "The number of tonk hets is ",$tonk_het_count," out of ",$tonk_total_count," sites, and the proportion is ",$tonk_het_count/$tonk_total_count,"\n";

print "The number of neme diffs is ",$nem_diffs," out of ",$nem_RAD_HiSeq_count," sites, and the proportion is ",$nem_diffs/$nem_RAD_HiSeq_count,"\n";
print "The number of nigr diffs is ",$nigra_diffs," out of ",$nigra_RAD_HiSeq_count," sites, and the proportion is ",$nigra_diffs/$nigra_RAD_HiSeq_count,"\n";
print "The number of tonk diffs is ",$tonk_diffs," out of ",$tonk_RAD_HiSeq_count," sites, and the proportion is ",$tonk_diffs/$tonk_RAD_HiSeq_count,"\n";



close DATAINPUT;
close OUTFILE3;
close OUTFILE4;
close OUTFILE5;
close OUTFILE6;

```

