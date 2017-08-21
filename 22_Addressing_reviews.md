# HiSeqX X:A polymorphism in M. nigra

Reviewer 2 was interested in knowning what the X:A ratio was in the female HiSeqX sample - very reasonable request.  I did this as follows:

```
Boot_from_tab_diverge_poly_2015.pl /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf_males_highest_depth_females_byvcf.tab 010 3_4_2 PF660_HiSeqX_nigra_poly_and_diverge.txt
```
in this directory:

```
/home/ben/2015_SulaRADtag/good_merged_samples
```

Another question concerned the number of heterozygous genotype calls in males on the X in the HiSeqX data.  This is also a very interesting request.  

This information is in these files:
```
-rw-rw-r-- 1 ben ben     970213 Nov 29  2016 /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/hetsites_nemtonk_chrX.vcf.gz.tab
-rw-rw-r-- 1 ben ben 1657155767 Nov 27  2016 /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf.gz.tab

```

Iassessed this using this script and commandline (with 011 and 110 settings for the sex to independently assay nem and tonk:

```
./Counts_chrX_het_sites_in_males.pl /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf.gz.tab 011 3_4 nemhetsites.out^C

```

Script:
```
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
my $counter=0;
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
my $number_of_genotypes=0;
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
			$number_of_genotypes+=1;
			$counter=0;
			foreach(@males){
				@bases=split('/',$temp[$_+$rhesus_and_begin[1]-1]);
				if($temp[$_+$rhesus_and_begin[1]-1] ne './.'){ # ignore sites with missing genotypes
					if(($bases[0] ne $bases[1])&&($bases[0] ne '*')&&($bases[1] ne '*')){
						$number_of_heterozygous_males_per_chrX_site{$temp[0].'_'.$temp[1]}+=1;
						$number_of_het_sites_per_male[$counter]+=1;
					}
				$counter+=1;
				}	
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
print "The number of genotyped sites on the X is ",$number_of_genotypes,"\n";



#foreach my $chr_and_position (keys %number_of_heterozygous_males_per_chrX_site ) {
#	print $chr_and_position,"\t",$number_of_heterozygous_males_per_chrX_site{$chr_and_position},"\n";
#}	

$counter=0;
foreach (@number_of_het_sites_per_male) {
	if(defined($number_of_het_sites_per_male[$counter])){
		print OUTFILE $namez[$counter],"\t",$number_of_het_sites_per_male[$counter],"\n";
	}
	else{
		print OUTFILE $namez[$counter],"\t0\n";
	}
	$counter+=1;	
}

print OUTFILE "The number_of_sites_with_at_least_one_het_male is ",$number_of_sites_with_at_least_one_het_male,"\n";
print OUTFILE "The number of genotyped sites on the X is ",$number_of_genotypes,"\n";



```
For tonkeaha PM592

The number_of_sites_with_at_least_one_het_male is 29208
The number of genotyped sites on the X is 58578113
Percentage: 0.05%
proportion: 0.0005
5 out of every 10,000 genotypes

For nemestrina PM664
The number_of_sites_with_at_least_one_het_male is 30422
The number of genotyped sites on the X is 58578113
Percentage: 0.05%
proportion: 0.0005
5 out of every 10,000 genotypes

For nigra female PF660 - heterozygosity is twice as high as the males
The number_of_sites_with_at_least_one_het_male is 61379
The number of genotyped sites on the X is 58578113
Percentage: 0.1%
proportion: 0.001
1 out of every 1,000 genotypes


