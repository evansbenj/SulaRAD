# Excluding sites with heterozygous genotypes on the X changes things somewhat, at least in the HiseqX data.

Well, after chasing down some of the reviewer concerns I now have identified a potentially serious problem.  When I looked at site patterns in the depth calls that had heterozygous calls on the non-PAR of the X in males, preliminary runs suggest that many more of the sites have ABBA patterns (tonk+nem gene flow) as compared to BABA patterns (nigra+nem gene flow). So these sites have high coverage reads in both males that are diverged from the outgroup and often from M. nigra.  A likely explanation is that they are from the Y.  This is spectacularly unfortunate because it means that all analyses that involve the X need to be repeated after excluding these sites.  This includes the phylogenetic estimation with the X, the RADseq analysis of polymorphism with Kai's models, the polymorphism tables, and the ABBABABA test on the X. This is surprising to me because I have previously done analyses where I excluded these sites.

To deal with this I need to identify sites that are heterozygous in the non-PAR of one or both males and delete them and do this separately in the HiSeq and RADseq data.

I have done this already in the HiSeqX data and listed them in this file:
```
/home/ben/2015_SulaRADtag/good_merged_samples/hets_on_one_or_both_male_Xs.txt
```
This file was created using a script called "Counts_chrX_het_sites_in_males.pl" which is below.  I ran it separately for tonk and nem on the hiSeqX diploid genotype file called:
```
/net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf.gz.tab
```

```
re#!/usr/bin/env perl
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
my $number_of_genotypes_not_in_PAR=0;
my %number_of_heterozygous_males_per_chrX_site;
my @number_of_het_sites_per_male;
my $number_of_sites_with_at_least_one_het_male=0;
my $number_of_sites_with_at_least_one_het_male_not_in_PAR=0;
my @position_of_het_sites;

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
			if($temp[1]>403495){
				$number_of_genotypes_not_in_PAR+=1;
			}

			$counter=0;
			foreach(@males){
				@bases=split('/',$temp[$_+$rhesus_and_begin[1]-1]);
				if($temp[$_+$rhesus_and_begin[1]-1] ne './.'){ # ignore sites with missing genotypes
					if(($bases[0] ne $bases[1])&&($bases[0] ne '*')&&($bases[1] ne '*')){
						$number_of_heterozygous_males_per_chrX_site{$temp[0].'_'.$temp[1]}+=1;
						$number_of_het_sites_per_male[$counter]+=1;
						push (@position_of_het_sites,$temp[1]);
					}
				$counter+=1;
				}	
			}
			if(defined($number_of_heterozygous_males_per_chrX_site{$temp[0].'_'.$temp[1]})){
				$number_of_sites_with_at_least_one_het_male+=1;
				if($temp[1]>403495){
					$number_of_sites_with_at_least_one_het_male_not_in_PAR+=1;
				}
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

print "The number_of_sites_with_at_least_one_het_male in the non-PAR is ",$number_of_sites_with_at_least_one_het_male_not_in_PAR,"\n";
print "The number of genotyped sites on the X  in the non-PAR is ",$number_of_genotypes_not_in_PAR,"\n";

foreach (@position_of_het_sites){
	print OUTFILE $_,"\n";
}

#foreach my $chr_and_position (keys %number_of_heterozygous_males_per_chrX_site ) {
#	print $chr_and_position,"\t",$number_of_heterozygous_males_per_chrX_site{$chr_and_position},"\n";
#}	

#$counter=0;
#foreach (@number_of_het_sites_per_male) {
#	if(defined($number_of_het_sites_per_male[$counter])){
#		print OUTFILE $namez[$counter],"\t",$number_of_het_sites_per_male[$counter],"\n";
#	}
#	else{
#		print OUTFILE $namez[$counter],"\t0\n";
#	}
#	$counter+=1;	
#}

print OUTFILE "The number_of_sites_with_at_least_one_het_male is ",$number_of_sites_with_at_least_one_het_male,"\n";
print OUTFILE "The number of genotyped sites on the X is ",$number_of_genotypes,"\n";

```


The worst part is that this will really take a lot of my time and a lot of time in general.

```
./Performs_ABBA_BABA_on_populations_onlychrX_haploid_alleles_haploid_outgroup_excludes_list_of_sites.pl /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf.all_with_one_highestdepth_allele.tab 000 3_4_1_2_3 H1nigra_H2tonk_H3nem_chrX_depth_excludeeithermalehetsites.jk H1nigra_H2tonk_H3nem_chrX_depth_excludeeithermalehetsites.stats hets_on_one_or_both_male_Xs.txt hets_on_one_or_both_male_Xs_patterns.out
```

In the RADseq data, this seems to be less of a problem. I found 317 sites with a heterozygous genotype in at least one male out of 230496 total genotypes. The position of these sites is in this file:
```
/home/ben/2015_SulaRADtag/good_merged_samples/fastq/positions_of_male_chrX_hets_in_RADseqdata.tab
```
And was generated in that directory using this command:
```
../Counts_chrX_het_sites_in_males.pl GenotypeVCFs_noBSQR_filtered_xDNA_only.vcf.gz.tab 1111100110000111100011100110010100000000 3_4 positions_of_male_chrX_hets_in_RADseqdata.tab
```

I am now making a new RADseq file with this command:
```
./Pulls_only_certain_lines_out.pl FINAL_RADseq_alldata_noBSQR_2016_haploiddepth_X.tab fastq/positions_of_male_chrX_hets_in_RADseqdata.tab FINAL_RADseq_alldata_noBSQR_2016_haploiddepth_X_nomaleXhets.tab
```

So the newly filtered RADseq file is:
```
FINAL_RADseq_alldata_noBSQR_2016_haploiddepth_X.tab
```


