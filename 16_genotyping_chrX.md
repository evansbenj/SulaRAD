# Genotyping chrX is difficult

Background: my preliminary results suggest that there are some issues with genotyping chrX in males.  I have tried doing this with Haplotype caller using the default setting for ploidy (=2) and setting ploidy=1.  The former has problems because many heterozygous calls are made.  The latter has problems because the genotypes have a bias towards reference sites.

As a solution, Janet Kelso from the MPI suggested I genotype male chrX sites using the highest frequenty SNP.  This sounds like a great idea to me!  I plan to do this for the female as well so that each chrX is treated the same. I will write a script that outputs a tab delimited file from a diploid vcf file with genotyping based on the AD annotation (AD is allele depth).  For sites with an equal frequency of the ref and alt SNP, I will select the genotype randomly (and also keep track of how frequently this happens.

# First output chrX

I have some vcf files made at the NYGenome center here (on iqaluk):
```
/work/ben/2015_SulaRADtag/vcf-constitutional
```
```
-rw-rw-r-- 1 ben ben 6238245944 Sep  2 00:32 nemestrina-PM664.g.vcf.gz
-rw-rw-r-- 1 ben ben 6125321940 Sep  2 00:32 nigra-PM664.g.vcf.gz
-rw-rw-r-- 1 ben ben 6455292312 Sep  2 00:34 tonkeana-PM592.g.vcf.gz
```
Index the vcf files:
```
tabix -p vcf nemestrina-PM664.g.vcf.gz
```
Export the chrX:

```
~/tabix-0.2.6/tabix -h nemestrina-PM664.g.vcf.gz chrX > nemHiSeqchrX.vcf 
~/tabix-0.2.6/tabix -h nigra-PM664.g.vcf.gz chrX > nigraHiSeqchrX.vcf 
~/tabix-0.2.6/tabix -h tonkeana-PM592.g.vcf.gz chrX > tonkHiSeqchrX.vcf 
```

Convert from gvcf to vcf
```
/work/ben/2015_SulaRADtag/gvcftools-0.16/bin/break_blocks --ref /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --region-file /work/ben/2015_SulaRADtag/vcf-constitutional/target_interval_list_allchrs.bed < nemHiSeqchrX.vcf > nemHiSeqchrX.vcf.noblock.vcf

/work/ben/2015_SulaRADtag/gvcftools-0.16/bin/break_blocks --ref /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --region-file /work/ben/2015_SulaRADtag/vcf-constitutional/target_interval_list_allchrs.bed < tonkHiSeqchrX.vcf > tonkHiSeqchrX.vcf.noblock.vcf

/work/ben/2015_SulaRADtag/gvcftools-0.16/bin/break_blocks --ref /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --region-file /work/ben/2015_SulaRADtag/vcf-constitutional/target_interval_list_allchrs.bed < nigraHiSeqchrX.vcf > nigraHiSeqchrX.vcf.noblock.vcf
```
merge the files

```
bgzip nemHiSeqchrX.vcf.noblock.vcf
tabix -p vcf nemHiSeqchrX.vcf.noblock.vcf.gz
bgzip nigraHiSeqchrX.vcf.noblock.vcf
tabix -p vcf nigraHiSeqchrX.vcf.noblock.vcf.gz
bgzip tonkHiSeqchrX.vcf.noblock.vcf
tabix -p vcf tonkHiSeqchrX.vcf.noblock.vcf.gz

export PERL5LIB=/work/ben/vcftools/src/perl

/work/ben/vcftools/bin/vcf-merge nemHiSeqchrX.vcf.noblock.vcf.gz tonkHiSeqchrX.vcf.noblock.vcf.gz nigraHiSeqchrX.vcf.noblock.vcf.gz | bgzip -c > nem_tonk_nigra_alldiploid_chrX.vcf.gz
```

All of this appears to be bad because positions with no data are called as reference sites for some stupid reason.  So I am going to try again using GATK.
```
/usr/lib/jvm/java-1.8.0-openjdk.x86_64/bin/java -Xmx32G -jar /work/ben/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /work/ben/2015_SulaRADtag/HiSeqX/Project_MEL_11554_B01_CUS_WGS.2016-07-27/rheMac2_YM/rheMac2.fa --variant nemHiSeqchrX.g.vcf --variant tonkHiSeqchrX.g.vcf --variant nigraHiSeqchrX.g.vcf --includeNonVariantSites -o nem_tonk_nigra_HiSeq_combined_alldiploid_chrX.vcf
```
Here is a script that converts a vcf file to a majority rule tab delimited file, with haploid genotypes for all (Genotypes_chrX_based_on_allelic_depth.pl:

```
#!/usr/bin/env perl
use strict;
use warnings;
use List::Util 'max';
use List::Util qw(shuffle);

# This program reads in a vcf file then genotypes chrX sequences
# based on the AD (allelic depth) annotation.

# It takes as input a vcf file and outputs a tab delimited file

my $outputfile = "chrX.tab";

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile   $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";


my $inputfile = "./temp.vcf";
unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file, jackass.\n";
	exit;
}


my $y;
my $x;
my @columns=();
my @fields;
my $AD;
my $GT;
my $counter=0;
my @genotypes;
my $genotypez;
my @alleledepth;
my $max;
my @maxcounter=();
my $counter2=0;
my @altalleles=();

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@columns=split("\t",$line);
		if(substr($columns[0],0,1) ne '#'){ # this is not a comment
			@fields=split(":",$columns[8]);
			$counter=0;
			$AD=0;
			$GT=0;
			# first find out where the AD and GT columns are
			foreach(@fields){
				if($_ eq 'AD'){
					$AD=$counter;
				}
				elsif($_ eq 'GT'){
					$GT=$counter;
				}
				$counter+=1;
			}
			# now print out genotypes
			# first check if we have no data for any individuals
			$genotypez=();
			for ($y = 9 ; $y <= $#columns; $y++ ) {
				@genotypes=split(":",$columns[$y]);
				$genotypez=$genotypez.$genotypes[$GT];
			}
			print "genotypez ",$genotypez,"\n";
			if(
				(index($genotypez,'0') != -1)||
				(index($genotypez,'1') != -1)||
				(index($genotypez,'2') != -1)){
				# there is at least one genotype
				# if $AD==0 then all individuals are ref
				if($AD==0){ # this probably never happens
					print OUTFILE $columns[0],"\t",$columns[1],"\t",$columns[3];
					for ($y = 9 ; $y <= $#columns; $y++ ) {
						@genotypes=split(":",$columns[$y]);
						if($genotypes[$GT] eq '.\/.'){
							#print ref 
							print OUTFILE "\t\.\/";
						}
						elsif($genotypes[$GT] eq '0/0'){
							#print ref 
							print OUTFILE "\t".$columns[3]."\/";
						}
						else{
							print "Something is weird with the invariant genotypes\n";
						}
					}	
					print OUTFILE "\n";
				}
				else{ # This is probably what happens all the time
					print "AD $AD\n";
					print OUTFILE $columns[0],"\t",$columns[1],"\t",$columns[3];
					for ($y = 9 ; $y <= $#columns; $y++ ) {
						@alleledepth=();
						@genotypes=();
						@genotypes=split(":",$columns[$y]);
						@alleledepth=split(",",$genotypes[$AD]);
						@maxcounter=();
						$counter2=0;
						$max=0;
						$max=max @alleledepth;
						# now cycle through each allele depth to find highest and see if there is a tie
						foreach my $alleledepth (@alleledepth){
							if($alleledepth == $max){
								push(@maxcounter,$counter2);
							}
							$counter2+=1;
						}	
						@maxcounter = shuffle @maxcounter;
						if($genotypes[$GT] eq './.'){
							print OUTFILE "\t\.\/";
						}
						elsif($maxcounter[0] eq '0'){
							print OUTFILE "\t".$columns[3]."\/";
						}
						else{
							@altalleles = split(",",$columns[4]);
							if($altalleles[$maxcounter[0]-1] ne '*'){
								print OUTFILE "\t",$altalleles[$maxcounter[0]-1]."\/";
							}
							else{
								print OUTFILE "\t\.\/";
							}	
						}
					}
					print OUTFILE "\n";	
				}
			}
		}# endif
		elsif(substr($columns[0],0,6) eq '#CHROM'){ # print the first line
			print OUTFILE "#CHROM	POS	REF";
				for ($y = 9 ; $y <= $#columns; $y++ ) {
					print OUTFILE "\t",$columns[$y];
				}
				print OUTFILE "\n";			
		}
}# end while
close DATAINPUT;
close OUTFILE;

```
