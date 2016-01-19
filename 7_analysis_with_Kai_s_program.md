# Analysis with Kai's program

The first step to analyze each species/population with Kai's program is to generate an input file.  I wrote a script to do this (19_generates_input_for_Kai_program.pl):

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;



#  This program reads in a tab delimited genotype file generated
#  by the perl program "17_adds_outgroup_to_lots_of_tab_files.pl"
#  or from vcftools vcf_to_tab
#  and generates an input file for analysis with Kai's program

# to execute type 19_Generates_input_for_Kai_program.pl 
# /home/ben/2015_SulaRADtag/good_merged_samples/tab_relative_to_genez/recal_51000plus.vcf.gz.tab_with_baboon.tab_and_human.tab 
# 1111100110000111100011100110010100000000 
# 3_6_22_23_25_26 nigra_kai_input.txt  
# where 1111100110000111100011100110010100000000 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) or is not (0) female ,and 4_6_22_23_25_26 refers to (i) the column that contains the 
# outgroup nucleotide (4 in this case), (ii) the column number of the first individual in the ingroup 
# (6 in this case), and (iii) the sample number that contain the data from the individuals you want to 
# include (22, 23, 25, and 26 in this case), which are the four nigra samples itemized below.

# IMPORTANT: (i) and (ii) are columns beginning with 1 but (iii) is based on the individual samples such
# as enumerated below

# Notes: 
# ##### nigra_PM1000 is actually nigrescens_PM1000 #####

# ##### these samples have very low coverage (<10x): ####
# ##### nigrescens_PF654_sorted (7.33X) ####
# ##### maura_PM613_sorted (8.65X) ####
# ##### ochreata_PM596_sorted (9.02X)####
# ##### nigra_660_sorted (9.61X) ####
# ##### togeanus_PF549 (9.63X) ####

# Here is the order of the samples for the SulaRad project:


#	1	brunescens_PF707_stampy_sorted
#	2	hecki_PF643_stampy_sorted
#	3	hecki_PF644_stampy_sorted
#	4	hecki_PF648_stampy_sorted
#	5	hecki_PF651_stampy_sorted
#	6	hecki_PM639_stampy_sorted
#	7	hecki_PM645_stampy_sorted
#	8	maura_PF615_stampy_sorted
#	9	maura_PF713_stampy_sorted
#	10	maura_PM613_stampy_sorted
#	11	maura_PM614_stampy_sorted
#	12	maura_PM616_stampy_sorted
#	13	maura_PM618_stampy_sorted
#	14	nem_Gumgum_stampy_sorted
#	15	nem_Kedurang_stampy_sorted
#	16	nem_Malay_stampy_sorted
#	17	nem_Ngasang_stampy_sorted
#	18	nem_PM664_stampy_sorted
#	19	nem_PM665_stampy_sorted
#	20	nem_Sukai_male_stampy_sorted
#	21	nem_pagensis_stampy_sorted
#	22	nigra_PF1001_stampy_sorted
#	23	nigra_PF660_stampy_sorted
#	24	nigrescens_PM1000_stampy_sorted
#	25	nigra_PM1003_stampy_sorted
#	26	nigrescens_PF654_stampy_sorted
#	27	ochreata_PF625_stampy_sorted
#	28	ochreata_PM571_stampy_sorted
#	29	ochreata_PM596_stampy_sorted
#	30	togeanus_PF549_stampy_sorted
#	31	togeanus_PM545_stampy_sorted
#	32	tonk_PF515_stampy_sorted
#	33	tonk_PM561_stampy_sorted
#	34	tonk_PM565_stampy_sorted
#	35	tonk_PM566_stampy_sorted
#	36	tonk_PM567_stampy_sorted
#	37	tonk_PM582_stampy_sorted
#	38	tonk_PM584_stampy_sorted
#	39	tonk_PM592_stampy_sorted
#	40	tonk_PM602_stampy_sorted

# tonk
# 32_33_34_35_36_37_38_39_40 
# hecki
# 2_3_4_5_6_7  

# maura
# 8_9_10_11_12_13 

# nigra 
# 22_23_25 

# nigrescens 
# 24_26 

# ochreata 
# 27_28_29 

# togeanus 
# 30_31 

# brunn 
# 1 

# borneo 
# 18_19_20_21 

# sumatra 
# 15_17 

# pagensis 
# 21 

# malay 
# 16 


my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $input3 = $ARGV[2];
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

my @sexes = split("",$ARGV[1]);
my @whotoinclude = split("_",$ARGV[2]);

my $number_of_individuals_genotyped=($#whotoinclude - 1);

print "The number of individuals to assess is ",$number_of_individuals_genotyped,"\n";

my $number_of_female_individuals_genotyped;
my $y;


for ($y = 2 ; $y <= $#whotoinclude ; $y++ ) {
	if($sexes[$whotoinclude[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}	


print "This includes ",$number_of_female_individuals_genotyped," female(s)\n";


my @temp;
my @temp1;
my $string;
my $x;
my @unique;
my %ahash=();
my %xhash=();
my $type1=0;
my $type2=0;

my $w;


# initialize the ahash 
for ($y = 2 ; $y <= ($number_of_individuals_genotyped*2) ; $y+=2 ) {
	$ahash{$y}[0]=$y;
	for ($x = 1 ; $x <= ($y+1) ; $x++ ) {
		$ahash{$y}[$x]=0;
	}
}
# initialize the xhash 
for ($y = 2 ; $y <= ($number_of_individuals_genotyped*2 - ($number_of_individuals_genotyped - $number_of_female_individuals_genotyped)) ; $y++ ) {
	$xhash{$y}[0]=$y;
	for ($x = 1 ; $x <= ($y+1) ; $x++ ) {
		$xhash{$y}[$x]=0;
	}
}


# Read in datainput file
while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line);
	if($temp[0] ne '#CHROM'){
		if(($temp[0] ne "chrX")&&($temp[0] ne "chrY")&&(length($temp[2]) == 1)){
			# load the autosomal data
			$string=();
			for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
				# load the first allele
				if($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '.'){
					$w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
					$string=$string.$w;
				}	
				# now load the second allele
				if($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] ne '.'){
					$w = uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]];
					$string=$string.$w;
				}			
			}
			if(defined($string)){
				@temp1=split('',$string);
				if($#temp1 >0){					
					# count up the number of each class of nucleotide.
					# the first class is G or C
					# the second class is A or T
					# so if there is a G and a C segregating, for example, this is a non-polymorphic site
					$type1=0;
					$type2=0;

					for ($y = 0 ; $y <= $#temp1 ; $y++ ) {
						if((uc $temp1[$y] eq "G")||(uc $temp1[$y] eq "C")){
							$type1+=1;
						}
						elsif((uc $temp1[$y] eq "A")||(uc $temp1[$y] eq "T")){
							$type2+=1;
						}
						else{
							print "Problemo1! ",$line,"g\n";
						}
					}
					$ahash{($#temp1+1)}[$type1+1]+=1;
				}		
			}	
		}
		elsif(($temp[0] eq "chrX")&&(length($temp[2]) == 1)){
		 # load the chrX data			
			$string=();
			for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
				# load both alleles if the individual is a female
				if($sexes[$whotoinclude[$y+2]-1] eq "1"){
					# load the first allele
					if($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '.'){
						$w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
						$string=$string.$w;
					}	
					# now load the second allele
					if($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] ne '.'){
						$w = uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]];
						$string=$string.$w;
					}	
				}
				# load one allele if the individual is a male
				elsif($sexes[$whotoinclude[$y+2]-1] eq "0"){
					# load only the first allele
					if($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '.'){
						$w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
						$string=$string.$w;
					}	
				}
				else{
					print "Something is wrong with figuring out what sex each individual is X ",$sexes[$whotoinclude[$y+2]-1]," ",$whotoinclude[$y+2],"\n";
				}
			} # end of cycling through each individual for chrX	
			if(defined($string)){
				@temp1=split('',$string);
				if($#temp1 >0){					
					# count up the number of each class of nucleotide.
					# the first class is G or C
					# the second class is A or T
					# so if there is a G and a C segregating, for example, this is a non-polymorphic site
					$type1=0;
					$type2=0;
		
					for ($y = 0 ; $y <= $#temp1 ; $y++ ) {
						if((uc $temp1[$y] eq "G")||(uc $temp1[$y] eq "C")){
							$type1+=1;
						}
						elsif((uc $temp1[$y] eq "A")||(uc $temp1[$y] eq "T")){
							$type2+=1;
						}
						else{
							print "Problemo2!\n";
						}
					}
					$xhash{($#temp1+1)}[$type1+1]+=1;
				}		
			}	
		}
	} # endif to check for first line
} # end while

close DATAINPUT;


print OUTFILE "X\n";
print OUTFILE ($number_of_individuals_genotyped*2 - ($number_of_individuals_genotyped - $number_of_female_individuals_genotyped) -1),"\n";
for ($y = 2 ; $y <= ($number_of_individuals_genotyped*2 - ($number_of_individuals_genotyped - $number_of_female_individuals_genotyped)) ; $y++ ) {
	print OUTFILE $xhash{$y}[0],"\t";
	for ($x = 1 ; $x <= $y ; $x++ ) {
		print OUTFILE $xhash{$y}[$x],"\t";
	}
	print OUTFILE $xhash{$y}[$y+1],"\n";
}

print OUTFILE "A\n";
print OUTFILE $number_of_individuals_genotyped,"\n";
for ($y = 2 ; $y <= $number_of_individuals_genotyped*2 ; $y+=2 ) {
	print OUTFILE $ahash{$y}[0],"\t";
	for ($x = 1 ; $x <= $y ; $x++ ) {
		print OUTFILE $ahash{$y}[$x],"\t";
	}
	print OUTFILE $ahash{$y}[$y+1],"\n";
}


close OUTFILE;
print "Done with input file 1\n";






sub in_array {
my ($arr,$search_for) = @_;
foreach my $value (@$arr) {
	return 1 if $value eq $search_for;
}
 	return 0;
}


```

We then run Kai's program using different models 100-500 times to ensure convergence on the maximum likelihood parameter values.  His program is run like this:

```
#!/bin/bash
for i in {501..1000}
do
#sqsub -r 7d -q serial -o myout ../data_2a_step_xa_m1 control_file_3_epoch $i
    /work/ben/2013.09.16_kai_program/data_2a_step_xa_m1 control_file_3_epoch $i
done
```

and uses a control file that looks like this:
```
outputFile: ../3epoch.txt
dataFile: ../../tonkeana_kai_input.txt
K: 200
nstep: 2
maxTA: 0.5 0.5
tau: 0.01
useNrSimplex: 0
nlopt_alg: NLOPT_LN_NELDERMEAD
initThetaRange:	1e-10	0.5
initGammaRange:	-50	50
initLambdaRange:	0.01	100
initRhoRange:		0.01	100
seed:
thetaOnLn: 1
lambdaOnLn: 1
rhoOnLn: 1
setBound: 1
rftol: 1e-15
maxeval: 100000
maxtime: 600
imprftol: 1e-15
nnoimp: 3
maximp: 100

```


Results can be summarized with Kai's java program 'Parsefolder' as follows:

```
java ParseFolder /work/ben/new_kai_program_2013.10.13/2015_Sularad/tonkeana/epoch3_full 1e-15 true
```

Occasionally there are problems with the outputfile because of an '-inf' value for the likelihood.  This can be solved by going to the directory and typing this:

```
sed -i -n '/inf/!p' *
```
which deletes any lines containing `inf`.

