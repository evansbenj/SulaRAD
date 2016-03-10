# Admixture analysis

I made a script to convert my tab delimited file to a geno file called "22_tab_to_geno.pl".  This file includes in the output only one SNP every 200 bp (so effectively one SNP per RADtag.  It also requires only 2 polymorphisms and that one must match the reference sequence. It also requires that at least half of the genotypes be not missing:

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;



#  This program reads in a tab delimited genotype file generated
#  by the perl program '17_adds_outgroup_to_lots_of_tab_files.pl'
#  or from vcftools vcf_to_tab
#  and generates a "geno" file as described here:
#  https://github.com/DReichLab/EIG/blob/master/CONVERTF/README
# specifically:
#EIGENSTRAT format: used by eigenstrat program
#  genotype file: see example.eigenstratgeno
#  snp file:      see example.snp (same as above)
#  indiv file:    see example.ind (same as above)
#Note that
#The genotype file contains 1 line per SNP.  
#  Each line contains 1 character per individual:
#  0 means zero copies of reference allele.
#  1 means one copy of reference allele.
#  2 means two copies of reference allele.
#  9 means missing data.

# commandline looks like this: 
# ./22_tab_to_geno.pl infile.tab 3_6 out.geno
# 3_6 refers to the column with the outgroup and the column where the ingroup start respectively


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
my $columns = $ARGV[1];
my $outputfile = $ARGV[2];

unless (open DATAINPUT, $inputfile) {
	print 'Can not find the input file.\n';
	exit;
}


unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}

my @temp;
my @temp1;
my @temp2;
my @temp3;

@temp=split('_',$ARGV[1]);
my $ref_seq=$temp[0]; # if in 3rd column, this should be a 3
my $start_of_ingroup=$temp[1]; # if in 6th column, this should be a 6
my $string;
my $string_no_9s;
my $bases;
my $x_uniq;
my $y_uniq;
my $y;
my $previous_chromosome="";
my $previous_position;
my $buffer = 200;

# Read in datainput file
while ( my $line = <DATAINPUT>) {
	chomp($line);
	#@temp=split(/[\/'\t']+/,$line);
	@temp=split('\t',$line);
	$string="";
	$string_no_9s="";
	$bases="";
	if($temp[0] ne '#CHROM'){
		if(length($temp[$ref_seq-1]) == 1){
			$bases=$bases.$temp[$ref_seq-1];
			# this will include only sites that match one ref base
			for ($y = ($start_of_ingroup-1); $y <= $#temp; $y++ ) {
				if($temp[$y] ne './.'){
					@temp1=split('/',$temp[$y]);
					if($temp1[0] eq $temp1[1]){
						if($temp1[0] eq $temp[$ref_seq-1]){
							# this is a homozygous ref genotype
							$string=$string.'2';
							$string_no_9s=$string_no_9s.'2';
							$bases=$bases.$temp1[0].$temp1[1];
						}
						else{
							# this is a homozygous alt seq
							$string=$string.'0';
							$string_no_9s=$string_no_9s.'0';
							$bases=$bases.$temp1[0].$temp1[1];							
						}	
					}
					elsif($temp1[0] ne $temp1[1]){
						# this is a heterozygous genotype
						$string=$string.'1';
						$string_no_9s=$string_no_9s.'1';
						$bases=$bases.$temp1[0].$temp1[1];						
					}	
				}
				else{
					$string=$string.'9';
				}
			}
			# test whether there is more than 2 SNPs in the data based on $bases
			# also test whether ingroup is homogenous based on $string_no_9s
			# and also require that at least 50% of the individuals have data
			# and also insert a requirement that the SNP not be on the 
			# same RADtag as the previous SNP
			@temp2=split('',$bases);
			@temp3=split('',$string_no_9s);
			$x_uniq = uniq @temp2;
			$y_uniq = uniq @temp3;
			if(($x_uniq == 2)&&($y_uniq != 1)&&(length($bases) >= length($string))){
				# if not, do not print the string
				# now check if the SNP is within $buffer bp of the previous SNP
				if($temp[0] eq $previous_chromosome){
					if($temp[1] > ($previous_position+$buffer)){
						print OUTFILE $string,"\n";
						$previous_position = $temp[1];
					}
				}
				else{
					print OUTFILE $string,"\n";
					$previous_chromosome = $temp[0];
					$previous_position = $temp[1];
				}	
			}
		}
	}
}	

close OUTFILE;
```

