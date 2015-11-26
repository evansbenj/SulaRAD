# Subsetting the data and popgen stats.

## Popgen stats

I've polished up a script that calculates the important population statistics.  This script is quite flexible and it also does bootstrapping to provide confidence intervals:

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;


#  This program reads in a tab delimited genotype file generated
#  by the perl program "Gets_outgroup_sequence_from_axt_files.pl"
#  or from vcftools
#  and bootstraps TajD and polymorphism statistics by autosomal DNA 
#  and by xDNA by resampling bases
#  with replacement. 

#	Because TajD does not require an outgroup, the data analyzed will be 
# 	different from other analyzes that require divergence data (such as 
#	pi/D or S/D and also analyses that require outgroup information (such
#	as the analysis of the derived AFS).

#	This analysis will include only positions that have genotype data for
# 	ALL individuals.

# I am going to try to make this program compatible with files with one or multiple outgroup columns

# It will also accomodate data from multiple species and calculate the stats only from selected columns

# to execute type Boot_from_tab_diverge_poly.pl inputfile.tab 1111100110000111100011100110010100000000 
# 3_4_22_23_25_26 nigra_poly_and_diverge.txt  
# where 1111100110000111100011100110010100000000 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) or is not (0) female ,and 3_4_22_23_25_26 refers to (i) the column that contains the 
# outgroup nucleotide (3 in this case), (ii) the column number of the first individual in the ingroup 
# (4 in this case), and (iii) the sample number that contain the data from the individuals you want to 
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


# for example, with a tab file with only the baboon sequence in the 4th column, here is the input command:

# tonk
# Boot_from_tab_diverge_poly_2015.pl final_round2_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_32_33_34_35_36_37_38_39_40 tonk_poly_and_diverge.txt
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_32_33_34_35_36_37_38_39_40 tonk_poly_and_diverge.txt

# hecki
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_2_3_4_5_6_7 hecki_poly_and_diverge_baboon_recal.txt

# maura
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_8_9_10_11_12_13 maura_poly_and_diverge_baboon_recal.txt

# nigra 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_22_23_25 nigra_poly_and_diverge_baboon_recal.txt

# nigrescens 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_24_26 nigresc_poly_and_diverge_baboon_recal.txt

# ochreata 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_27_28_29 och_poly_and_diverge_baboon_recal.txt

# togeanus 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_30_31 tog_poly_and_diverge_baboon_recal.txt

# brunn 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_1 brun_poly_and_diverge_baboon_recal.txt

# borneo 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_18_19_20_21 borneo_poly_and_diverge_baboon_recal.txt

# sumatra 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_15_17 suma_poly_and_diverge_baboon_recal.txt

# pagensis 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_21 pagensis_poly_and_diverge_baboon_recal.txt

# malay 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon.tab 1111100110000111100011100110010100000000 4_5_16 malay_poly_and_diverge_baboon_recal.txt


# now, with a tab file with human and baboon sequence in the 4th column, here is the input command with human as outgroup:

# tonk
# Boot_from_tab_diverge_poly_2015.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 4_6_32_33_34_35_36_37_38_39_40 tonk_poly_and_diverg_allsitese_humanout_1000boot.txt

# hecki
# Boot_from_tab_diverge_poly_2015.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 4_6_2_3_4_5_6_7 hecki_poly_and_diverg_allsitese_humanout_1000boot.txt

# maura
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 4_6_8_9_10_11_12_13 maura_poly_and_diverge_baboon_recal.txt

# nigra 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 4_6_22_23_25 nigra_poly_and_diverge_baboon_recal.txt

# nigrescens 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 4_6_24_26 nigresc_poly_and_diverge_baboon_recal.txt

# ochreata 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 4_6_27_28_29 och_poly_and_diverge_baboon_recal.txt

# togeanus 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 4_6_30_31 tog_poly_and_diverge_baboon_recal.txt

# brunn 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 4_6_1 brun_poly_and_diverge_baboon_recal.txt

# borneo 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 4_6_18_19_20_21 borneo_poly_and_diverge_baboon_recal.txt

# sumatra 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 4_6_15_17 suma_poly_and_diverge_baboon_recal.txt

# pagensis 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 4_6_21 pagensis_poly_and_diverge_baboon_recal.txt

# malay 
# Boot_from_tab_diverge_poly_2015.pl nonrecal_final_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 4_6_16 malay_poly_and_diverge_baboon_recal.txt


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


my @temp;
my @temp1;
my $previous= 0;
my $aDNA_segregating_sites=0;
my $aDNA_sites=0;
my $xDNA_segregating_sites=0;
my $xDNA_sites=0;
my $yDNA_segregating_sites=0;
my $yDNA_sites=0;
my $aDNA_divergence=0;
my $xDNA_divergence=0;
my $yDNA_divergence=0;
my $string;
my $m;
my $n;
my $number_of_bootstraps=10;
my $lower=int($number_of_bootstraps*0.025);
my $upper=int($number_of_bootstraps*0.975);
my $JC_divergence_aDNA;
my $JC_divergence_xDNA;
my $JC_divergence_yDNA;

print "bootlower ", $lower," bootupper ",$upper,"\n";

my $w;
my $y;
my $x;
my @unique;

my $number_of_individuals_genotyped=($#whotoinclude - 1);

print "The number of individuals to assess is ",$number_of_individuals_genotyped,"\n";

my $number_of_female_individuals_genotyped;
for ($y = 2 ; $y <= $#whotoinclude ; $y++ ) {
	if($sexes[$whotoinclude[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}	

print "This includes ",$number_of_female_individuals_genotyped," female(s)\n";

my $asum=0;
my $xsum=0;
my $ysum=0;
my $asum_baboons=0;
my $xsum_baboons=0;
my $ysum_baboons=0;

my @d_aDNA=();
my @d_xDNA=();
my @d_yDNA=();
my @pi_aDNA=();
my @pi_xDNA=();
my @pi_yDNA=();
my @S_aDNA=();
my @S_xDNA=();
my @S_yDNA=();
my $pi_counter=0;
my $diff=0;
my @slice;
my $x_uniq=0;
my $y_uniq=0;

# TAJIMA'S D VARIABLES based on http://en.wikipedia.org/wiki/Tajima's_D

my $n_aDNA=(2*$number_of_individuals_genotyped);  # this is the number of alleles required for aDNA
my $n_xDNA=(2*$number_of_female_individuals_genotyped)+($number_of_individuals_genotyped-$number_of_female_individuals_genotyped);  # this is the number of alleles required for xDNA
my $n_yDNA=$number_of_individuals_genotyped-$number_of_female_individuals_genotyped;  # this is the number of alleles required for yDNA

print "Num ausomomal alleles: ",$n_aDNA,"\n";
print "Num chrX alleles: ",$n_xDNA,"\n";
print "Num chrY alleles: ",$n_yDNA,"\n";
my $expected_number_of_adna_pairwise_comparisons=0;
my $expected_number_of_xdna_pairwise_comparisons=0;
my $expected_number_of_ydna_pairwise_comparisons=0;

for ($y = 1 ; $y < $n_aDNA ; $y++ ) {
	$expected_number_of_adna_pairwise_comparisons+=$y;	
}

for ($y = 1 ; $y < $n_xDNA ; $y++ ) {
	$expected_number_of_xdna_pairwise_comparisons+=$y;	
}

for ($y = 1 ; $y < $n_yDNA ; $y++ ) {
	$expected_number_of_ydna_pairwise_comparisons+=$y;	
}

my $TajD_aDNA;
my $TajD_xDNA;
my $TajD_yDNA;


# Calculate some constants for aDNA that depend only on the number of samples
my $a1_obs_aDNA;
	for ($y = 1 ; $y < $n_aDNA ; $y++ ) {
		$a1_obs_aDNA+= 1/$y;
	}
	
my $a2_obs_aDNA;
	for ($y = 1 ; $y < $n_aDNA ; $y++ ) {
		$a2_obs_aDNA+= 1/($y**2);
	}
my $b1_obs_aDNA = ($n_aDNA+1)/(3*($n_aDNA-1));
my $b2_obs_aDNA = (2*($n_aDNA**2+$n_aDNA+3))/(9*$n_aDNA*($n_aDNA-1));
my $c1_obs_aDNA = $b1_obs_aDNA - (1/$a1_obs_aDNA);
my $c2_obs_aDNA = $b2_obs_aDNA - ($n_aDNA+2)/($a1_obs_aDNA*$n_aDNA) + $a2_obs_aDNA/($a1_obs_aDNA**2);
my $e1_obs_aDNA = $c1_obs_aDNA/$a1_obs_aDNA;
my $e2_obs_aDNA = $c2_obs_aDNA/(($a1_obs_aDNA**2)+$a2_obs_aDNA);


# Calculate some constants for xDNA that depend only on the number of samples
my $a1_obs_xDNA;
	for ($y = 1 ; $y < $n_xDNA ; $y++ ) {
		$a1_obs_xDNA+= 1/$y;
	}
my $a2_obs_xDNA;
	for ($y = 1 ; $y < $n_xDNA ; $y++ ) {
		$a2_obs_xDNA+= 1/($y**2);
	}
my $b1_obs_xDNA = ($n_xDNA+1)/(3*($n_xDNA-1));
my $b2_obs_xDNA = (2*($n_xDNA**2+$n_xDNA+3))/(9*$n_xDNA*($n_xDNA-1));
my $c1_obs_xDNA = $b1_obs_xDNA - (1/$a1_obs_xDNA);
my $c2_obs_xDNA = $b2_obs_xDNA - ($n_xDNA+2)/($a1_obs_xDNA*$n_xDNA) + $a2_obs_xDNA/($a1_obs_xDNA**2);
my $e1_obs_xDNA = $c1_obs_xDNA/$a1_obs_xDNA;
my $e2_obs_xDNA = $c2_obs_xDNA/(($a1_obs_xDNA**2)+$a2_obs_xDNA);

# Calculate some constants for yDNA that depend only on the number of samples
my $a1_obs_yDNA;
my $a2_obs_yDNA;
my $b1_obs_yDNA;
my $b2_obs_yDNA;
my $c1_obs_yDNA;
my $c2_obs_yDNA;
my $e1_obs_yDNA;
my $e2_obs_yDNA;

if($n_yDNA>1){
	for ($y = 1 ; $y < $n_yDNA ; $y++ ) {
		$a1_obs_yDNA+= 1/$y;
	}

	for ($y = 1 ; $y < $n_yDNA ; $y++ ) {
		$a2_obs_yDNA+= 1/($y**2);
	}
	$b1_obs_yDNA = ($n_yDNA+1)/(3*($n_yDNA-1));
	$b2_obs_yDNA = (2*($n_yDNA**2+$n_yDNA+3))/(9*$n_yDNA*($n_yDNA-1));
	$c1_obs_yDNA = $b1_obs_yDNA - (1/$a1_obs_yDNA);
	$c2_obs_yDNA = $b2_obs_yDNA - ($n_yDNA+2)/($a1_obs_yDNA*$n_yDNA) + $a2_obs_yDNA/($a1_obs_yDNA**2);
	$e1_obs_yDNA = $c1_obs_yDNA/$a1_obs_yDNA;
	$e2_obs_yDNA = $c2_obs_yDNA/(($a1_obs_yDNA**2)+$a2_obs_yDNA);
}

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line);
	if($temp[0] ne '#CHROM'){
		#base the genomic location on the outgroup	
		# this could be changed later to always rely on the most closely related ingroup
		if(($temp[0] ne "chrX")&&($temp[0] ne "chrY")){
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
					if(($#temp1 == ($n_aDNA-1))&&((uc $temp[$whotoinclude[0]-1] eq "A")||(uc $temp[$whotoinclude[0]-1] eq "C")||(uc $temp[$whotoinclude[0]-1] eq "T")||(uc $temp[$whotoinclude[0]-1] eq "G"))){
						$x_uniq = uniq @temp1;
						if($x_uniq == 1){
									push(@pi_aDNA,0);
									push(@S_aDNA,0);
						}
						elsif($x_uniq == 2){
								$diff=0;
								for ($y = 0 ; $y < $#temp1 ; $y++ ) {
									for ($x = ($y+1) ; $x <= $#temp1 ; $x++ ) {
										if($temp1[$y] ne $temp1[$x]){
											$diff+=1;
										}
										$pi_counter+=1;
									}
								}
								if($pi_counter ne $expected_number_of_adna_pairwise_comparisons){
									print "problem with number of pairwise comparisons aDNA\n";
								}
								push(@pi_aDNA,$diff/$pi_counter);
								$diff=0;
								$aDNA_segregating_sites+=1;	
								push(@S_aDNA,1);					
								$pi_counter=0;	
						}
						if(($x_uniq == 1)||($x_uniq == 2)){
								$aDNA_sites+=1;
								if(uc $temp[$whotoinclude[0]-1] ne uc $temp1[0]){
									#print $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\n";
									#print uc $temp[$whotoinclude[0]-1]," hey ", uc $temp1[0],"\n";
									$aDNA_divergence+=1;
									push(@d_aDNA,1);
								}
								else{
									push(@d_aDNA,0);
								}
						}	
						if($aDNA_sites != ($#S_aDNA+1)){
							print $aDNA_sites,"\t",$#S_aDNA+1,"\t",$line,"\n";
						}
					}
				}	
		}
		elsif($temp[0] eq "chrX"){
					$string=();
					# for chrX, load both female alleles but only one male allele
					# for each column, we need to check if the individual is a female
					# first the female
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
				if(($#temp1 == ($n_xDNA-1))&&((uc $temp[$whotoinclude[0]-1] eq "A")||(uc $temp[$whotoinclude[0]-1] eq "C")||(uc $temp[$whotoinclude[0]-1] eq "T")||(uc $temp[$whotoinclude[0]-1] eq "G"))){
					$x_uniq = uniq @temp1;		
					if($x_uniq == 1){
						push(@pi_xDNA,0);
						push(@S_xDNA,0);
					}
					elsif($x_uniq == 2){
						$diff=0;
						for ($y = 0 ; $y < $#temp1 ; $y++ ) {
							for ($x = ($y+1) ; $x <= $#temp1 ; $x++ ) {
								if($temp1[$y] ne $temp1[$x]){
									$diff+=1;
								}
								$pi_counter+=1;
							}
						}
						if($pi_counter ne $expected_number_of_xdna_pairwise_comparisons){
							print "problem with number of pairwise comparisons XDNA\n";
						}
						push(@pi_xDNA,$diff/$pi_counter);
						$diff=0;
						$xDNA_segregating_sites+=1;	
						push(@S_xDNA,1);					
						$pi_counter=0;	
					}
					if(($x_uniq == 1)||($x_uniq == 2)){
						$xDNA_sites+=1;
						if(uc $temp[$whotoinclude[0]-1] ne uc $temp1[0]){
							$xDNA_divergence+=1;
							push(@d_xDNA,1);
						}
						else{
							push(@d_xDNA,0);
						}
					}	
					if($xDNA_sites != ($#S_xDNA+1)){
						print $xDNA_sites,"\t",$#S_xDNA+1,"\t",$line,"\n";
					}
				}
			}
		} # endifelse to check for aDNA, chrX	
		elsif($temp[0] eq "chrY"){
			$string=();
			# for chrY, load only the male allele
			# for each column, we need to check if the individual is a male
			for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
				# load one allele if the individual is a male
				if($sexes[$whotoinclude[$y+2]-1] eq "0"){
					# load only the first allele
					if($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '.'){
						$w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
						$string=$string.$w;
					}	
				}
				elsif($sexes[$whotoinclude[$y+2]-1] ne "1"){
					print "Something is wrong with figuring out what sex each individual is Y"
				}
			} # end of cycling through each individual for chrY	
			if(defined($string)){
				@temp1=split('',$string);			
				if(($#temp1 == ($n_yDNA-1))&&((uc $temp[$whotoinclude[0]-1] eq "A")||(uc $temp[$whotoinclude[0]-1] eq "C")||(uc $temp[$whotoinclude[0]-1] eq "T")||(uc $temp[$whotoinclude[0]-1] eq "G"))){
					$x_uniq = uniq @temp1;		
					if($x_uniq == 1){
						push(@pi_yDNA,0);
						push(@S_yDNA,0);
					}
					elsif($x_uniq == 2){
						$diff=0;
						for ($y = 0 ; $y < $#temp1 ; $y++ ) {
							for ($x = ($y+1) ; $x <= $#temp1 ; $x++ ) {
								if($temp1[$y] ne $temp1[$x]){
									$diff+=1;
								}
								$pi_counter+=1;
							}
						}
						if($pi_counter ne $expected_number_of_ydna_pairwise_comparisons){
							print "problem with number of pairwise comparisons yDNA\n";
						}
						push(@pi_yDNA,$diff/$pi_counter);
						$diff=0;
						$yDNA_segregating_sites+=1;	
						push(@S_yDNA,1);					
						$pi_counter=0;	
					}
					if(($x_uniq == 1)||($x_uniq == 2)){
						$yDNA_sites+=1;
						if(uc $temp[$whotoinclude[0]-1] ne uc $temp1[0]){
							$yDNA_divergence+=1;
							push(@d_yDNA,1);
						}
						else{
							push(@d_yDNA,0);
						}
					}	
					if($yDNA_sites != ($#S_yDNA+1)){
						print $yDNA_sites,"\t",$#S_yDNA+1,"\t",$line,"\n";
					}
				}
			}
		} # endifelse to check for aDNA, chrX, chrY	
	} # endif to check for first line
	elsif($temp[0] eq '#CHROM'){
		for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
			print "Individual ",$temp[$whotoinclude[$y+2]+$whotoinclude[1]-2]," is a ";
			if($sexes[$whotoinclude[$y+2]-1] == 1){
				print "female\n";
			} 
			elsif($sexes[$whotoinclude[$y+2]-1] == 0){
				print "male\n";
			} 
		}
	}	
} # end while


for ($y = 0 ; $y <= $#pi_aDNA ; $y++ ) {
	$asum+=$pi_aDNA[$y];
}
for ($y = 0 ; $y <= $#pi_xDNA ; $y++ ) {
	$xsum+=$pi_xDNA[$y];
}
for ($y = 0 ; $y <= $#pi_yDNA ; $y++ ) {
	$ysum+=$pi_yDNA[$y];
}


print $aDNA_sites," ",($#pi_aDNA+1)," ",($#S_aDNA+1),"\n";
print $xDNA_sites," ",($#pi_xDNA+1)," ",($#S_xDNA+1),"\n";
print $yDNA_sites," ",($#pi_yDNA+1)," ",($#S_yDNA+1),"\n";

my @aDNA_bootstrapped_indexes=();
my @xDNA_bootstrapped_indexes=();
# no bootstrap for the y - use coalescent simulations

my @S_aDNA_boot;
my @pi_aDNA_boot;
my @pi_aDNA_persite_boot;
my @d_aDNA_boot;

my @S_xDNA_boot;
my @pi_xDNA_boot;
my @pi_xDNA_persite_boot;
my @d_xDNA_boot;

my @aDNA_TajD_boot;
my @xDNA_TajD_boot;

my @JC_AP_divergence_aDNA_boot;
my @JC_AP_divergence_xDNA_boot;

my @thetapi_over_divergence_aDNA;
my @thetapi_over_divergence_xDNA;


for ($m = 0 ; $m < $number_of_bootstraps ; $m++ ) {
	# generate an array with bootstrapped indexes
	# first aDNA
	print "bootstrap ",$m,"\n";
	
	for ($n = 0 ; $n < $aDNA_sites ; $n++ ) {
		push(@aDNA_bootstrapped_indexes,int(rand($aDNA_sites)));
	} # end $n

	for ($n = 0 ; $n < $aDNA_sites ; $n++ ) {
		# calculate segregating sites and pi and d for bootstrap replicates
		$S_aDNA_boot[$m]+=$S_aDNA[$aDNA_bootstrapped_indexes[$n]];
		$pi_aDNA_boot[$m]+=$pi_aDNA[$aDNA_bootstrapped_indexes[$n]];
		$d_aDNA_boot[$m]+=$d_aDNA[$aDNA_bootstrapped_indexes[$n]];
	} # end $n
	
	
	push(@pi_aDNA_persite_boot,$pi_aDNA_boot[$m]/$aDNA_sites);

	# apply JC correction to divergence
	$JC_AP_divergence_aDNA_boot[$m]  = (-3/4)*log(1-(4/3)*($d_aDNA_boot[$m]/$aDNA_sites));
	# apply AP correction to divergence
	$JC_AP_divergence_aDNA_boot[$m]  = $JC_AP_divergence_aDNA_boot[$m] - ($pi_aDNA_boot[$m]/$aDNA_sites);

	# calculate pi/divergence for bootstrapped data
	push(@thetapi_over_divergence_aDNA,($pi_aDNA_boot[$m]/$aDNA_sites)/$JC_AP_divergence_aDNA_boot[$m]);

	# calculate TajD for bootstrapped data
	push(@aDNA_TajD_boot,($pi_aDNA_boot[$m] - ($S_aDNA_boot[$m]/$a1_obs_aDNA))/((($e1_obs_aDNA*$S_aDNA_boot[$m])+($e2_obs_aDNA*$S_aDNA_boot[$m]*($S_aDNA_boot[$m]-1)))**0.5));
   

	# now do xDNA
	for ($n = 0 ; $n < $xDNA_sites ; $n++ ) {
		push(@xDNA_bootstrapped_indexes,int(rand($xDNA_sites)));
	} # end $n

	for ($n = 0 ; $n < $xDNA_sites ; $n++ ) {
		# calculate segregating sites and pi for bootstrap replicates
		$S_xDNA_boot[$m]+=$S_xDNA[$xDNA_bootstrapped_indexes[$n]];
		$pi_xDNA_boot[$m]+=$pi_xDNA[$xDNA_bootstrapped_indexes[$n]];
		$d_xDNA_boot[$m]+=$d_xDNA[$xDNA_bootstrapped_indexes[$n]];
	} # end $n


	push(@pi_xDNA_persite_boot,$pi_xDNA_boot[$m]/$xDNA_sites);

	# apply JC correction to divergence
	$JC_AP_divergence_xDNA_boot[$m]  = (-3/4)*log(1-(4/3)*($d_xDNA_boot[$m]/$xDNA_sites));
	# apply AP correction to divergence
	$JC_AP_divergence_xDNA_boot[$m]  = $JC_AP_divergence_xDNA_boot[$m] - ($pi_xDNA_boot[$m]/$xDNA_sites);

	# calculate pi/divergence for bootstrapped data
	push(@thetapi_over_divergence_xDNA,($pi_xDNA_boot[$m]/$xDNA_sites)/$JC_AP_divergence_xDNA_boot[$m]);
	# calculate TajD for bootstrapped data
	push(@xDNA_TajD_boot,($pi_xDNA_boot[$m] - ($S_xDNA_boot[$m]/$a1_obs_xDNA))/((($e1_obs_xDNA*$S_xDNA_boot[$m])+($e2_obs_xDNA*$S_xDNA_boot[$m]*($S_xDNA_boot[$m]-1)))**0.5));
	# reset the indexes for next bootstrap
	@aDNA_bootstrapped_indexes=();
	@xDNA_bootstrapped_indexes=();
} # end $m bootstraps


#print "aDNA_boot @aDNA_TajD_boot \n";
#print "xDNA_boot @xDNA_TajD_boot \n";


# now get 95% CIs
@d_aDNA_boot = sort { $a <=> $b } @d_aDNA_boot;
@d_xDNA_boot = sort { $a <=> $b } @d_xDNA_boot;

@S_aDNA_boot = sort { $a <=> $b } @S_aDNA_boot;
@S_xDNA_boot = sort { $a <=> $b } @S_xDNA_boot;

@pi_aDNA_persite_boot = sort { $a <=> $b } @pi_aDNA_persite_boot;
@pi_xDNA_persite_boot = sort { $a <=> $b } @pi_xDNA_persite_boot;

@thetapi_over_divergence_aDNA = sort { $a <=> $b } @thetapi_over_divergence_aDNA;
@thetapi_over_divergence_xDNA = sort { $a <=> $b } @thetapi_over_divergence_xDNA;

@aDNA_TajD_boot = sort { $a <=> $b } @aDNA_TajD_boot;
@xDNA_TajD_boot = sort { $a <=> $b } @xDNA_TajD_boot;

@JC_AP_divergence_aDNA_boot  = sort { $a <=> $b } @JC_AP_divergence_aDNA_boot;
@JC_AP_divergence_xDNA_boot  = sort { $a <=> $b } @JC_AP_divergence_xDNA_boot;




print OUTFILE "aDNA\n";
print OUTFILE "# alleles\t",$n_aDNA,"\n";
print OUTFILE "# Sites\t",$aDNA_sites,"\n";
print OUTFILE "S\t",$aDNA_segregating_sites," (",$S_aDNA_boot[$lower]," - ",$S_aDNA_boot[$upper],")\n";
print OUTFILE "thetaW\t",sprintf("%.5f",$aDNA_segregating_sites/$a1_obs_aDNA/$aDNA_sites)," (",sprintf("%.5f",$S_aDNA_boot[$lower]/$a1_obs_aDNA/$aDNA_sites)," - ",sprintf("%.5f",$S_aDNA_boot[$upper]/$a1_obs_aDNA/$aDNA_sites),")\n";
print OUTFILE "pi\t",sprintf("%.5f",$asum/($#pi_aDNA+1))," (",sprintf("%.5f",$pi_aDNA_persite_boot[$lower])," - ",sprintf("%.5f",$pi_aDNA_persite_boot[$upper]),")\n";
print OUTFILE "d\t",sprintf("%.5f",$aDNA_divergence/$aDNA_sites)," (",sprintf("%.5f",$d_aDNA_boot[$lower]/$aDNA_sites)," - ",sprintf("%.5f",$d_aDNA_boot[$upper]/$aDNA_sites),")\n";
# apply JC correction
$JC_divergence_aDNA = (-3/4)*log(1-(4/3)*($aDNA_divergence/$aDNA_sites));
#print OUTFILE "d_JC\t",sprintf("%.5f",$JC_divergence_aDNA),"\n";
# apply correction for ancestral polymorphism
$JC_divergence_aDNA = $JC_divergence_aDNA- $asum/($#pi_aDNA+1);
print OUTFILE "d_JC_AP\t",sprintf("%.5f",$JC_divergence_aDNA)," (",sprintf("%.5f",$JC_AP_divergence_aDNA_boot[$lower])," - ",sprintf("%.5f",$JC_AP_divergence_aDNA_boot[$upper]),")\n";
print OUTFILE "pi/d_JC_AP\t",sprintf("%.3f",($asum/($#pi_aDNA+1))/$JC_divergence_aDNA)," (",sprintf("%.5f",$thetapi_over_divergence_aDNA[$lower])," - ",sprintf("%.5f",$thetapi_over_divergence_aDNA[$upper]),")\n";

# TajD aDNA
if($n_aDNA > 2){
	$TajD_aDNA = ($asum - ($aDNA_segregating_sites/$a1_obs_aDNA))/((($e1_obs_aDNA*$aDNA_segregating_sites)+($e2_obs_aDNA*$aDNA_segregating_sites*($aDNA_segregating_sites-1)))**0.5);
	print OUTFILE "TajD\t",sprintf("%.3f",$TajD_aDNA)," (",sprintf("%.3f",$aDNA_TajD_boot[$lower])," - ",sprintf("%.3f",$aDNA_TajD_boot[$upper]),")\n";
}
else{
 	$TajD_aDNA = "undefined";
 	print OUTFILE "TajD\tundefined\n";
}



print OUTFILE "xDNA\n";
print OUTFILE "# alleles\t",$n_xDNA,"\n";
print OUTFILE "# Sites\t",$xDNA_sites,"\n";
print OUTFILE "S\t",$xDNA_segregating_sites," (",$S_xDNA_boot[$lower]," - ",$S_xDNA_boot[$upper],")\n";
print OUTFILE "thetaW\t",sprintf("%.5f",$xDNA_segregating_sites/$a1_obs_xDNA/$xDNA_sites)," (",sprintf("%.5f",$S_xDNA_boot[$lower]/$a1_obs_xDNA/$xDNA_sites)," - ",sprintf("%.5f",$S_xDNA_boot[$upper]/$a1_obs_xDNA/$xDNA_sites),")\n";
print OUTFILE "pi\t",sprintf("%.5f",$xsum/($#pi_xDNA+1))," (",sprintf("%.5f",$pi_xDNA_persite_boot[$lower])," - ",sprintf("%.5f",$pi_xDNA_persite_boot[$upper]),")\n";
print OUTFILE "d\t",sprintf("%.5f",$xDNA_divergence/$xDNA_sites)," (",sprintf("%.5f",$d_xDNA_boot[$lower]/$xDNA_sites)," - ",sprintf("%.5f",$d_xDNA_boot[$upper]/$xDNA_sites),")\n";
# apply JC correction
$JC_divergence_xDNA = (-3/4)*log(1-(4/3)*($xDNA_divergence/$xDNA_sites));
#print OUTFILE "d_JC\t",sprintf("%.5f",$JC_divergence_xDNA),"\n";
# apply correction for ancestral polymorphism
$JC_divergence_xDNA = $JC_divergence_xDNA- $xsum/($#pi_xDNA+1);
print OUTFILE "d_JC_AP\t",sprintf("%.5f",$JC_divergence_xDNA)," (",sprintf("%.5f",$JC_AP_divergence_xDNA_boot[$lower])," - ",sprintf("%.5f",$JC_AP_divergence_xDNA_boot[$upper]),")\n";
print OUTFILE "pi/d_JC_AP\t",sprintf("%.3f",($xsum/($#pi_xDNA+1))/$JC_divergence_xDNA)," (",sprintf("%.5f",$thetapi_over_divergence_xDNA[$lower])," - ",sprintf("%.5f",$thetapi_over_divergence_xDNA[$upper]),")\n";
# TajD xDNA
if($n_xDNA > 2){
	$TajD_xDNA = ($xsum - ($xDNA_segregating_sites/$a1_obs_xDNA))/((($e1_obs_xDNA*$xDNA_segregating_sites)+($e2_obs_xDNA*$xDNA_segregating_sites*($xDNA_segregating_sites-1)))**0.5);
	print OUTFILE "TajD\t",sprintf("%.3f",$TajD_xDNA)," (",sprintf("%.3f",$xDNA_TajD_boot[$lower])," - ",sprintf("%.3f",$xDNA_TajD_boot[$upper]),")\n";
}
else{
	$TajD_xDNA = "undefined";
 	print OUTFILE "TajD\tundefined\n";
}

print OUTFILE "\n";
if(($number_of_individuals_genotyped-$number_of_female_individuals_genotyped)>0){
print OUTFILE "yDNA\n";
print OUTFILE "# alleles\t",$n_yDNA,"\n";
print OUTFILE "# Sites\t",$yDNA_sites,"\n";
print OUTFILE "S\t",$yDNA_segregating_sites,"\n";
print OUTFILE "thetaW\t",sprintf("%.5f",$yDNA_segregating_sites/$a1_obs_yDNA/$yDNA_sites),"\n";
print OUTFILE "pi\t",sprintf("%.5f",$ysum/($#pi_yDNA+1)),"\n";
print OUTFILE "d\t",sprintf("%.5f",$yDNA_divergence/$yDNA_sites),"\n";
# apply JC correction
$JC_divergence_yDNA = (-3/4)*log(1-(4/3)*($yDNA_divergence/$yDNA_sites));
#print OUTFILE "d_JC\t",sprintf("%.5f",$JC_divergence_yDNA),"\n";
# apply correction for ancestral polymorphism
$JC_divergence_yDNA = $JC_divergence_yDNA- $ysum/($#pi_yDNA+1);
print OUTFILE "d_JC_AP\t",sprintf("%.5f",$JC_divergence_yDNA),"\n";
print OUTFILE "pi/d_JC_AP\t",sprintf("%.3f",($ysum/($#pi_yDNA+1))/$JC_divergence_yDNA),"\n";
}

# TajD yDNA
if($n_yDNA>2){
	$TajD_yDNA = ($ysum - ($yDNA_segregating_sites/$a1_obs_yDNA))/((($e1_obs_yDNA*$yDNA_segregating_sites)+($e2_obs_yDNA*$yDNA_segregating_sites*($yDNA_segregating_sites-1)))**0.5);
	print OUTFILE "TajD\t",sprintf("%.3f",$TajD_yDNA),"\n";
}
else{
	$TajD_yDNA = "undefined";
	print OUTFILE "TajD\t",$TajD_yDNA,"\n";
}
print OUTFILE "\n";
print OUTFILE "pi_X/d_jc_ad_X/pi_a/d_jc_ad_a is ",(($xsum/($#pi_xDNA+1))/$JC_divergence_xDNA)/(($asum/($#pi_aDNA+1))/$JC_divergence_aDNA),"\n";
print OUTFILE "pi_Y/d_jc_ad_Y/pi_a/d_jc_ad_a is ",(($ysum/($#pi_yDNA+1))/$JC_divergence_yDNA)/(($asum/($#pi_aDNA+1))/$JC_divergence_aDNA),"\n";

print "95% CI for aDNA S\t",$S_aDNA_boot[$lower]," - ",$S_aDNA_boot[$upper],"\n";
print "95% CI for xDNA  S\t",$S_xDNA_boot[$lower]," - ",$S_xDNA_boot[$upper],"\n";

print "95% CI for aDNA  thetaW\t",sprintf("%.5f",$S_aDNA_boot[$lower]/$a1_obs_aDNA/$aDNA_sites)," - ",sprintf("%.5f",$S_aDNA_boot[$upper]/$a1_obs_aDNA/$aDNA_sites),"\n";
print "95% CI for xDNA  thetaW\t",sprintf("%.5f",$S_xDNA_boot[$lower]/$a1_obs_xDNA/$xDNA_sites)," - ",sprintf("%.5f",$S_xDNA_boot[$upper]/$a1_obs_xDNA/$xDNA_sites),"\n";

print "95% CI for aDNA pi per site\t",sprintf("%.5f",$pi_aDNA_persite_boot[$lower])," - ",sprintf("%.5f",$pi_aDNA_persite_boot[$upper]),"\n";
print "95% CI for xDNA  pi per site\t",sprintf("%.5f",$pi_xDNA_persite_boot[$lower])," - ",sprintf("%.5f",$pi_xDNA_persite_boot[$upper]),"\n";

print "95% CI for aDNA thetapi/divergence\t",sprintf("%.5f",$thetapi_over_divergence_aDNA[$lower])," - ",sprintf("%.5f",$thetapi_over_divergence_aDNA[$upper]),"\n";
print "95% CI for xDNA  thetapi/divergence\t",sprintf("%.5f",$thetapi_over_divergence_xDNA[$lower])," - ",sprintf("%.5f",$thetapi_over_divergence_xDNA[$upper]),"\n";

print "95% CI for aDNA TajD\t",sprintf("%.3f",$aDNA_TajD_boot[$lower])," - ",sprintf("%.3f",$aDNA_TajD_boot[$upper]),"\n";
print "95% CI for xDNA  TajD\t",sprintf("%.3f",$xDNA_TajD_boot[$lower])," - ",sprintf("%.3f",$xDNA_TajD_boot[$upper]),"\n";

print "95% CI for aDNA d\t",sprintf("%.5f",$d_aDNA_boot[$lower]/$aDNA_sites)," - ",sprintf("%.5f",$d_aDNA_boot[$upper]/$aDNA_sites),"\n";
print "95% CI for xDNA  d\t",sprintf("%.5f",$d_xDNA_boot[$lower]/$xDNA_sites)," - ",sprintf("%.5f",$d_xDNA_boot[$upper]/$xDNA_sites),"\n";

print "95% CI for JC_AP aDNA d\t",sprintf("%.5f",$JC_AP_divergence_aDNA_boot[$lower])," - ",sprintf("%.5f",$JC_AP_divergence_aDNA_boot[$upper]),"\n";
print "95% CI for JC_AP xDNA  d\t",sprintf("%.5f",$JC_AP_divergence_xDNA_boot[$lower])," - ",sprintf("%.5f",$JC_AP_divergence_xDNA_boot[$upper]),"\n";


close DATAINPUT;
close OUTFILE;
print "Done with input file 1\n";

```

And I have developed this script to split up the master files using the bedfiles I made for the 2014 paper.  I had to concatenate and sort some bed files to make a pooled file for the >51000 sites like this:

```bash
cat bedfile3_51000_to_101000.bed bedfile4_101000_to_151000.bed bedfile5_151000_and_higher.bed > bedfile_51000_and_higher.bed
sed -i -e 's/ /\t/g' bedfile_51000_and_higher.bed
sort -V -k 1,1 -k2,2 bedfile_51000_and_higher.bed > bedfile_51000_and_higher_sorted.bed
```

And now I used this script to split up the vcf files into different sections depending on the distance of the sites from genes:

```perl
#!/usr/bin/perl
# This script will read in a vcf file names and 
# make and execute a GATK commandline that marks INDELs
# in a new vcf file.  

my $status;
my $infile = "recal_stampy_allsites_round2_all_confident_sites.vcf";
my $outfile1 = "recal_plusminus_1000.vcf";
my $bedfile1 = "bedfile1_genes_plusminus_1000.bed";
my $outfile2 = "recal_1000_51000.vcf";
my $bedfile2 = "bedfile2_1001_to_51000.bed";
my $outfile3 = "recal_51000plus.vcf";
my $bedfile3 = "bedfile_51000_and_higher.bed";

my $commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
$commandline = $commandline." -o ".$outfile1." --variant ".$infile;
$commandline = $commandline." -L /home/ben/2015_SulaRADtag/bed_files_perfect/".$bedfile1;
$status = system($commandline);

$commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
$commandline = $commandline." -o ".$outfile2." --variant ".$infile;
$commandline = $commandline." -L /home/ben/2015_SulaRADtag/bed_files_perfect/".$bedfile2;
$status = system($commandline);

$commandline = "java -Xmx2g -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta"; 
$commandline = $commandline." -o ".$outfile3." --variant ".$infile;
$commandline = $commandline." -L /home/ben/2015_SulaRADtag/bed_files_perfect/".$bedfile3;
$status = system($commandline);
```

As previously, now we are ready to convert the vcf files to tab delimited files like this:

```bash
~/tabix-0.2.6/bgzip XXX.vcf
~/tabix-0.2.6/tabix -p vcf XXX.vcf.gz
zcat XXX.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > XXX.vcf.gz.tab

```

And this can be followed up by adding the outgroup sequences and calculating the popgen stats as above.
