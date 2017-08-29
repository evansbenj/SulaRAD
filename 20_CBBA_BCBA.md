# Counting CBBA_BCBA sites

working in this directory:
`
/net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07
`

I made two scripts to do this. For haploid chrX files (/Performs_CBBA_BCBA_on_populations_onlychrX_haploid_outgroup.pl), use this command:


`./Performs_CBBA_BCBA_on_populations_onlychrX_haploid_outgroup.pl nonrecal_filtered_chrX_final.vcf.gz_norepeat.vcf.all_with_one_highestdepth_allele.tab 000 3_4_1_2_3 H1nigra_H2tonk_H3nem_chrX_norepeat_CBBA_BCBA_highest.jk H1nigra_H2tonk_H3nem_chrX_norepeat_CBBA_BCBA_highest.stats`

I also did this after excluding male hets and with the baboon outgroup:
```
/net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/Performs_CBBA_BCBA_on_populations_onlychrX.pl /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_5X_filtered_chrX_final.vcf.gz_norepeat_nomalehets.vcf_haploidmales_and_female.vcf.tab 000 3_4_1_2_3 /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_5X_filtered_chrX_final.vcf.gz_norepeat_nomalehets.vcf_haploidmales_and_female.vcf.tab_CBBCBCBA.jk /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_5X_filtered_chrX_final.vcf.gz_norepeat_nomalehets.vcf_haploidmales_and_female.vcf.tab_CBBCBCBA.stats
```
```
/net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/Performs_CBBA_BCBA_on_populations_onlychrX.pl /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_5X_filtered_chrX_final.vcf.gz_norepeat_nomalehets.vcf_haploidmales_and_female_with_baboon.tab 000 3_4_1_2_3 /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_5X_filtered_chrX_final.vcf.gz_norepeat_nomalehets.vcf_haploidmales_and_female_with_baboon.tab.CBBCBCBA.jk /net/infofile4-inside/volume1/scratch/ben/2016_FINAL_Sulawesi_nem_WGS/Project_MEL_11554_B01_CUS_WGS.2016-10-07/nonrecal_5X_filtered_chrX_final.vcf.gz_norepeat_nomalehets.vcf_haploidmales_and_female_with_baboon.tab.CBBCBCBA.stats
```

and this script:

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use lib qw(~/perl_modules);
use List::MoreUtils qw/ uniq /;
use List::Util 'max';

#  This program reads in a tab delimited genotype file generated
#  by the perl program "Gets_outgroup_sequence_from_axt_files.pl"
#  or from vcftools
#  and performs the population ABBA-BABA test using one outgroup
#  sequence and one or more individual sequences from three other
#  taxa (H3, H2, and H1) using only chrX. 

#  This analysis will include all positions that have data from at least 
#  one individual from species H3, H2, and H1, including those that that 
#  missing data in some individuals.

# to execute type Performs_ABBA_BABA_on_populations_onlychrX.pl inputfile.tab 1111100110000111100011100110010100000000 
# 3_6_14-18-19-20_2-3-4-5-6-7_32-33-34-35-36-37-38-39-40 output1.txt output2.txt
# where 1111100110000111100011100110010100000000 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) or is not (0) female ,and 3_6 refers to (i) the column that contains the 
# outgroup nucleotide (3 in this case for rhesus) and (ii) the column number of the first individual in the ingroup 
# (6 in this case), and 14-18-19-20, 2-3-4-5-6-7, and 32-33-34-35-36-37-38-39-40 refer to H3, H1, and H2 samples as 
# itemized below (here, they are Borneo nemestrina, hecki, and tonkeana, respectively).

# output1.txt is the one that is fed into the jacknife.R script
# output2.txt has more information for each window, including BBAA sites, fdom and dxy

# IMPORTANT: (i) and (ii) are columns beginning with 1 but (iii) is based on the individual samples such
# as enumerated below



my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $input3 = $ARGV[2];
my $outputfile = $ARGV[3];
my $outputfile2 = $ARGV[4];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the genotype input file.\n";
	exit;
}

my @temp;
my $range;
my $line_number=0;
my $counter=0;


unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2\n";
	exit;
}
print "Creating output file: $outputfile2\n";


my @sexes = split("",$ARGV[1]);
my @whotoinclude = split("_",$ARGV[2]);
my @H3=split("-",$whotoinclude[2]);
my @H1=split("-",$whotoinclude[3]);
my @H2=split("-",$whotoinclude[4]);

if($#whotoinclude != 4){
	print "Problem with number of taxa in commandline\n";
}


my @temp1;
my $previous= 0;
my $string;
my $m;
my $n;
my $w;
my $y;
my $x;
my @unique;
my $x_uniq;

my $number_of_H3_genotyped=($#H3 + 1);
my $number_of_H2_genotyped=($#H2 + 1);
my $number_of_H1_genotyped=($#H1 + 1);

my $number_of_individuals_genotyped=$#H3+$#H2+$#H1+3;

print "The number of individuals to assess is ",$number_of_individuals_genotyped,"\n";

my $number_of_female_individuals_genotyped=0;
for ($y = 0 ; $y <= $#H1 ; $y++ ) {
	if($sexes[$H1[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}
for ($y = 0 ; $y <= $#H2 ; $y++ ) {
	if($sexes[$H2[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}	
for ($y = 0 ; $y <= $#H3 ; $y++ ) {
	if($sexes[$H3[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}
print "This includes ",$number_of_female_individuals_genotyped," female(s)\n";

my $sliding_window=5000000;
my $current_window=0;
my $window_counter=0;
my $current_chromosome="blah";


my %ABBA_hash;
my $ABBA_peak_hash;
my %BABA_hash;
my $BABA_peak_hash;
my %BBAA_hash;
my $A;
my $B;
my @allelez;
my $derived;
my $ancestral;
my $H3_derived_freq;
my $H1_derived_freq;
my $H2_derived_freq;
my $H3_ancestral_freq;
my $H1_ancestral_freq;
my $H2_ancestral_freq;
my $peak_H2_H3_derived_freq;
my $peak_H1_H3_derived_freq;
my @H3allelez;
my @H2allelez;
my @H1allelez;
my $diff;
my $num_comparisons;
my %H2_H3_pairwise_divergence_per_window;
my $diffH1H3;
my $num_comparisonsH1H3;
my %H1_H3_pairwise_divergence_per_window;
my %number_of_sites_per_window;
my %number_of_sites_per_windowH1H3;
my $ABBA_peak_hashH1H3;
my $BABA_peak_hashH1H3;
my %H1_pairwise_nucleotide_diversity_per_window;
my %H2_pairwise_nucleotide_diversity_per_window;
my %H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window;
my %H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window;


my $diffH2;
my $num_comparisonsH2;
my @H2_first_allelez;
my @H2_second_allelez;
my $aa;
my $bb;

my $diffH1;
my $num_comparisonsH1;
my @H1_first_allelez;
my @H1_second_allelez;

my %f_H2_H3;
my %f_H1_H3;
my %f_dm;
my %f_dm_counter;

my %f_H2_H3_counter;
my %f_H1_H3_counter;


my @diploidoutgroup;

my $H3string;
my $H2string;
my $H1string;
my @temp43;
my @temp42;
my @temp41;
my @H3_fourderivednucs;
my @H2_fourderivednucs;
my @H1_fourderivednucs;
my $H3max;
my $counter43;
my @H3_B;
my $H1_B;
my $H1_C;
my $H2_B;
my $H2_C;
my $H3_B;
my $H3_C;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] eq '#CHROM'){
		print "The outgroup sequence is ",$temp[$whotoinclude[0]-1],"\n";
		# here we need to randomly select a base from the outgroup genotype
		# but only if the outgroup genotype is a diploid genotype
		# @diploidoutgroup=split('/',$temp[$whotoinclude[0]-1]);
		# pick a random allele
		# $temp[$whotoinclude[0]-1] = $diploidoutgroup[rand @diploidoutgroup];		
		print " The outgroup variant is ",$temp[$whotoinclude[0]-1],"\n";
		print "H3 population consists of these individuals: ";
			for ($y = 0 ; $y <= $#H3; $y++ ) {
				print $temp[$H3[$y]+$whotoinclude[1]-2]," ";
			}	
		print "\n";
		print "H1 population consists of these individuals: ";
			for ($y = 0 ; $y <= $#H1; $y++ ) {
				print $temp[$H1[$y]+$whotoinclude[1]-2]," ";
			}	
		print "\n";
		print "H2 population consists of these individuals: ";
			for ($y = 0 ; $y <= $#H2; $y++ ) {
				print $temp[$H2[$y]+$whotoinclude[1]-2]," ";
			}	
		print "\n";
	}
	else{	
		# here we need to randomly select a base from the outgroup genotype
		# but only if the outgroup genotype is a diploid genotype
		# @diploidoutgroup=split('/',$temp[$whotoinclude[0]-1]);
		# pick a random allele
		#$temp[$whotoinclude[0]-1] = $diploidoutgroup[rand @diploidoutgroup];		
		
		# This is a correction to deal with the way that the XL genome has chr9 and 10 annotated
		#if($temp[0] eq "chr9_10L"){
		#	$temp[0] = "chr9and10L";
		#}
		#elsif($temp[0] eq "chr9_10S"){
		#	$temp[0] = "chr9and10S";
		#}

		if($temp[0] ne $current_chromosome){
			$current_chromosome = $temp[0];
			$current_window = 0;
			$window_counter+=1;
			$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
		}
		until($temp[1] < ($current_window+$sliding_window)){
			$current_window = $current_window+$sliding_window;
			$window_counter+=1;
			$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
		}
		if($temp[0] eq "chrX"){
			$string=();
			if((uc $temp[$whotoinclude[0]-1] eq "A")||(uc $temp[$whotoinclude[0]-1] eq "C")||(uc $temp[$whotoinclude[0]-1] eq "T")||(uc $temp[$whotoinclude[0]-1] eq "G")){
				# the outgroup is defined
				$string=();
				$A = uc $temp[$whotoinclude[0]-1];
				$string=$string.$A;
				# now calculate the frequency of the derived allele in the H3 data
				$H3_derived_freq=0;
				$H3_ancestral_freq=0;
				$derived=0;
				$ancestral=0;
				$H3string=();
				for ($y = 0 ; $y <= $#H3; $y++ ) {
					if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'A/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'T/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'G/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'C/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne './')
						&&(length($temp[$H3[$y]+$whotoinclude[1]-2]) <= 2)){	
						#print "yo",$H3[$y]," ",$temp[$H3[$y]+$whotoinclude[1]-2],"\n";				
						@allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						$H3string=$H3string.$allelez[0];
						# check if this individual is a female before considering the second X allele
						if($sexes[$H3[$y]-1] == 1){
							if($allelez[1] eq $A){
								$ancestral+=1;
							}
							else{
								$derived+=1;
							}
							$string=$string.$allelez[1];
							$H3string=$H3string.$allelez[1];
						}	
					}	
				}
				#print "hi @H3 ",$whotoinclude[1]-2," ",$H3string,"\n";
				if(($derived+$ancestral)>0){
					$H3_derived_freq=$derived/($derived+$ancestral);
					$H3_ancestral_freq=$ancestral/($derived+$ancestral);
				}
				# now calculate the frequency of the derived allele in the H1 data
				$H1_derived_freq=0;
				$H1_ancestral_freq=0;
				$derived=0;
				$ancestral=0;
				$H1string=();
				for ($y = 0 ; $y <= $#H1; $y++ ) {
					if(($temp[$H1[$y]+$whotoinclude[1]-2] ne './.')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/N')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/.')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'A/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'C/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'G/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'T/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/.')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne './')
						&&(length($temp[$H1[$y]+$whotoinclude[1]-2]) <= 2)){
						@allelez=split('/',$temp[$H1[$y]+$whotoinclude[1]-2]);
						#print "yo",$H1[$y]," ",$temp[$H1[$y]+$whotoinclude[1]-2],"\n";				
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						$H1string=$H1string.$allelez[0];
						# check if this individual is a female before considering the second X allele
						if($sexes[$H1[$y]-1] == 1){
							if($allelez[1] eq $A){
								$ancestral+=1;
							}
							else{
								$derived+=1;
							}
							$string=$string.$allelez[1];
							$H1string=$H1string.$allelez[1];
						}
					}
				}
				#print "hi @H1 ",$whotoinclude[1]-2," ",$H1string,"\n";
				if(($derived+$ancestral)>0){
					$H1_derived_freq=$derived/($derived+$ancestral);
					$H1_ancestral_freq=$ancestral/($derived+$ancestral);
				}
				# now calculate the frequency of the derived allele in the H2 data
				$H2_derived_freq=0;
				$H2_ancestral_freq=0;
				$derived=0;
				$ancestral=0;
				$H2string=();
				for ($y = 0 ; $y <= $#H2; $y++ ) {
					if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'A/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'T/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'C/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'G/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/.')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne './')
						&&(length($temp[$H2[$y]+$whotoinclude[1]-2]) <= 2)){
						#print "yo",$H2[$y]," ",$temp[$H2[$y]+$whotoinclude[1]-2],"\n";				
						@allelez=split('/',$temp[$H2[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						$H2string=$H2string.$allelez[0];
						# check if this individual is a female before considering the second X allele
						if($sexes[$H2[$y]-1] == 1){
							if($allelez[1] eq $A){
								$ancestral+=1;
							}
							else{
								$derived+=1;
							}
							$string=$string.$allelez[1];
							$H2string=$H2string.$allelez[1];
						}	
					}	
				}
				#print "hi @H2 ",$whotoinclude[1]-2," ",$H2string,"\n";
				if(($derived+$ancestral)>0){
					$H2_derived_freq=$derived/($derived+$ancestral);
					$H2_ancestral_freq=$ancestral/($derived+$ancestral);
				}
				# only consider sites for which there are data for H1, H2, and H3
				if((defined($string))&&(($H1_derived_freq>0)||($H1_ancestral_freq>0))&&(($H2_derived_freq>0)||($H2_ancestral_freq>0))&&
					(($H3_derived_freq>0)||($H3_ancestral_freq>0))){
					# Now calculate the average pairwise divergence between H2 and H3
					$diff=0;
					$num_comparisons=0;
					@H3allelez=();
					@H2allelez=();
					for ($y = 0 ; $y <= $#H3; $y++ ) {
						if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne './')
							&&(length($temp[$H3[$y]+$whotoinclude[1]-2]) <= 2)){
							for ($w = 0 ; $w <= $#H2; $w++ ) {
								if(($temp[$H2[$w]+$whotoinclude[1]-2] ne './.')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/N')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/*')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/.')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/*')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/.')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne './')
									&&(length($temp[$H2[$w]+$whotoinclude[1]-2]) <= 2)){
									# there are data for H2 and H3
									@H3allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
									@H2allelez=split('/',$temp[$H2[$w]+$whotoinclude[1]-2]);
									# if H3 is a male, only consider first base so delete the others
									if($sexes[$H3[$y]-1] == 0){
										@H3allelez = splice @H3allelez, 0, 1;
									}	

									foreach(@H3allelez){
										if($_ ne $H2allelez[0]){
											$diff+=1; 
										}
										# check if this individual is a female before considering the second X allele
										if($sexes[$H1[$w]-1] == 1){
											if($_ ne $H2allelez[1]){
												$diff+=1;
											}
											$num_comparisons+=2;	
										}
									}
								}
							}
						}					
					}
					# Now calculate the average pairwise divergence between H1 and H3
					$diffH1H3=0;
					$num_comparisonsH1H3=0;
					@H3allelez=();
					@H1allelez=();
					for ($y = 0 ; $y <= $#H3; $y++ ) {
						if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne './')
							&&(length($temp[$H3[$y]+$whotoinclude[1]-2]) <= 2)){
							for ($w = 0 ; $w <= $#H1; $w++ ) {
								if(($temp[$H1[$w]+$whotoinclude[1]-2] ne './.')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/N')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/.')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/.')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne './')
									&&(length($temp[$H1[$w]+$whotoinclude[1]-2]) <= 2)){
									# there are data for H1 and H3
									@H3allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
									@H1allelez=split('/',$temp[$H1[$w]+$whotoinclude[1]-2]);
									# if H3 is a male, only consider first base so delete the others
									if($sexes[$H3[$y]-1] == 0){
										@H3allelez = splice @H3allelez, 0, 1;
									}	

									foreach(@H3allelez){
										if($_ ne $H1allelez[0]){
											$diffH1H3+=1; 
										}
										# check if this individual is a female before considering the second X allele
										if($sexes[$H1[$w]-1] == 1){
											if($_ ne $H1allelez[1]){
												$diffH1H3+=1;
											}
											$num_comparisonsH1H3+=2;	
										}	
									}
								}
							}
						}					
					}
					# now calculate the average pairwise divergence for H2 and H3
					if($num_comparisons>0){
						# this does not depend on the position being polymorphic
							$number_of_sites_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=1;
							$H2_H3_pairwise_divergence_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diff/$num_comparisons);
							# later we will standardize this by the number of sites per window
					}
					# now calculate the average pairwise divergence for H1 and H3
					if($num_comparisonsH1H3>0){
						# this does not depend on the position being polymorphic
							$number_of_sites_per_windowH1H3{$window_counter."_".$current_chromosome."_".$current_window}+=1;
							$H1_H3_pairwise_divergence_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diffH1H3/$num_comparisonsH1H3);
							# later we will standardize this by the number of sites per window
					}
					# Now calculate the average pairwise nucleotide diversity within H2
					# this approach is justified here: http://binhe.org/2011/12/29/calculate-nucleotide-diversity-per-base-pair-summation-method/
					# but with a formula that is probably much quicker then the stuff below pi = (2j(n-j) / n(n-1) ), where j is the number of minor
					# alleles out of n alleles total.
					$diffH2=0;
					$num_comparisonsH2=0;
					@H2_first_allelez=();
					@H2_second_allelez=();
					for ($y = 0 ; $y < $#H2; $y++ ) {
						for ($w = ($y+1) ; $w <= $#H2; $w++ ) {
							if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne './')
								&&(length($temp[$H2[$y]+$whotoinclude[1]-2]) == 3) 
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne './')
								&&(length($temp[$H2[$w]+$whotoinclude[1]-2]) <= 2)){
								# there are data for both H2 alleles
								@H2_first_allelez=split('/',$temp[$H2[$y]+$whotoinclude[1]-2]);
								@H2_second_allelez=split('/',$temp[$H2[$w]+$whotoinclude[1]-2]);

								# check if each individual is male and adjust the alleles accordingly
								if($sexes[$H2[$y]-1] == 0){
									@H2_first_allelez = splice @H2_first_allelez, 0, 1;
								}	
								if($sexes[$H2[$w]-1] == 0){
									@H2_second_allelez = splice @H2_second_allelez, 0, 1;
								}	

								# combine the alleles into one array
								@H2_first_allelez = (@H2_first_allelez, @H2_second_allelez);
								# check that the array has between 2 and 4 elements 
								# (2 if only two males are compared)
								# (3 if a male and female genotype is compared)
								# (4 if two female genotypes are compared)
								if(($#H2_first_allelez > 3)||($#H2_first_allelez < 1)){
									print "Problem 2$#H2_first_allelez with number of alleles @H2_first_allelez @H2_second_allelez\n";
								}
								else{ # calculate the average pairwise diversity for this pair of genotypes
									for ($bb = 0 ; $bb < $#H2_first_allelez; $bb++) {
										for ($aa = ($bb+1); $aa <= $#H2_first_allelez; $aa++) {
											if($H2_first_allelez[$bb] ne $H2_first_allelez[$aa]){
												$diffH2+=1; 
											}
											$num_comparisonsH2+=1;
										}
									}		
								}
							}
						}					
					}
					# for some sites, there may be only one male individual, so no pi calculation is possible 
					if($num_comparisonsH2==0){
						for ($y = 0 ; $y <= $#H2; $y++ ) {
							if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne './')
								&&(length($temp[$H2[$y]+$whotoinclude[1]-2]) <= 2)){
								if($sexes[$H2[$y]-1] == 1){
								@H2_first_allelez=split('/',$temp[$H2[$y]+$whotoinclude[1]-2]);
									if($H2_first_allelez[0] ne $H2_first_allelez[1]){
										$diffH2+=1;
									}
									$num_comparisonsH2+=1;	
								}
							}
						}	
					}	
					# now tabulate the average pairwise nucleotide diversity within H2
					if($num_comparisonsH2>0){
							# this does not depend on the position being polymorphic
								$H2_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diffH2/$num_comparisonsH2);
								$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=1;
								# later we will standardize this by the number of sites per window
					}
					# Now calculate the average pairwise nucleotide diversity within H1
					$diffH1=0;
					$num_comparisonsH1=0;
					@H1_first_allelez=();
					@H1_second_allelez=();
					for ($y = 0 ; $y < $#H1; $y++ ) {
						for ($w = ($y+1) ; $w <= $#H1; $w++ ) {
							if(($temp[$H1[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne './')
								&&(length($temp[$H1[$y]+$whotoinclude[1]-2]) <= 2) 
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne './')
								&&(length($temp[$H1[$w]+$whotoinclude[1]-2]) <= 2)){
								# there are data for both H1 alleles
								@H1_first_allelez=split('/',$temp[$H1[$y]+$whotoinclude[1]-2]);
								@H1_second_allelez=split('/',$temp[$H1[$w]+$whotoinclude[1]-2]);
								# combine the alleles into one array
								@H1_first_allelez = (@H1_first_allelez, @H1_second_allelez);
								# check that the array has between 2 and 4 elements 
								# (2 if only two males are compared)
								# (3 if a male and female genotype is compared)
								# (4 if two female genotypes are compared)
								if(($#H1_first_allelez > 3)||($#H1_first_allelez < 1)){	
									print "Problem 1$#H1_first_allelez with number of alleles @H1_first_allelez @H1_second_allelez\n";
								}
								else{
									for ($bb = 0 ; $bb < $#H1_first_allelez; $bb++) {
										for ($aa = ($bb+1); $aa <= $#H1_first_allelez; $aa++) {
											if($H1_first_allelez[$bb] ne $H1_first_allelez[$aa]){
												$diffH1+=1; 
											}
											$num_comparisonsH1+=1;
										}
									}		
								}
							}
						}					
					}
					# for some sites, there may be only one individual with a genotype
					if($num_comparisonsH1==0){
						for ($y = 0 ; $y <= $#H1; $y++ ) {
							if(($temp[$H1[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne './')
								&&(length($temp[$H1[$y]+$whotoinclude[1]-2]) <= 2)){
								if($sexes[$H1[$y]-1] == 1){
									@H1_first_allelez=split('/',$temp[$H1[$y]+$whotoinclude[1]-2]);
									if($H1_first_allelez[0] ne $H1_first_allelez[1]){
										$diffH1+=1;
									}
									$num_comparisonsH1+=1;	
								}
							}	
						}	
					}	
					# now tabulate the average pairwise nucleotide diversity within H1
					if($num_comparisonsH1>0){
							# this does not depend on the position being polymorphic
								$H1_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diffH1/$num_comparisonsH1);
								$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=1;
								# later we will standardize this by the number of sites per window
					}
					# calculate the CBBA BCBA stats
					@temp1=split('',$string);
					@temp43=split('',$H3string);
					@temp42=split('',$H2string);
					@temp41=split('',$H1string);
					$x_uniq = uniq @temp1;
					$H1_B=0;
					$H1_C=0;
					$H2_B=0;
					$H2_C=0;
					$H3_B=0;
					$H3_C=0;
					#print "hey ",$line;
					if($x_uniq > 2){
						#print "hey ",$line;
						#$BBAA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_derived_freq*$H2_derived_freq*(1-$H3_derived_freq));
						#if((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq)>0) || (($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)>0)){
						if(($H1_derived_freq*$H2_derived_freq*$H3_derived_freq)>0){
						# I have changed this from the ABBA BABA if statement because we are looking at sites that have some
						# proportion ofderived nucleotides in H1, H2, and H3 	
							# this is a polymorphic position and it is an ABBA_BABA site   
							#if($current_window ==10000000 ){
							#	print $H1_derived_freq,"\t",$H2_derived_freq,"\t",$H3_derived_freq,"\n";
							#}
						    #print $line,"\n";
						    # Define Bs as the highest frequency nonancestral allele in H3
						    # Count up Bs and Cs in H3
						    $H3_fourderivednucs[0]=0;
						    $H3_fourderivednucs[1]=0;
						    $H3_fourderivednucs[2]=0;
						    $H3_fourderivednucs[3]=0;
						    foreach(@temp43){
						    	if(uc $_ ne uc $A){
						    		if(uc $_ eq 'A'){
						    			$H3_fourderivednucs[0]+=1;
						    		}
						    		elsif(uc $_ eq 'C'){
						    			$H3_fourderivednucs[1]+=1;
						    		}
						    		elsif(uc $_ eq 'T'){
						    			$H3_fourderivednucs[2]+=1;
						    		}
						    		elsif(uc $_ eq 'G'){
						    			$H3_fourderivednucs[3]+=1;
						    		}
						    	}
						    }	
						    $H3_fourderivednucs[0]=$H3_fourderivednucs[0]/($#temp43+1);	
						    $H3_fourderivednucs[1]=$H3_fourderivednucs[1]/($#temp43+1);	
						    $H3_fourderivednucs[2]=$H3_fourderivednucs[2]/($#temp43+1);	
						    $H3_fourderivednucs[3]=$H3_fourderivednucs[3]/($#temp43+1);

						    # Now I need to figure out which nucleotide has the max frequency in H3
						    # this will be the B nucleotide.  I think the others should be the C
						    # in case of a tie, perhaps I should make all of the highest ones the B?

						    $H3max = max(@H3_fourderivednucs);

						    $counter43=0;
						    @H3_B=();
						    foreach(@H3_fourderivednucs){
						    	if($_ == $H3max){
						    		push (@H3_B,$counter43);
						    	}
						    	$counter43+=1;
						    }
						    # Now the array @H3_B has the index of all of the max nucleotides

						    # calculate the sum of all derived allele frequencies in H3
						    $H3_C=$H3_fourderivednucs[0]+$H3_fourderivednucs[1]+$H3_fourderivednucs[2]+$H3_fourderivednucs[3];
						    # now calculate the sum of the frequencies of the B derived SNP only (including ties)
						    foreach(@H3_B){
						    	$H3_B+=$H3_fourderivednucs[$_];
						    }	
						    # Now, to get the C frequency, subtract the B freq from the total.
						    $H3_C=$H3_C-$H3_B;

						    # now get the B and C frequencies from H2
						    $H2_fourderivednucs[0]=0;
						    $H2_fourderivednucs[1]=0;
						    $H2_fourderivednucs[2]=0;
						    $H2_fourderivednucs[3]=0;
						    foreach(@temp42){
						    	if(uc $_ ne uc $A){
						    		if(uc $_ eq 'A'){
						    			$H2_fourderivednucs[0]+=1;
						    		}
						    		elsif(uc $_ eq 'C'){
						    			$H2_fourderivednucs[1]+=1;
						    		}
						    		elsif(uc $_ eq 'T'){
						    			$H2_fourderivednucs[2]+=1;
						    		}
						    		elsif(uc $_ eq 'G'){
						    			$H2_fourderivednucs[3]+=1;
						    		}
						    	}
						    }
						    $H2_fourderivednucs[0]=$H2_fourderivednucs[0]/($#temp42+1);	
						    $H2_fourderivednucs[1]=$H2_fourderivednucs[1]/($#temp42+1);	
						    $H2_fourderivednucs[2]=$H2_fourderivednucs[2]/($#temp42+1);	
						    $H2_fourderivednucs[3]=$H2_fourderivednucs[3]/($#temp42+1);	
						    # first calculate the sum of all derived allele frequencies in H2
						    $H2_C=$H2_fourderivednucs[0]+$H2_fourderivednucs[1]+$H2_fourderivednucs[2]+$H2_fourderivednucs[3];
						    # now calculate the sum of the frequencies of the B derived SNP only (including ties)
						    foreach(@H3_B){
						    	$H2_B+=$H2_fourderivednucs[$_];
						    }	
						    # Now, to get the C frequency, subtract the B freq from the total.
						    $H2_C=$H2_C-$H2_B;

						    # now get the B and C frequencies from H1
						    $H1_fourderivednucs[0]=0;
						    $H1_fourderivednucs[1]=0;
						    $H1_fourderivednucs[2]=0;
						    $H1_fourderivednucs[3]=0;
						    foreach(@temp41){
						    	if(uc $_ ne uc $A){
						    		if(uc $_ eq 'A'){
						    			$H1_fourderivednucs[0]+=1;
						    		}
						    		elsif(uc $_ eq 'C'){
						    			$H1_fourderivednucs[1]+=1;
						    		}
						    		elsif(uc $_ eq 'T'){
						    			$H1_fourderivednucs[2]+=1;
						    		}
						    		elsif(uc $_ eq 'G'){
						    			$H1_fourderivednucs[3]+=1;
						    		}
						    	}
						    }	
						    $H1_fourderivednucs[0]=$H1_fourderivednucs[0]/($#temp41+1);	
						    $H1_fourderivednucs[1]=$H1_fourderivednucs[1]/($#temp41+1);	
						    $H1_fourderivednucs[2]=$H1_fourderivednucs[2]/($#temp41+1);	
						    $H1_fourderivednucs[3]=$H1_fourderivednucs[3]/($#temp41+1);	
						    # first calculate the sum of all derived allele frequencies in H1
						    $H1_C=$H1_fourderivednucs[0]+$H1_fourderivednucs[1]+$H1_fourderivednucs[2]+$H1_fourderivednucs[3];
						    # now calculate the sum of the frequencies of the B derived SNP only (including ties)
						    foreach(@H3_B){
						    	$H1_B+=$H1_fourderivednucs[$_];
						    }	
						    # Now, to get the C frequenci, subtract the B freq from the total.
						    $H1_C=$H1_C-$H1_B;

						    # I am storing the CBBA proportions per site in the ABBA hash
						    # and the BCBA proportions per site in the BABA hash

						    # I will count Bs as the frequency in H1 or H2 of the most frequent derived nucleotide in H3
						    # including all of those that are tied in H3 for the most frequent
						    # C is the sum of the frequency in H1 or H2 of all other derived alleles 
						    # for H3, a B is just the sum of the frequency of all derived nucleotides (this is
						   	# the sum of Bs and Cs in H3 and is equal to the $H3_derived_freq calculated above)

						    # Cs are counted in H1 and H2 as the sum of each frequency of each derived nucleotide in H1 or H2 times 
						    # (1- the corresponding frequency in H3)

						    #if($H1_C ne $H2_C){
						    	print $line,"\n";
						    	print $H1string," ",$H2string," ",$H3string," ",$A," ",$H1_B," ",$H1_C," ",$H2_B," ",$H2_C," ",$H3_B," ",$H3_C,"\n";
							    #print $H1_C," ",$H2_C,"\n";
							#}

							$ABBA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_C*$H2_B*$H3_derived_freq);
							#$ABBA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_C);
							
							# $ABBA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq);
							
							$BABA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_B*$H2_C*$H3_derived_freq);
							#$BABA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H2_C);

							# $BABA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq);

							$BBAA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_C*$H2_C*$H3_derived_freq);

							if($H2_derived_freq > $H3_derived_freq){
								$peak_H2_H3_derived_freq = $H2_derived_freq;
							}
							else{
								$peak_H2_H3_derived_freq = $H3_derived_freq;	
							}
							$ABBA_peak_hash=((1-$H1_derived_freq)*$peak_H2_H3_derived_freq*$peak_H2_H3_derived_freq);
							$BABA_peak_hash=($H1_derived_freq*(1-$peak_H2_H3_derived_freq)*$peak_H2_H3_derived_freq);
							# here we are calculating stats for f for H2 and H3.
							# we need to do this for each site because some of them can be undefined and should
							# not be included in the average
							if(($ABBA_peak_hash - $BABA_peak_hash) != 0){
								$f_H2_H3{$window_counter."_".$current_chromosome."_".$current_window} += ((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq))-(($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)))/ ($ABBA_peak_hash - $BABA_peak_hash);
								$f_H2_H3_counter{$window_counter."_".$current_chromosome."_".$current_window} += 1;
							}	
							# here we are calculating stats for f and the assignment of H1 and H2 needs to be switched.
							if($H1_derived_freq > $H3_derived_freq){
								$peak_H1_H3_derived_freq = $H1_derived_freq;
							}
							else{
								$peak_H1_H3_derived_freq = $H3_derived_freq;	
							}
							$BABA_peak_hashH1H3=($peak_H1_H3_derived_freq*(1-$H2_derived_freq)*$peak_H1_H3_derived_freq);
							$ABBA_peak_hashH1H3=((1-$peak_H1_H3_derived_freq)*$H2_derived_freq*$peak_H1_H3_derived_freq);
							if(($BABA_peak_hashH1H3 - $ABBA_peak_hashH1H3) != 0){
								$f_H1_H3{$window_counter."_".$current_chromosome."_".$current_window} += ((($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq))-(((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq)))/ ($BABA_peak_hashH1H3 - $ABBA_peak_hashH1H3);
								$f_H1_H3_counter{$window_counter."_".$current_chromosome."_".$current_window} += 1;
							}
							if(($H2_derived_freq >= $H1_derived_freq)&&(($ABBA_peak_hash - $BABA_peak_hash)!=0)){
								$f_dm{$window_counter."_".$current_chromosome."_".$current_window} += ((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq))-(($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)))/ ($ABBA_peak_hash - $BABA_peak_hash);
								$f_dm_counter{$window_counter."_".$current_chromosome."_".$current_window} += 1;
							}
							elsif(($ABBA_peak_hashH1H3 - $BABA_peak_hashH1H3)!=0){
								$f_dm{$window_counter."_".$current_chromosome."_".$current_window} += ((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq))-(($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)))/ -($ABBA_peak_hashH1H3 - $BABA_peak_hashH1H3);	
								$f_dm_counter{$window_counter."_".$current_chromosome."_".$current_window} += 1;
							}
	
						}
					}
				}	
			}
		}
	}
}


close DATAINPUT;

# now merge the hash keys
my @common_keys = ();

foreach (keys %ABBA_hash) {
	push(@common_keys, $_);
}

foreach (keys %BABA_hash) {
	push(@common_keys, $_) unless exists $ABBA_hash{$_};
}

my @common_keys2 = ();

foreach (keys %ABBA_hash) {
	push(@common_keys2, $_);
}

foreach (keys %BBAA_hash) {
	push(@common_keys2, $_) unless exists $ABBA_hash{$_};
}


my @out = keys %{{map {($_ => 1)} (@common_keys, @common_keys2)}};

@out = map  { $_->[0] }
             sort { $a->[1] <=> $b->[1] }
             map  { [$_, $_=~/(\d+)/] }
                 @out;



foreach (@out) {
	@temp1=split('_',$_);
	print OUTFILE $temp1[1],"\t",$temp1[2]+1,"\t",$temp1[2]+$sliding_window,"\t";
	if(defined($ABBA_hash{$_})){
		print OUTFILE $ABBA_hash{$_},"\t";
	}
	else{
		print OUTFILE "0\t";
	}
	if(defined($BABA_hash{$_})){
		print OUTFILE $BABA_hash{$_},"\t0\t0\t0\t0\n";
	}
	else{
		print OUTFILE "0\t0\t0\t0\t0\n";
	}
}
if($#out == -1){
	print OUTFILE "chr0\t1\t500000\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
}

close OUTFILE;


print OUTFILE2 "chromosome\tbegin\tend\tCBBA\tBCBA\tCCBA\tD\tfdH2H3\tfH1H3\tf_dm\tdH2H3\tnum_sites_per_windowH2H3\tdH1H3\tnum_sites_per_windowH1H3\tH2pi\tnumsitesH2pi\tH1pi\tnumsitesH1pi\n";
foreach (@out) {
	@temp1=split('_',$_);
	print OUTFILE2 $temp1[1],"\t",$temp1[2]+1,"\t",$temp1[2]+$sliding_window,"\t";
	if(defined($ABBA_hash{$_})){
		print OUTFILE2 $ABBA_hash{$_},"\t";
	}
	else{
		print OUTFILE2 "0\t";
	}
	if(defined($BABA_hash{$_})){
		print OUTFILE2 $BABA_hash{$_},"\t";
	}
	else{
		print OUTFILE2 "0\t";
	}
	if(defined($BBAA_hash{$_})){
		print OUTFILE2 $BBAA_hash{$_},"\t";
	}
	else{
		print OUTFILE2 "0\t";
	}
	#print D for this window
	if((defined($ABBA_hash{$_}))&&(defined($BBAA_hash{$_}))){
		if(($ABBA_hash{$_}+$BABA_hash{$_})>0){
			print OUTFILE2 ($ABBA_hash{$_}-$BABA_hash{$_}) / ($ABBA_hash{$_}+$BABA_hash{$_}),"\t";
		}
		else{
			print OUTFILE2 "NAN\t";
		}
	}
	else{
		print OUTFILE2 "NAN\t";
	}

	#print fd H2H3 for this window
	if(defined($f_H2_H3_counter{$_})){
		print OUTFILE2 $f_H2_H3{$_}/$f_H2_H3_counter{$_},"\t";
	}
	else{
		print OUTFILE2 "NAN\t";
	}
	#print fd H1H3 for this window
	if(defined($f_H1_H3_counter{$_})){
		print OUTFILE2 $f_H1_H3{$_}/$f_H1_H3_counter{$_},"\t";
	}
	else{
		print OUTFILE2 "NAN\t";
	}	
	#print f_dm for this window
	if(defined($f_dm{$_})){
		print OUTFILE2 $f_dm{$_}/$f_dm_counter{$_},"\t";
	}
	else{
		print OUTFILE2 "NAN\t";
	}	
	#print average H2H3 pairwise divergence for this window and number of H2H3 sites
	if(defined($number_of_sites_per_window{$_})){
		if($number_of_sites_per_window{$_}>0){
			print OUTFILE2 ($H2_H3_pairwise_divergence_per_window{$_}/$number_of_sites_per_window{$_}),"\t",$number_of_sites_per_window{$_},"\t";
		}
		else{
			print OUTFILE2 "NAN\t",$number_of_sites_per_window{$_},"\t";
		}
	}
	else{
		print OUTFILE2 "NAN\tNAN\t";
	}
	#print average H1H3 pairwise divergence for this window and number of H1H3 sites
	if(defined($number_of_sites_per_windowH1H3{$_})){
		if($number_of_sites_per_windowH1H3{$_} > 0){
			print OUTFILE2 ($H1_H3_pairwise_divergence_per_window{$_}/$number_of_sites_per_windowH1H3{$_}),"\t",$number_of_sites_per_windowH1H3{$_},"\t";
		}
		else{
			print OUTFILE2 "NAN\t",$number_of_sites_per_windowH1H3{$_},"\t";
		}
	}
	else{
		print OUTFILE2 "NAN\tNAN\t";
	}
	#print H2 pairwise nucleotide diversity per site for this window and number of H2 pairwise nucleotide diversity sites
	if(defined($H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_})){
		if($H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}>0){
			print OUTFILE2 ($H2_pairwise_nucleotide_diversity_per_window{$_}/$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}),"\t",$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\t";
		}
		else{
			print OUTFILE2 "NAN\t",$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\t";
		}
	}
	else{
		print OUTFILE2 "NAN\tNAN\t";
	}
	#print H1 pairwise nucleotide diversity per site for this window and number of H1 pairwise nucleotide diversity sites
	if(defined($H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_})){
		if($H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}>0){
			print OUTFILE2 ($H1_pairwise_nucleotide_diversity_per_window{$_}/$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}),"\t",$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\n";
		}
		else{
			print OUTFILE2 "NAN\t",$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\n";
		}
	}
	else{
		print OUTFILE2 "NAN\tNAN\n";
	}
}
if($#out == -1){
	print OUTFILE2 "chr0\t1\t500000\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
}

close OUTFILE2;



```

For diploid autosomal data (Performs_CBBA_BCBA_on_populations_autosomes_haploid_outgroup.pl):

```
#!/usr/bin/env perl
use strict;
use warnings;
use lib qw(~/perl_modules);
use List::MoreUtils qw/ uniq /;
use List::Util 'max';

#  This program reads in a tab delimited genotype file generated
#  by the perl program "Gets_outgroup_sequence_from_axt_files.pl"
#  or from vcftools
#  and performs the population ABBA-BABA test using one outgroup
#  sequence and one or more individual sequences from three other
#  taxa (H3, H2, and H1) using only chrX. 

#  This analysis will include all positions that have data from at least 
#  one individual from species H3, H2, and H1, including those that that 
#  missing data in some individuals.

# to execute type Performs_ABBA_BABA_on_populations_onlychrX.pl inputfile.tab 1111100110000111100011100110010100000000 
# 3_6_14-18-19-20_2-3-4-5-6-7_32-33-34-35-36-37-38-39-40 output1.txt output2.txt
# where 1111100110000111100011100110010100000000 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) or is not (0) female ,and 3_6 refers to (i) the column that contains the 
# outgroup nucleotide (3 in this case for rhesus) and (ii) the column number of the first individual in the ingroup 
# (6 in this case), and 14-18-19-20, 2-3-4-5-6-7, and 32-33-34-35-36-37-38-39-40 refer to H3, H1, and H2 samples as 
# itemized below (here, they are Borneo nemestrina, hecki, and tonkeana, respectively).

# output1.txt is the one that is fed into the jacknife.R script
# output2.txt has more information for each window, including BBAA sites, fdom and dxy

# IMPORTANT: (i) and (ii) are columns beginning with 1 but (iii) is based on the individual samples such
# as enumerated below




my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $input3 = $ARGV[2];
my $outputfile = $ARGV[3];
my $outputfile2 = $ARGV[4];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the genotype input file.\n";
	exit;
}

my @temp;
my $range;
my $line_number=0;
my $counter=0;


unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2\n";
	exit;
}
print "Creating output file: $outputfile2\n";


my @sexes = split("",$ARGV[1]);
my @whotoinclude = split("_",$ARGV[2]);
my @H3=split("-",$whotoinclude[2]);
my @H1=split("-",$whotoinclude[3]);
my @H2=split("-",$whotoinclude[4]);

if($#whotoinclude != 4){
	print "Problem with number of taxa in commandline\n";
}


my @temp1;
my $previous= 0;
my $string;
my $m;
my $n;
my $w;
my $y;
my $x;
my @unique;
my $x_uniq;

my $number_of_H3_genotyped=($#H3 + 1);
my $number_of_H2_genotyped=($#H2 + 1);
my $number_of_H1_genotyped=($#H1 + 1);

my $number_of_individuals_genotyped=$#H3+$#H2+$#H1+3;

print "The number of individuals to assess is ",$number_of_individuals_genotyped,"\n";

my $number_of_female_individuals_genotyped=0;
for ($y = 0 ; $y <= $#H1 ; $y++ ) {
	if($sexes[$H1[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}
for ($y = 0 ; $y <= $#H2 ; $y++ ) {
	if($sexes[$H2[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}	
for ($y = 0 ; $y <= $#H3 ; $y++ ) {
	if($sexes[$H3[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}
print "This includes ",$number_of_female_individuals_genotyped," female(s)\n";

my $sliding_window=5000000;
my $current_window=0;
my $window_counter=0;
my $current_chromosome="blah";


my %ABBA_hash;
my $ABBA_peak_hash;
my %BABA_hash;
my $BABA_peak_hash;
my %BBAA_hash;
my $A;
my $B;
my @allelez;
my $derived;
my $ancestral;
my $H3_derived_freq;
my $H1_derived_freq;
my $H2_derived_freq;
my $H3_ancestral_freq;
my $H1_ancestral_freq;
my $H2_ancestral_freq;
my $peak_H2_H3_derived_freq;
my $peak_H1_H3_derived_freq;
my @H3allelez;
my @H2allelez;
my @H1allelez;
my $diff;
my $num_comparisons;
my %H2_H3_pairwise_divergence_per_window;
my $diffH1H3;
my $num_comparisonsH1H3;
my %H1_H3_pairwise_divergence_per_window;
my %number_of_sites_per_window;
my %number_of_sites_per_windowH1H3;
my $ABBA_peak_hashH1H3;
my $BABA_peak_hashH1H3;
my %H1_pairwise_nucleotide_diversity_per_window;
my %H2_pairwise_nucleotide_diversity_per_window;
my %H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window;
my %H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window;


my $diffH2;
my $num_comparisonsH2;
my @H2_first_allelez;
my @H2_second_allelez;
my $aa;
my $bb;

my $diffH1;
my $num_comparisonsH1;
my @H1_first_allelez;
my @H1_second_allelez;

my %f_H2_H3;
my %f_H1_H3;
my %f_dm;
my %f_dm_counter;

my %f_H2_H3_counter;
my %f_H1_H3_counter;


my @diploidoutgroup;

my $H3string;
my $H2string;
my $H1string;
my @temp43;
my @temp42;
my @temp41;
my @H3_fourderivednucs;
my @H2_fourderivednucs;
my @H1_fourderivednucs;
my $H3max;
my $counter43;
my @H3_B;
my $H1_B;
my $H1_C;
my $H2_B;
my $H2_C;
my $H3_B;
my $H3_C;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] eq '#CHROM'){
		print "The outgroup sequence is ",$temp[$whotoinclude[0]-1],"\n";
		# here we need to randomly select a base from the outgroup genotype
		# but only if the outgroup genotype is a diploid genotype
		# @diploidoutgroup=split('/',$temp[$whotoinclude[0]-1]);
		# pick a random allele
		# $temp[$whotoinclude[0]-1] = $diploidoutgroup[rand @diploidoutgroup];		
		print " The outgroup variant is ",$temp[$whotoinclude[0]-1],"\n";
		print "H3 population consists of these individuals: ";
			for ($y = 0 ; $y <= $#H3; $y++ ) {
				print $temp[$H3[$y]+$whotoinclude[1]-2]," ";
			}	
		print "\n";
		print "H1 population consists of these individuals: ";
			for ($y = 0 ; $y <= $#H1; $y++ ) {
				print $temp[$H1[$y]+$whotoinclude[1]-2]," ";
			}	
		print "\n";
		print "H2 population consists of these individuals: ";
			for ($y = 0 ; $y <= $#H2; $y++ ) {
				print $temp[$H2[$y]+$whotoinclude[1]-2]," ";
			}	
		print "\n";
	}
	else{	
		# here we need to randomly select a base from the outgroup genotype
		# but only if the outgroup genotype is a diploid genotype
		#@diploidoutgroup=split('/',$temp[$whotoinclude[0]-1]);
		# pick a random allele
		#$temp[$whotoinclude[0]-1] = $diploidoutgroup[rand @diploidoutgroup];		
		
		# This is a correction to deal with the way that the XL genome has chr9 and 10 annotated
		#if($temp[0] eq "chr9_10L"){
		#	$temp[0] = "chr9and10L";
		#}
		#elsif($temp[0] eq "chr9_10S"){
		#	$temp[0] = "chr9and10S";
		#}

		if($temp[0] ne $current_chromosome){
			$current_chromosome = $temp[0];
			$current_window = 0;
			$window_counter+=1;
			$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
		}
		until($temp[1] < ($current_window+$sliding_window)){
			$current_window = $current_window+$sliding_window;
			$window_counter+=1;
			$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
			$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}=0;
		}
		if(($temp[0] ne "chrX")&&($temp[0] ne "chrY")&&($temp[0] ne "chrM")){
			$string=();
			if((uc $temp[$whotoinclude[0]-1] eq "A")||(uc $temp[$whotoinclude[0]-1] eq "C")||(uc $temp[$whotoinclude[0]-1] eq "T")||(uc $temp[$whotoinclude[0]-1] eq "G")){
				# the outgroup is defined
				$string=();
				$A = uc $temp[$whotoinclude[0]-1];
				$string=$string.$A;
				# now calculate the frequency of the derived allele in the H3 data
				$H3_derived_freq=0;
				$H3_ancestral_freq=0;
				$derived=0;
				$ancestral=0;
				$H3string=();
				for ($y = 0 ; $y <= $#H3; $y++ ) {
					if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'A/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'T/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'G/*')
						&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'C/*')
						&&(length($temp[$H3[$y]+$whotoinclude[1]-2]) == 3)){	
						#print "yo",$H3[$y]," ",$temp[$H3[$y]+$whotoinclude[1]-2],"\n";				
						@allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						$H3string=$H3string.$allelez[0];
						if($allelez[1] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[1];
						$H3string=$H3string.$allelez[1];
					}	
				}
				#print "hi @H3 ",$whotoinclude[1]-2," ",$H3string,"\n";
				if(($derived+$ancestral)>0){
					$H3_derived_freq=$derived/($derived+$ancestral);
					$H3_ancestral_freq=$ancestral/($derived+$ancestral);
				}
				# now calculate the frequency of the derived allele in the H1 data
				$H1_derived_freq=0;
				$H1_ancestral_freq=0;
				$derived=0;
				$ancestral=0;
				$H1string=();
				for ($y = 0 ; $y <= $#H1; $y++ ) {
					if(($temp[$H1[$y]+$whotoinclude[1]-2] ne './.')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/N')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/.')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'A/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'C/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'G/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'T/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/*')
						&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/.')
						&&(length($temp[$H1[$y]+$whotoinclude[1]-2]) == 3)){
						@allelez=split('/',$temp[$H1[$y]+$whotoinclude[1]-2]);
						#print "yo",$H1[$y]," ",$temp[$H1[$y]+$whotoinclude[1]-2],"\n";				
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						$H1string=$H1string.$allelez[0];
						if($allelez[1] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[1];
						$H1string=$H1string.$allelez[1];
					}
				}
				#print "hi @H1 ",$whotoinclude[1]-2," ",$H1string,"\n";
				if(($derived+$ancestral)>0){
					$H1_derived_freq=$derived/($derived+$ancestral);
					$H1_ancestral_freq=$ancestral/($derived+$ancestral);
				}
				# now calculate the frequency of the derived allele in the H2 data
				$H2_derived_freq=0;
				$H2_ancestral_freq=0;
				$derived=0;
				$ancestral=0;
				$H2string=();
				for ($y = 0 ; $y <= $#H2; $y++ ) {
					if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'A/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'T/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'C/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'G/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')
						&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/.')
						&&(length($temp[$H2[$y]+$whotoinclude[1]-2]) == 3)){
						#print "yo",$H2[$y]," ",$temp[$H2[$y]+$whotoinclude[1]-2],"\n";				
						@allelez=split('/',$temp[$H2[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						$H2string=$H2string.$allelez[0];
						if($allelez[1] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[1];
						$H2string=$H2string.$allelez[1];
					}	
				}
				#print "hi @H2 ",$whotoinclude[1]-2," ",$H2string,"\n";
				if(($derived+$ancestral)>0){
					$H2_derived_freq=$derived/($derived+$ancestral);
					$H2_ancestral_freq=$ancestral/($derived+$ancestral);
				}
				# only consider sites for which there are data for H1, H2, and H3
				if((defined($string))&&(($H1_derived_freq>0)||($H1_ancestral_freq>0))&&(($H2_derived_freq>0)||($H2_ancestral_freq>0))&&
					(($H3_derived_freq>0)||($H3_ancestral_freq>0))){
					# Now calculate the average pairwise divergence between H2 and H3
					$diff=0;
					$num_comparisons=0;
					@H3allelez=();
					@H2allelez=();
					for ($y = 0 ; $y <= $#H3; $y++ ) {
						if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne './')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/')
							&&(length($temp[$H3[$y]+$whotoinclude[1]-2]) == 3)){
							for ($w = 0 ; $w <= $#H2; $w++ ) {
								if(($temp[$H2[$w]+$whotoinclude[1]-2] ne './.')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne './')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/N')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/*')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/.')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/*')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/.')
									&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/')
									&&(length($temp[$H2[$w]+$whotoinclude[1]-2]) == 3)){
									# there are data for H2 and H3
									@H3allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
									@H2allelez=split('/',$temp[$H2[$w]+$whotoinclude[1]-2]);
									foreach(@H3allelez){
										if($_ ne $H2allelez[0]){
											$diff+=1; 
										}
										if($_ ne $H2allelez[1]){
											$diff+=1;
										}
										$num_comparisons+=2;	
									}
								}
							}
						}					
					}
					# Now calculate the average pairwise divergence between H1 and H3
					$diffH1H3=0;
					$num_comparisonsH1H3=0;
					@H3allelez=();
					@H1allelez=();
					for ($y = 0 ; $y <= $#H3; $y++ ) {
						if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')
							&&($temp[$H3[$y]+$whotoinclude[1]-2] ne './')
							&&(length($temp[$H3[$y]+$whotoinclude[1]-2]) == 3)){
							for ($w = 0 ; $w <= $#H1; $w++ ) {
								if(($temp[$H1[$w]+$whotoinclude[1]-2] ne './.')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/N')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/.')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/*')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/.')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/')
									&&($temp[$H1[$w]+$whotoinclude[1]-2] ne './')
									&&(length($temp[$H1[$w]+$whotoinclude[1]-2]) == 3)){
									# there are data for H1 and H3
									@H3allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
									@H1allelez=split('/',$temp[$H1[$w]+$whotoinclude[1]-2]);
									foreach(@H3allelez){
										if($_ ne $H1allelez[0]){
											$diffH1H3+=1; 
										}
										if($_ ne $H1allelez[1]){
											$diffH1H3+=1;
										}
										$num_comparisonsH1H3+=2;	
									}
								}
							}
						}					
					}
					# now calculate the average pairwise divergence for H2 and H3
					if($num_comparisons>0){
						# this does not depend on the position being polymorphic
							$number_of_sites_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=1;
							$H2_H3_pairwise_divergence_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diff/$num_comparisons);
							# later we will standardize this by the number of sites per window
					}
					# now calculate the average pairwise divergence for H1 and H3
					if($num_comparisonsH1H3>0){
						# this does not depend on the position being polymorphic
							$number_of_sites_per_windowH1H3{$window_counter."_".$current_chromosome."_".$current_window}+=1;
							$H1_H3_pairwise_divergence_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diffH1H3/$num_comparisonsH1H3);
							# later we will standardize this by the number of sites per window
					}
					# Now calculate the average pairwise nucleotide diversity within H2
					# this approach is justified here: http://binhe.org/2011/12/29/calculate-nucleotide-diversity-per-base-pair-summation-method/
					# but with a formula that is probably much quicker then the stuff below pi = (2j(n-j) / n(n-1) ), where j is the number of minor
					# alleles out of n alleles total.
					$diffH2=0;
					$num_comparisonsH2=0;
					@H2_first_allelez=();
					@H2_second_allelez=();
					for ($y = 0 ; $y < $#H2; $y++ ) {
						for ($w = ($y+1) ; $w <= $#H2; $w++ ) {
							if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne './')
								&&(length($temp[$H2[$y]+$whotoinclude[1]-2]) == 3) 
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/')
								&&($temp[$H2[$w]+$whotoinclude[1]-2] ne './')
								&&(length($temp[$H2[$w]+$whotoinclude[1]-2]) == 3)){
								# there are data for both H2 alleles
								@H2_first_allelez=split('/',$temp[$H2[$y]+$whotoinclude[1]-2]);
								@H2_second_allelez=split('/',$temp[$H2[$w]+$whotoinclude[1]-2]);
								# combine the alleles into one array
								@H2_first_allelez = (@H2_first_allelez, @H2_second_allelez);
								# check that the array has 4 elements 
								if($#H2_first_allelez != 3){
									print "Problem 2$#H2_first_allelez with number of alleles @H2_first_allelez @H2_second_allelez\n";
								}
								else{ # calculate the average pairwise diversity for this pair of genotypes
									for ($bb = 0 ; $bb < $#H2_first_allelez; $bb++) {
										for ($aa = ($bb+1); $aa <= $#H2_first_allelez; $aa++) {
											if($H2_first_allelez[$bb] ne $H2_first_allelez[$aa]){
												$diffH2+=1; 
											}
											$num_comparisonsH2+=1;
										}
									}		
								}
							}
						}					
					}
					# for some sites, there may be only one male individual, so no pi calculation is possible 
					if($num_comparisonsH2==0){
						for ($y = 0 ; $y <= $#H2; $y++ ) {
							if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/')
								&&($temp[$H2[$y]+$whotoinclude[1]-2] ne './')
								&&(length($temp[$H2[$y]+$whotoinclude[1]-2]) == 3)){
								if($sexes[$H2[$y]-1] == 1){
								@H2_first_allelez=split('/',$temp[$H2[$y]+$whotoinclude[1]-2]);
									if($H2_first_allelez[0] ne $H2_first_allelez[1]){
										$diffH2+=1;
									}
									$num_comparisonsH2+=1;	
								}
							}
						}	
					}	
					# now tabulate the average pairwise nucleotide diversity within H2
					if($num_comparisonsH2>0){
							# this does not depend on the position being polymorphic
								$H2_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diffH2/$num_comparisonsH2);
								$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=1;
								# later we will standardize this by the number of sites per window
					}
					# Now calculate the average pairwise nucleotide diversity within H1
					$diffH1=0;
					$num_comparisonsH1=0;
					@H1_first_allelez=();
					@H1_second_allelez=();
					for ($y = 0 ; $y < $#H1; $y++ ) {
						for ($w = ($y+1) ; $w <= $#H1; $w++ ) {
							if(($temp[$H1[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne './')
								&&(length($temp[$H1[$y]+$whotoinclude[1]-2]) == 3) 
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/.')
								&&($temp[$H1[$w]+$whotoinclude[1]-2] ne './')
								&&(length($temp[$H1[$w]+$whotoinclude[1]-2]) == 3)){
								# there are data for both H1 alleles
								@H1_first_allelez=split('/',$temp[$H1[$y]+$whotoinclude[1]-2]);
								@H1_second_allelez=split('/',$temp[$H1[$w]+$whotoinclude[1]-2]);
								# combine the alleles into one array
								@H1_first_allelez = (@H1_first_allelez, @H1_second_allelez);
								# check that the array has 4 elements 
								if($#H1_first_allelez != 3){	
									print "Problem 1$#H1_first_allelez with number of alleles @H1_first_allelez @H1_second_allelez\n";
								}
								else{
									for ($bb = 0 ; $bb < $#H1_first_allelez; $bb++) {
										for ($aa = ($bb+1); $aa <= $#H1_first_allelez; $aa++) {
											if($H1_first_allelez[$bb] ne $H1_first_allelez[$aa]){
												$diffH1+=1; 
											}
											$num_comparisonsH1+=1;
										}
									}		
								}
							}
						}					
					}
					# for some sites, there may be only one individual with a genotype
					if($num_comparisonsH1==0){
						for ($y = 0 ; $y <= $#H1; $y++ ) {
							if(($temp[$H1[$y]+$whotoinclude[1]-2] ne './.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/N')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/.')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/*')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/')
								&&($temp[$H1[$y]+$whotoinclude[1]-2] ne './')
								&&(length($temp[$H1[$y]+$whotoinclude[1]-2]) == 3)){
								@H1_first_allelez=split('/',$temp[$H1[$y]+$whotoinclude[1]-2]);
								if($H1_first_allelez[0] ne $H1_first_allelez[1]){
									$diffH1+=1;
								}
								$num_comparisonsH1+=1;	
							}
						}	
					}	
					# now tabulate the average pairwise nucleotide diversity within H1
					if($num_comparisonsH1>0){
							# this does not depend on the position being polymorphic
								$H1_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diffH1/$num_comparisonsH1);
								$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=1;
								# later we will standardize this by the number of sites per window
					}
					# calculate the CBBA BCBA stats
					@temp1=split('',$string);
					@temp43=split('',$H3string);
					@temp42=split('',$H2string);
					@temp41=split('',$H1string);
					$x_uniq = uniq @temp1;
					$H1_B=0;
					$H1_C=0;
					$H2_B=0;
					$H2_C=0;
					$H3_B=0;
					$H3_C=0;

					if($x_uniq > 2){
						#$BBAA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_derived_freq*$H2_derived_freq*(1-$H3_derived_freq));
						#if((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq)>0) || (($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)>0)){
						if(($H1_derived_freq*$H2_derived_freq*$H3_derived_freq)>0){
						# I have changed this from the ABBA BABA if statement because we are looking at sites that have some
						# proportion ofderived nucleotides in H1, H2, and H3 	
							# this is a polymorphic position and it is an ABBA_BABA site   
							#if($current_window ==10000000 ){
							#	print $H1_derived_freq,"\t",$H2_derived_freq,"\t",$H3_derived_freq,"\n";
							#}
						    #print $line,"\n";
						    # Define Bs as the highest frequency nonancestral allele in H3
						    # Count up Bs and Cs in H3
						    $H3_fourderivednucs[0]=0;
						    $H3_fourderivednucs[1]=0;
						    $H3_fourderivednucs[2]=0;
						    $H3_fourderivednucs[3]=0;
						    foreach(@temp43){
						    	if(uc $_ ne uc $A){
						    		if(uc $_ eq 'A'){
						    			$H3_fourderivednucs[0]+=1;
						    		}
						    		elsif(uc $_ eq 'C'){
						    			$H3_fourderivednucs[1]+=1;
						    		}
						    		elsif(uc $_ eq 'T'){
						    			$H3_fourderivednucs[2]+=1;
						    		}
						    		elsif(uc $_ eq 'G'){
						    			$H3_fourderivednucs[3]+=1;
						    		}
						    	}
						    }	
						    $H3_fourderivednucs[0]=$H3_fourderivednucs[0]/($#temp43+1);	
						    $H3_fourderivednucs[1]=$H3_fourderivednucs[1]/($#temp43+1);	
						    $H3_fourderivednucs[2]=$H3_fourderivednucs[2]/($#temp43+1);	
						    $H3_fourderivednucs[3]=$H3_fourderivednucs[3]/($#temp43+1);

						    # Now I need to figure out which nucleotide has the max frequency in H3
						    # this will be the B nucleotide.  I think the others should be the C
						    # in case of a tie, perhaps I should make all of the highest ones the B?

						    $H3max = max(@H3_fourderivednucs);

						    $counter43=0;
						    @H3_B=();
						    foreach(@H3_fourderivednucs){
						    	if($_ == $H3max){
						    		push (@H3_B,$counter43);
						    	}
						    	$counter43+=1;
						    }
						    # Now the array @H3_B has the index of all of the max nucleotides

						    # calculate the sum of all derived allele frequencies in H3
						    $H3_C=$H3_fourderivednucs[0]+$H3_fourderivednucs[1]+$H3_fourderivednucs[2]+$H3_fourderivednucs[3];
						    # now calculate the sum of the frequencies of the B derived SNP only (including ties)
						    foreach(@H3_B){
						    	$H3_B+=$H3_fourderivednucs[$_];
						    }	
						    # Now, to get the C frequency, subtract the B freq from the total.
						    $H3_C=$H3_C-$H3_B;

						    # now get the B and C frequencies from H2
						    $H2_fourderivednucs[0]=0;
						    $H2_fourderivednucs[1]=0;
						    $H2_fourderivednucs[2]=0;
						    $H2_fourderivednucs[3]=0;
						    foreach(@temp42){
						    	if(uc $_ ne uc $A){
						    		if(uc $_ eq 'A'){
						    			$H2_fourderivednucs[0]+=1;
						    		}
						    		elsif(uc $_ eq 'C'){
						    			$H2_fourderivednucs[1]+=1;
						    		}
						    		elsif(uc $_ eq 'T'){
						    			$H2_fourderivednucs[2]+=1;
						    		}
						    		elsif(uc $_ eq 'G'){
						    			$H2_fourderivednucs[3]+=1;
						    		}
						    	}
						    }
						    $H2_fourderivednucs[0]=$H2_fourderivednucs[0]/($#temp42+1);	
						    $H2_fourderivednucs[1]=$H2_fourderivednucs[1]/($#temp42+1);	
						    $H2_fourderivednucs[2]=$H2_fourderivednucs[2]/($#temp42+1);	
						    $H2_fourderivednucs[3]=$H2_fourderivednucs[3]/($#temp42+1);	
						    # first calculate the sum of all derived allele frequencies in H2
						    $H2_C=$H2_fourderivednucs[0]+$H2_fourderivednucs[1]+$H2_fourderivednucs[2]+$H2_fourderivednucs[3];
						    # now calculate the sum of the frequencies of the B derived SNP only (including ties)
						    foreach(@H3_B){
						    	$H2_B+=$H2_fourderivednucs[$_];
						    }	
						    # Now, to get the C frequency, subtract the B freq from the total.
						    $H2_C=$H2_C-$H2_B;

						    # now get the B and C frequencies from H1
						    $H1_fourderivednucs[0]=0;
						    $H1_fourderivednucs[1]=0;
						    $H1_fourderivednucs[2]=0;
						    $H1_fourderivednucs[3]=0;
						    foreach(@temp41){
						    	if(uc $_ ne uc $A){
						    		if(uc $_ eq 'A'){
						    			$H1_fourderivednucs[0]+=1;
						    		}
						    		elsif(uc $_ eq 'C'){
						    			$H1_fourderivednucs[1]+=1;
						    		}
						    		elsif(uc $_ eq 'T'){
						    			$H1_fourderivednucs[2]+=1;
						    		}
						    		elsif(uc $_ eq 'G'){
						    			$H1_fourderivednucs[3]+=1;
						    		}
						    	}
						    }	
						    $H1_fourderivednucs[0]=$H1_fourderivednucs[0]/($#temp41+1);	
						    $H1_fourderivednucs[1]=$H1_fourderivednucs[1]/($#temp41+1);	
						    $H1_fourderivednucs[2]=$H1_fourderivednucs[2]/($#temp41+1);	
						    $H1_fourderivednucs[3]=$H1_fourderivednucs[3]/($#temp41+1);	
						    # first calculate the sum of all derived allele frequencies in H1
						    $H1_C=$H1_fourderivednucs[0]+$H1_fourderivednucs[1]+$H1_fourderivednucs[2]+$H1_fourderivednucs[3];
						    # now calculate the sum of the frequencies of the B derived SNP only (including ties)
						    foreach(@H3_B){
						    	$H1_B+=$H1_fourderivednucs[$_];
						    }	
						    # Now, to get the C frequenci, subtract the B freq from the total.
						    $H1_C=$H1_C-$H1_B;

						    # I am storing the CBBA proportions per site in the ABBA hash
						    # and the BCBA proportions per site in the BABA hash

						    # I will count Bs as the frequency in H1 or H2 of the most frequent derived nucleotide in H3
						    # including all of those that are tied in H3 for the most frequent
						    # C is the sum of the frequency in H1 or H2 of all other derived alleles 
						    # for H3, a B is just the sum of the frequency of all derived nucleotides (this is
						   	# the sum of Bs and Cs in H3 and is equal to the $H3_derived_freq calculated above)

						    # Cs are counted in H1 and H2 as the sum of each frequency of each derived nucleotide in H1 or H2 times 
						    # (1- the corresponding frequency in H3)

						    #if($H1_C ne $H2_C){
						    	print $H1string," ",$H2string," ",$H3string," ",$A," ",$H1_B," ",$H1_C," ",$H2_B," ",$H2_C," ",$H3_B," ",$H3_C,"\n";
							    #print $H1_C," ",$H2_C,"\n";
							#}

							$ABBA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_C*$H2_B*$H3_derived_freq);
							#$ABBA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_C);
							
							# $ABBA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq);
							
							$BABA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_B*$H2_C*$H3_derived_freq);
							#$BABA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H2_C);

							# $BABA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq);

							$BBAA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_C*$H2_C*$H3_derived_freq);

							if($H2_derived_freq > $H3_derived_freq){
								$peak_H2_H3_derived_freq = $H2_derived_freq;
							}
							else{
								$peak_H2_H3_derived_freq = $H3_derived_freq;	
							}
							$ABBA_peak_hash=((1-$H1_derived_freq)*$peak_H2_H3_derived_freq*$peak_H2_H3_derived_freq);
							$BABA_peak_hash=($H1_derived_freq*(1-$peak_H2_H3_derived_freq)*$peak_H2_H3_derived_freq);
							# here we are calculating stats for f for H2 and H3.
							# we need to do this for each site because some of them can be undefined and should
							# not be included in the average
							if(($ABBA_peak_hash - $BABA_peak_hash) != 0){
								$f_H2_H3{$window_counter."_".$current_chromosome."_".$current_window} += ((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq))-(($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)))/ ($ABBA_peak_hash - $BABA_peak_hash);
								$f_H2_H3_counter{$window_counter."_".$current_chromosome."_".$current_window} += 1;
							}	
							# here we are calculating stats for f and the assignment of H1 and H2 needs to be switched.
							if($H1_derived_freq > $H3_derived_freq){
								$peak_H1_H3_derived_freq = $H1_derived_freq;
							}
							else{
								$peak_H1_H3_derived_freq = $H3_derived_freq;	
							}
							$BABA_peak_hashH1H3=($peak_H1_H3_derived_freq*(1-$H2_derived_freq)*$peak_H1_H3_derived_freq);
							$ABBA_peak_hashH1H3=((1-$peak_H1_H3_derived_freq)*$H2_derived_freq*$peak_H1_H3_derived_freq);
							if(($BABA_peak_hashH1H3 - $ABBA_peak_hashH1H3) != 0){
								$f_H1_H3{$window_counter."_".$current_chromosome."_".$current_window} += ((($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq))-(((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq)))/ ($BABA_peak_hashH1H3 - $ABBA_peak_hashH1H3);
								$f_H1_H3_counter{$window_counter."_".$current_chromosome."_".$current_window} += 1;
							}
							if(($H2_derived_freq >= $H1_derived_freq)&&(($ABBA_peak_hash - $BABA_peak_hash)!=0)){
								$f_dm{$window_counter."_".$current_chromosome."_".$current_window} += ((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq))-(($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)))/ ($ABBA_peak_hash - $BABA_peak_hash);
								$f_dm_counter{$window_counter."_".$current_chromosome."_".$current_window} += 1;
							}
							elsif(($ABBA_peak_hashH1H3 - $BABA_peak_hashH1H3)!=0){
								$f_dm{$window_counter."_".$current_chromosome."_".$current_window} += ((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq))-(($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)))/ -($ABBA_peak_hashH1H3 - $BABA_peak_hashH1H3);	
								$f_dm_counter{$window_counter."_".$current_chromosome."_".$current_window} += 1;
							}
	
						}
					}
				}	
			}
		}
	}
}


close DATAINPUT;

# now merge the hash keys
my @common_keys = ();

foreach (keys %ABBA_hash) {
	push(@common_keys, $_);
}

foreach (keys %BABA_hash) {
	push(@common_keys, $_) unless exists $ABBA_hash{$_};
}

my @common_keys2 = ();

foreach (keys %ABBA_hash) {
	push(@common_keys2, $_);
}

foreach (keys %BBAA_hash) {
	push(@common_keys2, $_) unless exists $ABBA_hash{$_};
}


my @out = keys %{{map {($_ => 1)} (@common_keys, @common_keys2)}};

@out = map  { $_->[0] }
             sort { $a->[1] <=> $b->[1] }
             map  { [$_, $_=~/(\d+)/] }
                 @out;



foreach (@out) {
	@temp1=split('_',$_);
	print OUTFILE $temp1[1],"\t",$temp1[2]+1,"\t",$temp1[2]+$sliding_window,"\t";
	if(defined($ABBA_hash{$_})){
		print OUTFILE $ABBA_hash{$_},"\t";
	}
	else{
		print OUTFILE "0\t";
	}
	if(defined($BABA_hash{$_})){
		print OUTFILE $BABA_hash{$_},"\t0\t0\t0\t0\n";
	}
	else{
		print OUTFILE "0\t0\t0\t0\t0\n";
	}
}
if($#out == -1){
	print OUTFILE "chr0\t1\t500000\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
}

close OUTFILE;


print OUTFILE2 "chromosome\tbegin\tend\tCBBA\tBCBA\tCCBA\tD\tfdH2H3\tfH1H3\tf_dm\tdH2H3\tnum_sites_per_windowH2H3\tdH1H3\tnum_sites_per_windowH1H3\tH2pi\tnumsitesH2pi\tH1pi\tnumsitesH1pi\n";
foreach (@out) {
	@temp1=split('_',$_);
	print OUTFILE2 $temp1[1],"\t",$temp1[2]+1,"\t",$temp1[2]+$sliding_window,"\t";
	if(defined($ABBA_hash{$_})){
		print OUTFILE2 $ABBA_hash{$_},"\t";
	}
	else{
		print OUTFILE2 "0\t";
	}
	if(defined($BABA_hash{$_})){
		print OUTFILE2 $BABA_hash{$_},"\t";
	}
	else{
		print OUTFILE2 "0\t";
	}
	if(defined($BBAA_hash{$_})){
		print OUTFILE2 $BBAA_hash{$_},"\t";
	}
	else{
		print OUTFILE2 "0\t";
	}
	#print D for this window
	if((defined($ABBA_hash{$_}))&&(defined($BBAA_hash{$_}))){
		if(($ABBA_hash{$_}+$BABA_hash{$_})>0){
			print OUTFILE2 ($ABBA_hash{$_}-$BABA_hash{$_}) / ($ABBA_hash{$_}+$BABA_hash{$_}),"\t";
		}
		else{
			print OUTFILE2 "NAN\t";
		}
	}
	else{
		print OUTFILE2 "NAN\t";
	}

	#print fd H2H3 for this window
	if(defined($f_H2_H3_counter{$_})){
		print OUTFILE2 $f_H2_H3{$_}/$f_H2_H3_counter{$_},"\t";
	}
	else{
		print OUTFILE2 "NAN\t";
	}
	#print fd H1H3 for this window
	if(defined($f_H1_H3_counter{$_})){
		print OUTFILE2 $f_H1_H3{$_}/$f_H1_H3_counter{$_},"\t";
	}
	else{
		print OUTFILE2 "NAN\t";
	}	
	#print f_dm for this window
	if(defined($f_dm{$_})){
		print OUTFILE2 $f_dm{$_}/$f_dm_counter{$_},"\t";
	}
	else{
		print OUTFILE2 "NAN\t";
	}	
	#print average H2H3 pairwise divergence for this window and number of H2H3 sites
	if(defined($number_of_sites_per_window{$_})){
		if($number_of_sites_per_window{$_}>0){
			print OUTFILE2 ($H2_H3_pairwise_divergence_per_window{$_}/$number_of_sites_per_window{$_}),"\t",$number_of_sites_per_window{$_},"\t";
		}
		else{
			print OUTFILE2 "NAN\t",$number_of_sites_per_window{$_},"\t";
		}
	}
	else{
		print OUTFILE2 "NAN\t",$number_of_sites_per_window{$_},"\t";
	}
	#print average H1H3 pairwise divergence for this window and number of H1H3 sites
	if(defined($number_of_sites_per_windowH1H3{$_})){
		if($number_of_sites_per_windowH1H3{$_} > 0){
			print OUTFILE2 ($H1_H3_pairwise_divergence_per_window{$_}/$number_of_sites_per_windowH1H3{$_}),"\t",$number_of_sites_per_windowH1H3{$_},"\t";
		}
		else{
			print OUTFILE2 "NAN\t",$number_of_sites_per_windowH1H3{$_},"\t";
		}
	}
	else{
		print OUTFILE2 "NAN\t",$number_of_sites_per_windowH1H3{$_},"\t";
	}
	#print H2 pairwise nucleotide diversity per site for this window and number of H2 pairwise nucleotide diversity sites
	if(defined($H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_})){
		if($H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}>0){
			print OUTFILE2 ($H2_pairwise_nucleotide_diversity_per_window{$_}/$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}),"\t",$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\t";
		}
		else{
			print OUTFILE2 "NAN\t",$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\t";
		}
	}
	else{
		print OUTFILE2 "NAN\t",$H2_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\t";
	}
	#print H1 pairwise nucleotide diversity per site for this window and number of H1 pairwise nucleotide diversity sites
	if(defined($H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_})){
		if($H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}>0){
			print OUTFILE2 ($H1_pairwise_nucleotide_diversity_per_window{$_}/$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_}),"\t",$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\n";
		}
		else{
			print OUTFILE2 "NAN\t",$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\n";
		}
	}
	else{
		print OUTFILE2 "NAN\t",$H1_number_of_sites_for_pairwise_nucleotide_diversity_per_window{$_},"\n";
	}
}
if($#out == -1){
	print OUTFILE2 "chr0\t1\t500000\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
}

close OUTFILE2;

```
