# Performing ABBABABA on populations

I have learned lots more about ABBABABA tests by reading the cool paper by Martin, Davey, and Jiggins in MBE.  I wrote a script to perform the ABBABABA test on population samples and also calculate the f_d statistic and average pairwise divergence between H2 and H3 as discussed in that paper.  I called it "Performs_ABBA_BABA_on_populations.pl":

```perl 
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;


#  This program reads in a tab delimited genotype file generated
#  by the perl program "Gets_outgroup_sequence_from_axt_files.pl"
#  or from vcftools
#  and performs the population ABBA-BABA test using one outgroup
#  sequence and one or more individual sequences from three other
#  taxa (H3, H2, and H1). 

#  This analysis will include all positions that have data from at least 
#  one individual from species H3, H2, and H1, including those that that 
#  missing data in some individuals.

# to execute type Performs_ABBA_BABA_on_popoulations.pl inputfile.tab 1111100110000111100011100110010100000000 
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

# born
# 14-18-19-20
# sum
# 15-16-17-21
# tonk
# 32-33-34-35-36-37-38-39-40
# heck
# 2-3-4-5-6-7_
# maura
# 8-9-10-11-12-13
# nigrescens
# 24-26
# nigra
# 22-23-25
# togeanus
# 30-31
# och_brun
# 27-28-29-1


# for example, with a tab file with the rhesus, baboon and human sequence , here is the input command using the rhesus as the outgroup:

# tonk_heck
# perl Performs_ABBA_BABA_on_populations.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 3_6_14-18-19-20_2-3-4-5-6-7_32-33-34-35-36-37-38-39-40 heck_tonk_born_pop.abbababa
# perl Performs_ABBA_BABA_on_populations.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 3_6_14-18-19-20_32-33-34-35-36-37-38-39-40_2-3-4-5-6-7 tonk_heck_born_pop.abbababa

# nigra_och/brun
# perl Performs_ABBA_BABA_on_populations.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 3_6_14-18-19-20_22-23-25_27-28-29-1 nigra_ochbrun_borneo_pop.abbababa
# perl Performs_ABBA_BABA_on_populations.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 3_6_15-16-17-21_22-23-25_27-28-29-1 nigra_ochbrun_sum_pop.abbababa

# heck_maura
# perl Performs_ABBA_BABA_on_populations.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 3_6_14-18-19-20_2-3-4-5-6-7_8-9-10-11-12-13 hecki_maura_borneo_pop.abbababa
# perl Performs_ABBA_BABA_on_populations.pl final_round2_filtered.vcf.gz_with_baboon_and_human.tab 1111100110000111100011100110010100000000 3_6_15-16-17-21_2-3-4-5-6-7_8-9-10-11-12-13 hecki_maura_sum_pop.abbababa



my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $input3 = $ARGV[2];
my $outputfile = $ARGV[3];
my $outputfile2 = $ARGV[4];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

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


my @temp;
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

my $number_of_female_individuals_genotyped;
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
my %ABBA_peak_hash;
my %BABA_hash;
my %BABA_peak_hash;
my %BBAA_hash;
my $A;
my $B;
my @allelez;
my $derived;
my $ancestral;
my $H3_derived_freq;
my $H1_derived_freq;
my $H2_derived_freq;
my $peak_H2_H3_derived_freq;
my @H3allelez;
my @H2allelez;
my $diff;
my $num_comparisons;
my %H2_H3_pairwise_divergence_per_window;
my %number_of_sites_per_window;


while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] eq '#CHROM'){
		print "The outgroup sequence is ",$temp[$whotoinclude[0]-1],"\n";
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
		if($temp[0] ne $current_chromosome){
			$current_chromosome = $temp[0];
			$current_window = 0;
			$window_counter+=1;
		}
		until($temp[1] < ($current_window+$sliding_window)){
			$current_window = $current_window+$sliding_window;
			$window_counter+=1;
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
				$derived=0;
				$ancestral=0;
				for ($y = 0 ; $y <= $#H3; $y++ ) {
					if($temp[$H3[$y]+$whotoinclude[1]-2] ne './.'){
						@allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						if($allelez[1] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[1];
					}	
				}
				if(($derived+$ancestral)>0){
					$H3_derived_freq=$derived/($derived+$ancestral);
				}
				# now calculate the frequency of the derived allele in the H1 data
				$H1_derived_freq=0;
				$derived=0;
				$ancestral=0;
				for ($y = 0 ; $y <= $#H1; $y++ ) {
					if($temp[$H1[$y]+$whotoinclude[1]-2] ne './.'){
						@allelez=split('/',$temp[$H1[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						if($allelez[1] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[1];
					}
				}
				if(($derived+$ancestral)>0){
					$H1_derived_freq=$derived/($derived+$ancestral);
				}
				# now calculate the frequency of the derived allele in the H2 data
				$H2_derived_freq=0;
				$derived=0;
				$ancestral=0;
				for ($y = 0 ; $y <= $#H2; $y++ ) {
					if($temp[$H2[$y]+$whotoinclude[1]-2] ne './.'){
						@allelez=split('/',$temp[$H2[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						if($allelez[1] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[1];
					}	
				}
				if(($derived+$ancestral)>0){
					$H2_derived_freq=$derived/($derived+$ancestral);
				}
				# Now calculate the average pairwise divergence between H2 and H3
				$diff=0;
				$num_comparisons=0;
				@H3allelez=();
				@H2allelez=();
				for ($y = 0 ; $y <= $#H3; $y++ ) {
					if($temp[$H3[$y]+$whotoinclude[1]-2] ne './.'){
						for ($w = 0 ; $w <= $#H2; $w++ ) {
							if($temp[$H2[$w]+$whotoinclude[1]-2] ne './.'){
								# there are data for H2 and H3
								@H3allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
								@H2allelez=split('/',$temp[$H3[$w]+$whotoinclude[1]-2]);
							}
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
				# now calculate the average pairwise divergence
				if($num_comparisons>0){
					if((defined($string))&&($#H1>=0)&&($#H2>=0)&&($#H3>=0)){
						$number_of_sites_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=1;
						$H2_H3_pairwise_divergence_per_window{$window_counter."_".$current_chromosome."_".$current_window}+=($diff/$num_comparisons);
						# later we will standardize this by the number of sites per window
						# this is important because non-diverged sites are an important part of divergence.	
					}
				}
				# Now calculate the ABBA_BABA site weighting
				if((defined($string))&&($#H1>=0)&&($#H2>=0)&&($#H3>=0)){
					@temp1=split('',$string);
					$x_uniq = uniq @temp1;
					if($x_uniq == 2){
						if(
							(((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq)>0)
						||
							(($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)>0)
						||
							(($H1_derived_freq*$H2_derived_freq*(1-$H3_derived_freq))>0)
							){
							# this is a polymorphic position and it is an ABBA_BABA site.  
								$ABBA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq);
								$BABA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq);
								$BBAA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_derived_freq*$H2_derived_freq*(1-$H3_derived_freq));
								if($H2_derived_freq > $H3_derived_freq){
									$peak_H2_H3_derived_freq = $H2_derived_freq;
								}
								else{
									$peak_H2_H3_derived_freq = $H3_derived_freq;	
								}
								$ABBA_peak_hash{$window_counter."_".$current_chromosome."_".$current_window}+=((1-$H1_derived_freq)*$peak_H2_H3_derived_freq*$peak_H2_H3_derived_freq);
								$BABA_peak_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_derived_freq*(1-$peak_H2_H3_derived_freq)*$peak_H2_H3_derived_freq);

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

@common_keys = map  { $_->[0] }
             sort { $a->[1] <=> $b->[1] }
             map  { [$_, $_=~/(\d+)/] }
                 @common_keys;
#@common_keys = sort @common_keys;


foreach (@common_keys) {
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
if($#common_keys == -1){
	print OUTFILE "chr0\t1\t500000\t0\t0\t0\t0\t0\t0\n";
}

close OUTFILE;


print OUTFILE2 "chromosome\tbegin\tend\tABBA\tBABA\tBBAA\tD\tfd\tdxy\tnum_sites_per_window\n";
foreach (@common_keys) {
	@temp1=split('_',$_);
	print OUTFILE2 $temp1[1],"\t",$temp1[2]+1,"\t",$temp1[2]+$sliding_window,"\t";
	if(defined($ABBA_hash{$_})){
		print OUTFILE2 $ABBA_hash{$_},"\t";
	}
	else{
		print OUTFILE2 "0\t";
	}
	if(defined($BABA_hash{$_})){
		print OUTFILE2 $BABA_hash{$_},"\t0\t0\t0\t0\t";
	}
	else{
		print OUTFILE2 "0\t0\t0\t0\t0\t";
	}
	if(defined($BBAA_hash{$_})){
		print OUTFILE2 $BBAA_hash{$_},"\t";
	}
	else{
		print OUTFILE2 "0\t";
	}
	#print D for this window
	if(($ABBA_hash{$_}+$BABA_hash{$_})>0){
		print OUTFILE2 ($ABBA_hash{$_}-$BABA_hash{$_})/($ABBA_hash{$_}+$BABA_hash{$_}),"\t";
	}
	else{
		print OUTFILE2 "NAN\t;"
	}

	#print fd for this window
	if(($ABBA_peak_hash{$_}-$BABA_peak_hash{$_})>0){
		print OUTFILE2 ($ABBA_hash{$_}-$BABA_hash{$_})/($ABBA_peak_hash{$_}-$BABA_peak_hash{$_}),"\t";
	}
	else{
		print OUTFILE2 "NAN\t;"
	}
	#print average pairwise divergence for this window
	if($number_of_sites_per_window{$_}>0){
		print OUTFILE2 ($H2_H3_pairwise_divergence_per_window{$_}/$number_of_sites_per_window{$_}),"\t",$number_of_sites_per_window{$_},"\n";
	}
	else{
		print OUTFILE2 "NAN\n;"
	}
}
if($#common_keys == -1){
	print OUTFILE2 "chr0\t1\t500000\t0\t0\t0\t0\t0\t0\t0\n";
}

close OUTFILE2;
```

This script generates two files - one is analyzed by the R script from ANGSD called 'jacknife.R' and the other can be fed into an as yet not written R script to analyze cool stuff in windows that is in the file 'stats'.

And I wrote a wrapper that runs this program with each of the Sulawesi macaques as H2, the other Sulawesi macaques as H1, and nemestrina as H3.  I called this "Wrapper_for_Performs_ABBA_BABA_populations_most_inclusive.pl":

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;


#  This program reads in a tab delimited genotype file generated
#  by the perl program "Gets_outgroup_sequence_from_axt_files.pl"
#  or from vcftools
#  and performs the population ABBA-BABA test using one outgroup
#  sequence and one or more individual sequences from three other
#  taxa (H3, H2, and H1). 

#  This analysis will include all positions that have data from at least 
#  one individual from species H3, H2, and H1, including those that that 
#  missing data in some individuals.

# to execute type Performs_ABBA_BABA_on_popoulations.pl inputfile.tab 1111100110000111100011100110010100000000 
# 3_6_14-18-19-20_2-3-4-5-6-7_32-33-34-35-36-37-38-39-40 output.txt  
# where 1111100110000111100011100110010100000000 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) or is not (0) female ,and 3_6 refers to (i) the column that contains the 
# outgroup nucleotide (3 in this case for rhesus) and (ii) the column number of the first individual in the ingroup 
# (6 in this case), and 14-18-19-20, 2-3-4-5-6-7, and 32-33-34-35-36-37-38-39-40 refer to H3, H1, and H2 samples as 
# itemized below (here, they are Borneo nemestrina, hecki, and tonkeana, respectively).

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

# born
# 14-18-19-20
# sum
# 15-16-17-21
# tonk
# 32-33-34-35-36-37-38-39-40
# heck
# 2-3-4-5-6-7_
# maura
# 8-9-10-11-12-13
# nigrescens
# 24-26
# nigra
# 22-23-25
# togeanus
# 30-31
# och_brun
# 27-28-29-1

my @borneo = ("nem_Gumgum_stampy_sorted","nem_PM664_stampy_sorted","nem_PM665_stampy_sorted","nem_Sukai_male_stampy_sorted");
my $borneo_numbers = "14-18-19-20";
my @sumatra = ("nem_Kedurang_stampy_sorted","nem_Malay_stampy_sorted","nem_Ngasang_stampy_sorted");
my @sumatra_numbers=("15-16-17");
my $pagensis = "nem_pagensis_stampy_sorted";
my $pagensis_number=21;
my @sumatra_pagensis_numbers=("15-16-17-21");

#my @H3_names=("borneo","sumatra","sumatra_pagensis","pagensis");
#my @H3_numbers=("14-18-19-20","15-16-17","15-16-17-21","21");
my @H3_names=("nem");
my @H3_numbers=("14-15-16-17-18-19-20-21");

my @sulawesi = ("brun","hecki","maura","nigra","nigrescens","ochreata","togeanus","tonk");
my @sulawesi_numbers=("1","2-3-4-5-6-7","8-9-10-11-12-13","22-23-25","24-26","27-28-29","30-31","32-33-34-35-36-37-38-39-40");

my @H1_numbers;
my $H1_name;
my $H1_sample_list;

my $sliding_window=5000000;
my $outgroup_number=3; # for the tab file with human and baboon outgroup, 3=rhesus, 4=human, and 5=baboon
my $ingroup_column_begin_number=6; # for the tab file with human and baboon outgroup, this is 6
my $infile_tab = "final_round2_filtered.vcf.gz_with_baboon_and_human.tab";
my $commandline;
my $status;
my $x;
my $y;
my $z;



for ($x = 0 ; $x <= $#sulawesi_numbers ; $x++ ) { # This is the H2 Sulawesi sample
	# first make the H1 name_and_sample list
	$H1_name="";
	$H1_sample_list="";
	@H1_numbers=();	
	for ($z = 0; $z <= $#sulawesi_numbers ; $z++ ) {
		if($z != $x){
			push(@H1_numbers, $z);
		}
	}	
	for ($z = 0; $z < $#H1_numbers ; $z++ ) { # This is the H1 Sulawesi sample
		$H1_sample_list = $H1_sample_list.$sulawesi_numbers[$H1_numbers[$z]]."-";
		$H1_name = $H1_name.$sulawesi[$H1_numbers[$z]]."-";
	}
	$H1_sample_list = $H1_sample_list.$sulawesi_numbers[$H1_numbers[$#H1_numbers]];
	$H1_name = $H1_name.$sulawesi[$H1_numbers[$#H1_numbers]];
	print "H1_name ",$H1_name,"\n";
	print "H1_sample_list ",$H1_sample_list,"\n";
	print "H2_name ",$sulawesi_numbers[$x],"\n";
	print "H2_sample_list ",$sulawesi[$x],"\n";

	for ($y = 0 ; $y <= $#H3_numbers ; $y++ ) { # This is the nemestrina samples
			$commandline = "perl Performs_ABBA_BABA_on_populations.pl ".$infile_tab." 1111100110000111100011100110010100000000 ".$outgroup_number."_".$ingroup_column_begin_number."_".$H3_numbers[$y]."_";
			$commandline = $commandline.$H1_sample_list."_".$sulawesi_numbers[$x]." ".$sulawesi_numbers[$x]."_".$H1_sample_list."_".$H3_numbers[$y].".abbababa ".$sulawesi_numbers[$x]."_".$H1_sample_list."_".$H3_numbers[$y].".stats";
			print $commandline,"\n";
			$status = system($commandline);
			# now make a file with the names of the taxa
			my $outputfile2 = $sulawesi_numbers[$x]."_".$H1_sample_list."_".$H3_numbers[$y].".names";
			unless (open(OUTFILE2, ">$outputfile2"))  {
				print "I can\'t write to $outputfile2\n";
				exit;
			}
			print "Creating output file: $outputfile2\n";
			print OUTFILE2 $H3_names[$y],"\n",$sulawesi[$x],"\n",$H1_name,"\n";
			close OUTFILE2;
			$commandline = "Rscript jackKnife.R file=".$sulawesi_numbers[$x]."_".$H1_sample_list."_".$H3_numbers[$y].".abbababa indNames=".$sulawesi_numbers[$x]."_".$H1_sample_list."_".$H3_numbers[$y].".names outfile=".$sulawesi_numbers[$x]."_".$H1_sample_list."_".$H3_numbers[$y].".out";
			print $commandline,"\n";
			$status = system($commandline);
	}
}

```
