I have not been good with documenting the stuff I have done lately.  Here is another effort.

First, I have made some new Wrapper scripts to generate f_dm stat files using each Sulawesi species as H2 and all other Sulawesi species as H1, and only the Borneo samples as H3.  This needs to be done separately for the autosomes and the X.  Here is the autosomal one (Wrapper_for_Performs_ABBA_BABA_Borneo_vs_allsulawesi_populations_mirror.pl):

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

my @H3_names=("borneo");
my @H3_numbers=("14-18-19-20");
#my @sulawesi = ("brun","hecki","maura","nigra","nigrescens","ochreata","togeanus","tonk");
#my @sulawesi_numbers=("1","2-3-4-5-6-7","8-9-10-11-12-13","22-23-25","24-26","27-28-29","30-31","32-33-34-35-36-37-38-39-40");

my @sulawesi = ("tonk","togeanus","ochreata","nigrescens","nigra","maura","hecki","brun");
my @sulawesi_numbers=("32-33-34-35-36-37-38-39-40","30-31","27-28-29","24-26","22-23-25","8-9-10-11-12-13","2-3-4-5-6-7","1");

my $sliding_window=5000000;
my $outgroup_number=3; # for the tab file with human and baboon outgroup, 3=rhesus, 4=human, and 5=baboon
my $ingroup_column_begin_number=6; # for the tab file with human and baboon outgroup, this is 6
my $infile_tab = "final_round2_filtered.vcf.gz_with_baboon_and_human.tab";
my $commandline;
my $status;
my $x;
my $y;
my $z;

my $Hi_string;


for ($x = 0 ; $x <= $#sulawesi_numbers ; $x++ ) { # This is the first (H1) Sulawesi sample
	for ($y = 0 ; $y <= $#H3_numbers ; $y++ ) { # This is each of the nem populations
		$commandline = "perl Performs_ABBA_BABA_on_populations.pl ".$infile_tab." 1111100110000111100011100110010100000000 ".$outgroup_number."_".$ingroup_column_begin_number."_".$H3_numbers[$y]."_";
		# this is H1
		$Hi_string=();
		for ($z = 0 ; $z <= $#sulawesi_numbers ; $z++ ) { # This is the second (H2) Sulawesi sample
			if(($z != $x)&&($z != $#sulawesi_numbers)&&($x != $#sulawesi_numbers)){
				$commandline = $commandline.$sulawesi_numbers[$z]."-";
				$Hi_string=$Hi_string.$sulawesi_numbers[$z]."-";
			}	
			elsif(($z != $x)&&($z == $#sulawesi_numbers)){
				$commandline = $commandline.$sulawesi_numbers[$z];
				$Hi_string=$Hi_string.$sulawesi_numbers[$z];
			}
			elsif(($z != $x)&&($x == $#sulawesi_numbers)){
				if($z < ($#sulawesi_numbers-1)){
					$commandline = $commandline.$sulawesi_numbers[$z]."-";
					$Hi_string=$Hi_string.$sulawesi_numbers[$z]."-";
				}
				elsif($z == ($#sulawesi_numbers-1)){
					$commandline = $commandline.$sulawesi_numbers[$z];
					$Hi_string=$Hi_string.$sulawesi_numbers[$z];
				}
			}
		}
		# now add H2
		$commandline = $commandline."_".$sulawesi_numbers[$x]." ".$Hi_string."_".$sulawesi_numbers[$x]."_".$H3_numbers[$y].".abbababa"." ".$Hi_string."_".$sulawesi_numbers[$x]."_".$H3_numbers[$y].".stats";
			
		print $commandline,"\n";
		$status = system($commandline);
		# now make a file with the names of the taxa
		my $outputfile2 = $Hi_string."_".$sulawesi_numbers[$x]."_".$H3_numbers[$y].".names";
		unless (open(OUTFILE2, ">$outputfile2"))  {
			print "I can\'t write to $outputfile2\n";
			exit;
		}
		print "Creating output file: $outputfile2\n";
		print OUTFILE2 $H3_names[$y],"\n",$Hi_string,"\n",$sulawesi[$x],"\n";
		close OUTFILE2;
		$commandline = "Rscript jackKnife.R file=".$Hi_string."_".$sulawesi_numbers[$x]."_".$H3_numbers[$y].".abbababa indNames=".$Hi_string."_".$sulawesi_numbers[$x]."_".$H3_numbers[$y].".names outfile=".$Hi_string."_".$sulawesi_numbers[$x]."_".$H3_numbers[$y].".out";
		print $commandline,"\n";
		$status = system($commandline);
	}
}

```

and here is the chrX one (Wrapper_for_Performs_ABBA_BABA_Borneo_vs_allsulawesi_populations_onlychrX_mirror.pl) DONT USE THIS BECAUSE IT DOESNT HAVE THE CORRECT CHRX SPECIFIC ABBABABA SCRIPT, WHICH IS NOW BELOW:

``` perl
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

my @H3_names=("borneo","sumatra","sumatra_pagensis","pagensis");
my @H3_numbers=("14-18-19-20","15-16-17","15-16-17-21","21");
#my @sulawesi = ("brun","hecki","maura","nigra","nigrescens","ochreata","togeanus","tonk");
#my @sulawesi_numbers=("1","2-3-4-5-6-7","8-9-10-11-12-13","22-23-25","24-26","27-28-29","30-31","32-33-34-35-36-37-38-39-40");

my @sulawesi = ("tonk","togeanus","ochreata","nigrescens","nigra","maura","hecki","brun");
my @sulawesi_numbers=("32-33-34-35-36-37-38-39-40","30-31","27-28-29","24-26","22-23-25","8-9-10-11-12-13","2-3-4-5-6-7","1");

my $sliding_window=5000000;
my $outgroup_number=3; # for the tab file with human and baboon outgroup, 3=rhesus, 4=human, and 5=baboon
my $ingroup_column_begin_number=6; # for the tab file with human and baboon outgroup, this is 6
my $infile_tab = "final_round2_filtered.vcf.gz_with_baboon_and_human.tab";
my $commandline;
my $status;
my $x;
my $y;
my $z;

my $Hi_string;

for ($x = 0 ; $x <= $#sulawesi_numbers ; $x++ ) { # This is the first (H1) Sulawesi sample
	for ($y = 0 ; $y <= $#H3_numbers ; $y++ ) { # This is each of the nem populations
		$commandline = "perl Performs_ABBA_BABA_on_populations_onlychrX.pl ".$infile_tab." 1111100110000111100011100110010100000000 ".$outgroup_number."_".$ingroup_column_begin_number."_".$H3_numbers[$y]."_";
		# this is H1
		$Hi_string="";
		for ($z = 0 ; $z <= $#sulawesi_numbers ; $z++ ) { # This is the second (H2) Sulawesi sample
			if(($z != $x)&&($z != $#sulawesi_numbers)&&($x != $#sulawesi_numbers)){
				$commandline = $commandline.$sulawesi_numbers[$z]."-";
				$Hi_string=$Hi_string.$sulawesi_numbers[$z]."-";
			}	
			elsif(($z != $x)&&($z == $#sulawesi_numbers)){
				$commandline = $commandline.$sulawesi_numbers[$z];
				$Hi_string=$Hi_string.$sulawesi_numbers[$z];
			}
			elsif(($z != $x)&&($x == $#sulawesi_numbers)){
				if($z < ($#sulawesi_numbers-1)){
					$commandline = $commandline.$sulawesi_numbers[$z]."-";
					$Hi_string=$Hi_string.$sulawesi_numbers[$z]."-";
				}
				elsif($z == ($#sulawesi_numbers-1)){
					$commandline = $commandline.$sulawesi_numbers[$z];
					$Hi_string=$Hi_string.$sulawesi_numbers[$z];
				}
			}
		}
		# now add H2
		$commandline = $commandline."_".$sulawesi_numbers[$x]." ".$Hi_string."_".$sulawesi_numbers[$x]."_".$H3_numbers[$y].".onlychrX_abbababa"." ".$Hi_string."_".$sulawesi_numbers[$x]."_".$H3_numbers[$y].".onlychrX_stats";
			
		print $commandline,"\n";
		$status = system($commandline);
		# now make a file with the names of the taxa
		my $outputfile2 = $Hi_string."_".$sulawesi_numbers[$x]."_".$H3_numbers[$y].".onlychrX_names";
		unless (open(OUTFILE2, ">$outputfile2"))  {
			print "I can\'t write to $outputfile2\n";
			exit;
		}
		print "Creating output file: $outputfile2\n";
		print OUTFILE2 $H3_names[$y],"\n",$Hi_string,"\n",$sulawesi[$x],"\n";
		close OUTFILE2;
		$commandline = "Rscript jackKnife.R file=".$Hi_string."_".$sulawesi_numbers[$x]."_".$H3_numbers[$y].".onlychrX_abbababa indNames=".$Hi_string."_".$sulawesi_numbers[$x]."_".$H3_numbers[$y].".onlychrX_names outfile=".$Hi_string."_".$sulawesi_numbers[$x]."_".$H3_numbers[$y].".onlychrX_out";
		print $commandline,"\n";
		$status = system($commandline);
	}
}

```

I modified the `Performs_ABBA_BABA_on_populations.pl` script to work on chrX (`Performs_ABBA_BABA_on_populations_onlychrX.pl`).  This has several important differences and should accomodate chrX of males being represented like this `T/` or like this `T/T` even though the latter indicates a hemizygous genotype.  Here is the script:

```
#!/usr/bin/env perl
use strict;
use warnings;
#use lib qw(/home/ben/perl_modules/Number-Range-0.12/lib/);
use List::MoreUtils qw/ uniq /;
#use Number::Range;


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
				for ($y = 0 ; $y <= $#H3; $y++ ) {
					if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')&&($temp[$H3[$y]+$whotoinclude[1]-2] ne './')&&
						($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/')&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/')&&
						($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')&&
						($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')&&
						(
						((length($temp[$H3[$y]+$whotoinclude[1]-2]) == 3)&&($sexes[$H3[$y]-1] == 1))
							||
						((length($temp[$H3[$y]+$whotoinclude[1]-2]) == 2)&&($sexes[$H3[$y]-1] == 0))
						)
						){						
						@allelez=split('/',$temp[$H3[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						# check if this individual is a female before considering the second X allele
						if($sexes[$H3[$y]-1] == 1){
							if($allelez[1] eq $A){
								$ancestral+=1;
							}
							else{
								$derived+=1;
							}
							$string=$string.$allelez[1];
						}
					}	
				}
				if(($derived+$ancestral)>0){
					$H3_derived_freq=$derived/($derived+$ancestral);
					$H3_ancestral_freq=$ancestral/($derived+$ancestral);
				}
				# now calculate the frequency of the derived allele in the H1 data
				$H1_derived_freq=0;
				$H1_ancestral_freq=0;
				$derived=0;
				$ancestral=0;
				for ($y = 0 ; $y <= $#H1; $y++ ) {
					if(($temp[$H1[$y]+$whotoinclude[1]-2] ne './.')&&($temp[$H1[$y]+$whotoinclude[1]-2] ne './')&&
						($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/')&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/')&&
						($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/.')&&
						($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/N')&&
						(
						((length($temp[$H1[$y]+$whotoinclude[1]-2]) == 3)&&($sexes[$H1[$y]-1] == 1))
							||
						((length($temp[$H1[$y]+$whotoinclude[1]-2]) == 2)&&($sexes[$H1[$y]-1] == 0))
						)
						){
						@allelez=split('/',$temp[$H1[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						# check if this individual is a female before considering the second X allele
						if($sexes[$H1[$y]-1] == 1){
							if($allelez[1] eq $A){
								$ancestral+=1;
							}
							else{
								$derived+=1;
							}
							$string=$string.$allelez[1];
						}	
					}
				}
				if(($derived+$ancestral)>0){
					$H1_derived_freq=$derived/($derived+$ancestral);
					$H1_ancestral_freq=$ancestral/($derived+$ancestral);
				}
				# now calculate the frequency of the derived allele in the H2 data
				$H2_derived_freq=0;
				$H2_ancestral_freq=0;
				$derived=0;
				$ancestral=0;
				for ($y = 0 ; $y <= $#H2; $y++ ) {
					if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')&&($temp[$H2[$y]+$whotoinclude[1]-2] ne './')&&
						($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/')&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/')&&
						($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')&&
						($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')&&
						(
						((length($temp[$H2[$y]+$whotoinclude[1]-2]) == 3)&&($sexes[$H2[$y]-1] == 1))
							||
						((length($temp[$H2[$y]+$whotoinclude[1]-2]) == 2)&&($sexes[$H2[$y]-1] == 0))
						)
						){
						@allelez=split('/',$temp[$H2[$y]+$whotoinclude[1]-2]);
						if($allelez[0] eq $A){
							$ancestral+=1;
						}
						else{
							$derived+=1;
						}
						$string=$string.$allelez[0];
						# check if this individual is a female before considering the second X allele
						if($sexes[$H2[$y]-1] == 1){
							if($allelez[1] eq $A){
								$ancestral+=1;
							}
							else{
								$derived+=1;
							}
							$string=$string.$allelez[1];
						}	
					}	
				}
				if(($derived+$ancestral)>0){
					$H2_derived_freq=$derived/($derived+$ancestral);
					$H2_ancestral_freq=$ancestral/($derived+$ancestral);
				}
				# only consider sites for which there are data for H1, H2, and H3
				if((defined($string))&&(($H1_derived_freq>0)||($H1_ancestral_freq>0))&&(($H2_derived_freq>0)||($H2_ancestral_freq>0))&&(($H3_derived_freq>0)||($H3_ancestral_freq>0))){
					# Now calculate the average pairwise divergence between H2 and H3
					$diff=0;
					$num_comparisons=0;
					@H3allelez=();
					@H2allelez=();
					for ($y = 0 ; $y <= $#H3; $y++ ) {
						if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')&&($temp[$H3[$y]+$whotoinclude[1]-2] ne './')&&
							($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/')&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/')&&
							($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')&&
							($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')&&
							(
							((length($temp[$H3[$y]+$whotoinclude[1]-2]) == 3)&&($sexes[$H3[$y]-1] == 1))
								||
							((length($temp[$H3[$y]+$whotoinclude[1]-2]) == 2)&&($sexes[$H3[$y]-1] == 0))
							)
							){
							for ($w = 0 ; $w <= $#H2; $w++ ) {
								if(($temp[$H2[$w]+$whotoinclude[1]-2] ne './.')&&($temp[$H2[$w]+$whotoinclude[1]-2] ne './')&&
									($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/')&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/')&&
									($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/.')&&
									($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/N')&&
									(
									((length($temp[$H2[$w]+$whotoinclude[1]-2]) == 3)&&($sexes[$H2[$w]-1] == 1))
										||
									((length($temp[$H2[$w]+$whotoinclude[1]-2]) == 2)&&($sexes[$H2[$w]-1] == 0))
									)
									){
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
										$num_comparisons+=1;
										# check if this individual is a female before considering the second X allele
										if($sexes[$H2[$w]-1] == 1){
											if($_ ne $H2allelez[1]){
												$diff+=1;
											}
											$num_comparisons+=1;
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
						if(($temp[$H3[$y]+$whotoinclude[1]-2] ne './.')&&($temp[$H3[$y]+$whotoinclude[1]-2] ne './')&&
							($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/')&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/')&&
							($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/.')&&
							($temp[$H3[$y]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H3[$y]+$whotoinclude[1]-2] ne 'N/N')&&
							(
							((length($temp[$H3[$y]+$whotoinclude[1]-2]) == 3)&&($sexes[$H3[$y]-1] == 1))
								||
							((length($temp[$H3[$y]+$whotoinclude[1]-2]) == 2)&&($sexes[$H3[$y]-1] == 0))
							)
							){
							for ($w = 0 ; $w <= $#H1; $w++ ) {
								if(($temp[$H1[$w]+$whotoinclude[1]-2] ne './.')&&($temp[$H1[$w]+$whotoinclude[1]-2] ne './')&&
									($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/')&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/')&&
									($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/.')&&
									($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/N')&&
									(
									((length($temp[$H1[$w]+$whotoinclude[1]-2]) == 3)&&($sexes[$H1[$w]-1] == 1))
										||
									((length($temp[$H1[$w]+$whotoinclude[1]-2]) == 2)&&($sexes[$H1[$w]-1] == 0))
									)
									){
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
										$num_comparisonsH1H3+=1;
										# check if this individual is a female before considering the second X allele
										if($sexes[$H1[$w]-1] == 1){
											if($_ ne $H1allelez[1]){
												$diffH1H3+=1;
											}
											$num_comparisonsH1H3+=1;
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
							if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')&&($temp[$H2[$y]+$whotoinclude[1]-2] ne './')&&
								($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/')&&
								($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/')&&
								($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')&&
								(
								((length($temp[$H2[$y]+$whotoinclude[1]-2]) == 3)&&($sexes[$H2[$y]-1] == 1))
									||
								((length($temp[$H2[$y]+$whotoinclude[1]-2]) == 2)&&($sexes[$H2[$y]-1] == 0))
								)
								&&
								($temp[$H2[$w]+$whotoinclude[1]-2] ne './.')&&($temp[$H2[$w]+$whotoinclude[1]-2] ne './')&&
								($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/')&&
								($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/N')&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/')&&
								($temp[$H2[$w]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H2[$w]+$whotoinclude[1]-2] ne 'N/.')&&
								(
								((length($temp[$H2[$w]+$whotoinclude[1]-2]) == 3)&&($sexes[$H2[$w]-1] == 1))
									||
								((length($temp[$H2[$w]+$whotoinclude[1]-2]) == 2)&&($sexes[$H2[$w]-1] == 0))
								)


								){
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
							if(($temp[$H2[$y]+$whotoinclude[1]-2] ne './.')&&($temp[$H2[$y]+$whotoinclude[1]-2] ne './')&&
								($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/')&&
								($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/N')&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/')&&
								($temp[$H2[$y]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H2[$y]+$whotoinclude[1]-2] ne 'N/.')&&
								(
								((length($temp[$H2[$y]+$whotoinclude[1]-2]) == 3)&&($sexes[$H2[$y]-1] == 1))
									||
								((length($temp[$H2[$y]+$whotoinclude[1]-2]) == 2)&&($sexes[$H2[$y]-1] == 0))
								)
								){
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
							if(($temp[$H1[$y]+$whotoinclude[1]-2] ne './.')&&($temp[$H1[$y]+$whotoinclude[1]-2] ne './')&&
								($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/')&&
								($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/N')&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/')&&
								($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/.')&&
								(
								((length($temp[$H2[$y]+$whotoinclude[1]-2]) == 3)&&($sexes[$H2[$y]-1] == 1))
									||
								((length($temp[$H2[$y]+$whotoinclude[1]-2]) == 2)&&($sexes[$H2[$y]-1] == 0))
								)

								&& 
								($temp[$H1[$w]+$whotoinclude[1]-2] ne './.')&&($temp[$H1[$w]+$whotoinclude[1]-2] ne './')&&
								($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/')&&
								($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/N')&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/')&&
								($temp[$H1[$w]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H1[$w]+$whotoinclude[1]-2] ne 'N/.')&&
								(
								((length($temp[$H1[$w]+$whotoinclude[1]-2]) == 3)&&($sexes[$H1[$w]-1] == 1))
									||
								((length($temp[$H1[$w]+$whotoinclude[1]-2]) == 2)&&($sexes[$H1[$w]-1] == 0))
								)
								){
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
					# for some sites, there may be only one male individual, so no pi calculation is possible 
					if($num_comparisonsH1==0){
						for ($y = 0 ; $y <= $#H1; $y++ ) {
							if(($temp[$H1[$y]+$whotoinclude[1]-2] ne './.')&&($temp[$H1[$y]+$whotoinclude[1]-2] ne './')&&
								($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/*')&&($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/')&&
								($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/N')&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/')&&
								($temp[$H1[$y]+$whotoinclude[1]-2] ne '*/.')&&($temp[$H1[$y]+$whotoinclude[1]-2] ne 'N/.')&&
								(
								((length($temp[$H1[$y]+$whotoinclude[1]-2]) == 3)&&($sexes[$H1[$y]-1] == 1))
									||
								((length($temp[$H1[$y]+$whotoinclude[1]-2]) == 2)&&($sexes[$H1[$y]-1] == 0))
								)
								){
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
					# calculate the ABBA BABBA stats, plus more
					@temp1=split('',$string);
					$x_uniq = uniq @temp1;
					if($x_uniq == 2){
						$BBAA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_derived_freq*$H2_derived_freq*(1-$H3_derived_freq));
						if((((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq)>0) || (($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq)>0)){
							# this is a polymorphic position and it is an ABBA_BABA site   
							#if($current_window ==10000000 ){
							#	print $H1_derived_freq,"\t",$H2_derived_freq,"\t",$H3_derived_freq,"\n";
							#}
						    print $line,"\n";

							$ABBA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=((1-$H1_derived_freq)*$H2_derived_freq*$H3_derived_freq);
							$BABA_hash{$window_counter."_".$current_chromosome."_".$current_window}+=($H1_derived_freq*(1-$H2_derived_freq)*$H3_derived_freq);
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


print OUTFILE2 "chromosome\tbegin\tend\tABBA\tBABA\tBBAA\tD\tfdH2H3\tfH1H3\tf_dm\tdH2H3\tnum_sites_per_windowH2H3\tdH1H3\tnum_sites_per_windowH1H3\tH2pi\tnumsitesH2pi\tH1pi\tnumsitesH1pi\n";
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


And I generated a new R script to plot the 'stats' files generated by the scripts above (genome_plot.R). Theend of this script calculates the standard error using the block jackknife approach.  Here it is:
``` R
library (ggplot2)
setwd('/projects/SulaRADTAG/perl_scripts/sliding_window/genome_plot')

#pdf("sliding_plot.pdf",w=8, h=4, version="1.4", bg="transparent")
#data<-read.table("1_2-3-4-5-6-7_14-18-19-20.stats",header=TRUE)

#pdf("H1_allSulawesi_H2_maura_H3_Borneo.pdf",w=7, h=1.3, version="1.4", bg="transparent")
#data<-read.table("1-2-3-4-5-6-7-22-23-25-24-26-27-28-29-30-31-32-33-34-35-36-37-38-39-40_8-9-10-11-12-13_14-18-19-20.stats",header=TRUE)
 
#pdf("H1_allSulawesi_H2_nigra_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("1-2-3-4-5-6-7-8-9-10-11-12-13-24-26-27-28-29-30-31-32-33-34-35-36-37-38-39-40_22-23-25_14-18-19-20.stats",header=TRUE)

#pdf("H1_allSulawesi_H2_nigrescens_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("1-2-3-4-5-6-7-8-9-10-11-12-13-22-23-25-27-28-29-30-31-32-33-34-35-36-37-38-39-40_24-26_14-18-19-20.stats",header=TRUE)

#pdf("H1_allSulawesi_H2_hecki_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("1-8-9-10-11-12-13-22-23-25-24-26-27-28-29-30-31-32-33-34-35-36-37-38-39-40_2-3-4-5-6-7_14-18-19-20.stats",header=TRUE)

#pdf("H1_allSulawesi_H2_togeanus_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("32-33-34-35-36-37-38-39-40-27-28-29-24-26-22-23-25-8-9-10-11-12-13-2-3-4-5-6-7-1_30-31_14-18-19-20.stats",header=TRUE)

#pdf("H1_allSulawesi_H2_ochreata_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("32-33-34-35-36-37-38-39-40-30-31-24-26-22-23-25-8-9-10-11-12-13-2-3-4-5-6-7-1_27-28-29_14-18-19-20.stats",header=TRUE)

#pdf("H1_allSulawesi_H2_brun_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("32-33-34-35-36-37-38-39-40-30-31-27-28-29-24-26-22-23-25-8-9-10-11-12-13-2-3-4-5-6-7_1_14-18-19-20.stats",header=TRUE)

#pdf("H1_allSulawesi_H2_tonkeana_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("30-31-27-28-29-24-26-22-23-25-8-9-10-11-12-13-2-3-4-5-6-7-1_32-33-34-35-36-37-38-39-40_14-18-19-20.stats",header=TRUE)

pdf("H1_tonk_H2_nigra_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("32-33-34-35-36-37-38-39-40_22-23-25_14-18-19-20.stats",header=TRUE)
data<-read.table("H3born_H1tonk_H2nigra_5mil_windows_onlyRAD.stat",header=TRUE)

#pdf("H1_maura_H2_nigra_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("8-9-10-11-12-13_22-23-25_14-18-19-20.stats",header=TRUE)

#pdf("H1_hecki_H2_nigra_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("2-3-4-5-6-7_22-23-25_14-18-19-20.stats",header=TRUE)

#pdf("H1_nigrescens_H2_nigra_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("24-26_22-23-25_14-18-19-20.stats",header=TRUE)

#pdf("H1_ochreata_H2_nigra_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("27-28-29_22-23-25_14-18-19-20.stats",header=TRUE)

#pdf("H1_togeanus_H2_nigra_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("30-31_22-23-25_14-18-19-20.stats",header=TRUE)

#pdf("H1_brunn_H2_nigra_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("1_22-23-25_14-18-19-20.stats",header=TRUE)

#pdf("H1_maura_H2_heck_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("8-9-10-11-12-13_2-3-4-5-6-7_14-18-19-20.stats",header=TRUE)

#pdf("H1_tonk_H2_maura_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("32-33-34-35-36-37-38-39-40_8-9-10-11-12-13_14-18-19-20.stats",header=TRUE)

#pdf("H1_tonk_H2_hecki_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("32-33-34-35-36-37-38-39-40_2-3-4-5-6-7_14-18-19-20.stats",header=TRUE)

#pdf("H1_maura_H2_nigra_H3_Borneo.pdf",w=7, h=1.5, version="1.4", bg="transparent")
#data<-read.table("8-9-11-12-13_22-25_14-18-19-20.stats",header=TRUE)
#data$chromosome <- factor(data$chromosome,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chrX","chrY","chrM"),labels = c(1:20,"X","Y","M"))

pdf("H1_tonk_H2_nigra_H3_Borneo_HiSeqRAD.pdf",w=8, h=2.5, version="1.4", bg="transparent")
data<-read.table("combined.stats_noX",header=TRUE)
data$chromosome <- factor(data$chromosome,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chrX","chrY","chrM"),labels = c(1:20,"X","Y","M"))


# Make a new column withonly zeros that will be used for coloring dots
data$colors <- ">0.005"
# make values of the color column equal to one for low dH2H3
data$colors[!is.na(data$f_dm) & abs(data$dH2H3) <= 0.005] <- "<=0.005"
# plot the data
#my_y_title <- expression(italic("f[dm]")
d<-ggplot(data=data, aes(x=end,y=-f_dm, group=chromosome, colour=colors)) + geom_point(aes(colour = factor(colors)), alpha = 0.8) + 
  # combine the chromosomes
  facet_grid(~ chromosome, scales = "free_x", space = "free_x", margins = FALSE, drop = TRUE) + 
  # get rid of the gray background
  theme_bw() +
  # add an x-axis and y-axis label
  labs(x=expression(Chromosomes), y=expression(f[dm])) +
  # get rid of the space between the panels
  theme(panel.margin = unit(0, "lines")) + 
  # get rid of the x axis
  scale_x_discrete() +
  # rename the legend and make the dots specific colors
  scale_colour_manual(values = c("red","blue"), name = "divergence") +
  # add some horizontal blocks
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.7, ymax = 1.0), alpha = 0.002, linetype=0, fill="yellow") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -0.7, ymax = -1.0), alpha = 0.002, linetype=0, fill="yellow")+
  # add a horizontal line
  geom_hline(yintercept=0)
d
dev.off()




library(resample)
# first make a column with f_dm values weighted by the number of informative
# ABBA and BABA positions
# make a new column with the total number of ABBA_BABA sites
data$ABBA_BABA_sites<-data$BABA+data$ABBA
# calculate the mean number of ABBA_BABA sites
mean_number_of_ABBA_BABA<-mean(data$ABBA_BABA_sites)
# make a colum with teh f_dm values weighted by the number of ABBA_BABA sites
data$weighted_f_dm<-data$f_dm*(data$ABBA_BABA_sites/mean_number_of_ABBA_BABA)
jk<-jackknife(data, mean(weighted_f_dm, na.rm=TRUE), trace=TRUE)
var_replicates<-var(jk$replicates)
sterr_f_dm<-sqrt(length(which(!is.na(data$f_dm)))*var_replicates)
mean(data$weighted_f_dm, na.rm=TRUE)
sterr_f_dm
test_stat<-mean(data$weighted_f_dm, na.rm=TRUE)/sterr_f_dm
test_stat

data$weighted_D<-data$D*(data$ABBA_BABA_sites/mean_number_of_ABBA_BABA)
jk_D<-jackknife(data, mean(weighted_D, na.rm=TRUE), trace=TRUE)
var_replicates_D<-var(jk_D$replicates)
sterr_D<-sqrt(length(which(!is.na(data$weighted_D)))*var_replicates_D)
mean(data$weighted_D, na.rm=TRUE)
sterr_D
test_stat<-mean(data$weighted_D, na.rm=TRUE)/sterr_D
test_stat

```
