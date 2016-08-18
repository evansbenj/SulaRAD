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

and here is the chrX one (Wrapper_for_Performs_ABBA_BABA_Borneo_vs_allsulawesi_populations_onlychrX_mirror.pl):

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
