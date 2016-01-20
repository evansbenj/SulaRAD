# Automating addition of outgroups and calculation of popgen stats for each species

The script that adds the outgroup species can be executed for each tab file (there are 3 - one including genes plus 1000 bp on each end, one from 1001-51000 bp from genes, and one with the other bits). This needs to be done first for baboons, and then for humans and then a sed command to fix the header.  And for each of the three recal and non recal files (6 total). This can be done by pipine moving the tab delimited files into a folder and then globing them into the script I wrote for adding the outgroup.

```perl
#!/usr/bin/env perl
use strict;
use warnings;

# this program will read in all tab delimited files in a folder and
# add outgroup sequences to them using the script
# 16_Gets_outgroup_sequence_from_axt_files_NEW2015.pl
# this should be executed from within the folder containing the axt files

# once this is done for baboons, make a symbolic link to the baboon files and this 
# script in the folder containing the human axt files, and then run it again
# to add the human outgroup

my $status;
my @tabfiles = glob("recal*.vcf.gz.tab");
my $commandline;

foreach(@tabfiles){
	$commandline = "16_Gets_outgroup_sequence_from_axt_files_NEW2015.pl ".$_." ".$_."_with_baboon.tab";
	$status = system($commandline);

}

```



And then we need to calculate the popgen stats and boostrap CIs for each of the species/populations for each of these files.  There are the following species:
* tonkeana
* hecki
* nigra
* nigrescent
* maura
* togeanus
* brunnescens
* ochreata
* nemestrina (Borneo)
* nemestrina (Sumatra)
* pagensis
* nemestrina (Malaysia)

I plan to do this using only the humans as outgroups, the justification being that I can get divergence data for the Y chromosome as well. Using the baboon I would not be able to do this because there is no baboon yDNA sequence.

Here is a program that automates the execution of the Boot_from_tab_diverge_poly_2015.pl script on the previous page: (18_calculate_popgenstats.pl):

```perl
#!/usr/bin/env perl
use strict;
use warnings;

# this program will read in all tab delimited files in a folder and
# calculate popgen stats plus boostraps for all species for each
# file

my $status;
my @tabfiles = glob("/home/ben/2015_SulaRADtag/good_merged_samples/tab_relative_to_genez/*with_baboon.tab_and_human.tab");
my $bin_sex = "1111100110000111100011100110010100000000";
my $commandline;
my @what_what;
my @species = ("tonk", "heck", "nigra", "nigres", "och", "maura", "brun", "tog","nemBor", "nemSum", "pag", "Malay");
my @numberz = ("32_33_34_35_36_37_38_39_40","2_3_4_5_6_7","22_23_25","24_26","27_28_29","8_9_10_11_12_13","1","30_31","18_19_20_21","15_17","21","16");
my $y;

foreach(@tabfiles){
	@what_what = split(".vcf.gz.",$_);
	for ($y = 0 ; $y <= $#species ; $y++ ) {
		$commandline = "Boot_from_tab_diverge_poly_2015.pl ".$_." ".$bin_sex." 4_6_".$numberz[$y]." ".$what_what[0]."_".$species[$y]."_boot.poly";
		$status = system($commandline);
	}
}

```

