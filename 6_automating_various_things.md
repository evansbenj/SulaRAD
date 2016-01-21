# Popgen stats


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

