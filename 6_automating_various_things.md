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

Once we have the output for each species from the script above, I summarized them using this script:

```perl
#!/usr/bin/env perl
use strict;
use warnings;

# this program will read in the output of 
# Boot_from_tab_diverge_poly_2015.pl and do several things
# First, it will make an almost publication ready table 
# for each species
# which could go in SupplInfo
# and then it will make a new file for ploting
# an overlay graph that can be fed into R for a nice 
# figure in the main text


my $status;
#my @tabfiles = glob("/home/ben/2015_SulaRADtag/good_merged_samples/tab_relative_to_genez/recal*poly*");
my @tabfiles = glob("./recal/recal*poly");
my $commandline;
my @what_what;
my @species = ("tonk", "heck", "nigra", "nigres", "och", "maura", "brun", "tog","nemBor", "nemSum", "pag", "Malay");
my $y;
my @line;
my %hash_data;
my $wherewhere;
my @wherewhere=("aDNA","xDNA","yDNA","ratio_of_X_to_A_pi_over_d");
my $file_path;
my @file_path_parts;
my $wherewherewhere;
my @wherewherewhere=("plusminus","1000_to_50000","51000plus");
my $species;
my @polystatnamez;
my $polystatnamez;
my $inputfile;
my $file_path_parts;


foreach(@tabfiles){
	@polystatnamez=();
	$inputfile=$_;
	unless (open DATAINPUT, $inputfile) {
		print "Can not find the input file.\n";
		exit;
	}
	$file_path=$_;
	@file_path_parts=split("_",$file_path);
	if($file_path_parts[1] eq "1000"){
		$wherewherewhere = "1000_to_50000";
	}
	elsif($file_path_parts[1] eq "51000plus"){
		$wherewherewhere = "51000plus";
	}	
	elsif($file_path_parts[1] eq "plusminus"){
		$wherewherewhere = "plusminus";
	}

	foreach(@species){
		if($file_path =~ /$_/){
			$species = $_;
			print "The species is ",$species,"\n";
		}
	}

	while ( my $line = <DATAINPUT>) {
		@line=split(/[\t\n]/,$line);
		if(defined($line[0])){
			if ($line[0] =~ /DNA/) {
				$wherewhere = $line[0];
			}
			elsif($line[0] =~ /ratio_of_X_to_A_pi_over_d/){
				$wherewhere = $line[0];
			}
			elsif($line[0] ne "\n"){
				$hash_data{$species."_".$wherewhere."_".$wherewherewhere."_".$line[0]}=$line[1];
				if(($wherewhere =~ /aDNA/)||($wherewhere =~ /ratio_of_X_to_A_pi_over_d/)){
					push(@polystatnamez,$line[0]);
				}
			}
		}
	}
}

# now print out several datasets in a tab delimited format

foreach(@species){
	my $outputfile = $_."poly.tab";
	unless (open(OUTFILE, ">$outputfile"))  {
		print "I can\'t write to $outputfile\n";
		exit;
	}
	print "Creating output file: $outputfile\n";
	$species=$_;
	print OUTFILE "$species\n";

	foreach(@wherewhere){
		$wherewhere=$_;
		print OUTFILE $wherewhere,"\t";
		if($wherewhere eq "aDNA"){
			foreach(@wherewherewhere){
				print OUTFILE $_,"\t";
			}
		}
		print OUTFILE "\n";	
		if($wherewhere ne "ratio_of_X_to_A_pi_over_d"){
			foreach(@polystatnamez){
				$polystatnamez=$_;
				if(($polystatnamez ne "pi_X/d_jc_ad_X/pi_a/d_jc_ad_a")&&($polystatnamez ne "pi_Y/d_jc_ad_Y/pi_a/d_jc_ad_a")){
					print OUTFILE $polystatnamez,"\t";
					foreach(@wherewherewhere){
						if(defined($hash_data{$species."_".$wherewhere."_".$_."_".$polystatnamez})){
							print OUTFILE $hash_data{$species."_".$wherewhere."_".$_."_".$polystatnamez},"\t";
						}
						else{
							print OUTFILE "NA\t";
						}
					}
					print OUTFILE "\n";
				}
			}
		}
		else{
			foreach(@polystatnamez){
				$polystatnamez=$_;
				if(($polystatnamez eq "pi_X/d_jc_ad_X/pi_a/d_jc_ad_a")||($polystatnamez eq "pi_Y/d_jc_ad_Y/pi_a/d_jc_ad_a")){
					print OUTFILE $polystatnamez,"\t";
					foreach(@wherewherewhere){
						if(defined($hash_data{$species."_".$wherewhere."_".$_."_".$polystatnamez})){
							print OUTFILE $hash_data{$species."_".$wherewhere."_".$_."_".$polystatnamez},"\t";
						}
						else{
							print OUTFILE "NA\t";
						}
					}
					print OUTFILE "\n";
				}
			}

		}	
	}
	close OUTFILE;
}	


# now print out the x/a ratio for each genomic interval (near genes, medium distance 
# from genes, far from genes for plotting with R

my @bits_n_pieces;
my @bits_n_pieces2;

my $outputfile = "X_to_A.tab";
	unless (open(OUTFILE, ">$outputfile"))  {
		print "I can\'t write to $outputfile\n";
		exit;
	}
	print "Creating output file: $outputfile\n";

foreach(@species){
	$species=$_;
	foreach(@wherewherewhere){
		$wherewherewhere=$_;
		@bits_n_pieces=();
		if(defined($hash_data{$species."_ratio_of_X_to_A_pi_over_d_".$_."_pi_X/d_jc_ad_X/pi_a/d_jc_ad_a"})){
			@bits_n_pieces=split(/[\s+()\-]+/,$hash_data{$species."_ratio_of_X_to_A_pi_over_d_".$_."_pi_X/d_jc_ad_X/pi_a/d_jc_ad_a"});
		}
		if(defined($bits_n_pieces[0])){
			print OUTFILE "$species\t$wherewherewhere\t",sprintf("%.3f",$bits_n_pieces[0]),"\t";
		}
		else{
			print OUTFILE "$species\t$wherewherewhere\tNA\t";
		}
		if(defined($bits_n_pieces[1])){
			print OUTFILE sprintf("%.3f",$bits_n_pieces[1]),"\t";
		}
		else{
			print OUTFILE "NA\t";
		}
		if(defined($bits_n_pieces[2])){
			print OUTFILE sprintf("%.3f",$bits_n_pieces[2]),"\n";
		}
		else{
			print OUTFILE "NA\n";
		}
	}
}	
close OUTFILE;

```

and then printed a pretty plot using these R commands:

```R
library (ggplot2)

pdf("recal_X_to_A_plot.pdf",w=8, h=4, version="1.4", bg="transparent")

#par(mfrow=c(3,2),mar=c(4, 4, 0.5, 1))


data<-read.table("recal_X_to_A.tab",header=FALSE)
pd <- position_dodge(0.3)
data$V2 <- factor(data$V2,levels=c("plusminus","1000_to_50000","51000plus"))
data$V1 <- factor(data$V1,levels=c("nemestrina(Malaysia)","nemestrina(Sumatra)","nemestrina(Borneo)","pagensis","nigra","nigrescens","hecki","tonkeana","togeanus","ochreata","brunnescens","maura"))
ggplot(data, aes(x=V2,y=V3,colour=V1,group=V1))+geom_errorbar(aes(ymin=V4,ymax=V5),width=.1, position=pd) + geom_line(size=.6, position=pd, alpha=0.3)+ geom_point(position=pd) + theme_bw() + geom_hline(yintercept=0.75) + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + ylab("X/A") + xlab("Distance relative to genes") + guides(color=guide_legend(title="Species")) + scale_x_discrete(labels = c("<1000bp","1000-51000bp",">51000bp"))
dev.off()
```
