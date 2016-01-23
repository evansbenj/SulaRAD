# Subsetting the data and popgen stats.
I wrote a script to split up the master vcf file using the bedfiles I made for the 2014 paper.  I had to concatenate and sort some bed files to make a pooled file for the >51000 sites like this:

```bash
cat bedfile3_51000_to_101000.bed bedfile4_101000_to_151000.bed bedfile5_151000_and_higher.bed > bedfile_51000_and_higher.bed
sed -i -e 's/ /\t/g' bedfile_51000_and_higher.bed
sort -V -k 1,1 -k2,2 bedfile_51000_and_higher.bed > bedfile_51000_and_higher_sorted.bed
```
But this did not work so I instead ended up using the Galaxy webserver to merge the beds for regions >51000 bp away from genes from the 2014 paper.

And now I used this script to split up the vcf files into different sections depending on the distance of the sites from genes and then make tab delimited files (15_Executes_GATK_commands_SelectVariants_output_bed.pl):

```perl
#!/usr/bin/perl
# This script will read in a vcf file names and 
# make and execute a GATK commandline that divide up
# a vcf file into bunch of 
# new vcf files based on some bed files.

# then it will convert them to tab delimited format 


my $status;
my $infile = "final_round2_filtered.vcf";
#my $infile = "nonrecal_final_filtered.vcf";
my $outfile1 = "recal_plusminus_1000.vcf";
#my $outfile1 = "nonrecal_plusminus_1000.vcf";
my $bedfile1 = "bedfile1_genes_plusminus_1000.bed";
my $outfile2 = "recal_1000_51000.vcf";
#my $outfile2 = "nonrecal_1000_51000.vcf";
my $bedfile2 = "bedfile2_1001_to_51000.bed";
my $outfile3 = "recal_51000plus.vcf";
#my $outfile3 = "nonrecal_51000plus.vcf";
my $bedfile3 = "51000_and_more_galaxy.bed";

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


$commandline = "~/tabix-0.2.6/bgzip ".$outfile1;
$status = system($commandline);
$commandline = "~/tabix-0.2.6/tabix -p vcf ".$outfile1.".gz";
$status = system($commandline);
$commandline = "zcat ".$outfile1.".gz | /usr/local/vcftools/src/perl/vcf-to-tab > ".$outfile1.".gz.tab";
$status = system($commandline);

$commandline = "~/tabix-0.2.6/bgzip ".$outfile2;
$status = system($commandline);
$commandline = "~/tabix-0.2.6/tabix -p vcf ".$outfile2.".gz";
$status = system($commandline);
$commandline = "zcat ".$outfile2.".gz | /usr/local/vcftools/src/perl/vcf-to-tab > ".$outfile2.".gz.tab";
$status = system($commandline);

$commandline = "~/tabix-0.2.6/bgzip ".$outfile3;
$status = system($commandline);
$commandline = "~/tabix-0.2.6/tabix -p vcf ".$outfile3.".gz";
$status = system($commandline);
$commandline = "zcat ".$outfile3.".gz | /usr/local/vcftools/src/perl/vcf-to-tab > ".$outfile3.".gz.tab";
$status = system($commandline);

```

## Add outgroup sequences

And this can be followed up by adding the outgroup sequences using the script below and calculating the popgen stats using the script after that.

First, convert the vcf file to a tab delimited file like this:

``` 
~/tabix-0.2.6/bgzip final_filtered.vcf
~/tabix-0.2.6/tabix -p vcf final_filtered.vcf.gz
zcat final_filtered.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > final_filtered.vcf.gz.tab
```

And then use this script to add to this tab deimited file data from baboons (`16_Gets_outgroup_sequence_from_axt_files_NEW2015.pl`). This needs to be executed from with a directory that has axt alignment files. I also batch processed this script with another script called `17_adds_outgroup_to_lots_of_tab_files.pl`.  Here is the first script:

``` perl
#!/usr/bin/env perl
use strict;
use warnings;


# This program reads in a tab delimited file created by 
# vcftools and then extracts an outgroup sequence
# from all axt files that are in the directory

# It will then make a new file that has the outgroup
# sequence inserted in a column after the reference sequence.

# run it like this:
# Gets_outgroup_sequence_from_axt_files_NEW2015.pl in_tabfile out_tabfile

# the main concern here is that theaxt files have gaps inserted so that the 
# number of bases don't necessarily match the difference between the coordinates
# although hopefully the coordinates match.  So I need to count the gaps in the
# rhesus sequence and adjust the coordinates of the outgroup sequence appropriately

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";


my @files = glob("*.axt");

my @macaque_coordinates;
my @ingroup;
my @outgroup;
my $coordinate;
my %macaque_coordinates_key;
my $line;
my @line;
my @axt;
my $y;
my $switch;
my $source;
my $begin;
my $chr;


while ( my $line = <DATAINPUT>) {
	@macaque_coordinates=split("	",$line);
	if ($macaque_coordinates[0] ne "#CHROM"){
		#$macaque_coordinates_key{$macaque_coordinates[0]."_".$macaque_coordinates[1]}=$macaque_coordinates[2];
		$macaque_coordinates_key{$macaque_coordinates[0]."_".$macaque_coordinates[1]}="N";
	}
}

close DATAINPUT;
print "Done with input file 1\n";

foreach(@files){
	print $_,"\n";
	unless (open DATAINPUT1, $_) {
		print "Can not find the axt files, jackass.\n";
		exit;
	}

	LINE: while ( my $line1 = <DATAINPUT1>) {	
		@axt=split(" ",$line1);
		if(defined($axt[1])){
			if($axt[1] =~ /^chr/){
				$switch=1;
				$chr=$axt[1];
				$begin=$axt[2];
				next LINE;
			}
		}
		elsif((defined($axt[0]))&&($switch == 1)){
			# we need to find out where the gaps are in the reference sequence
			@ingroup=split("",$line1);	
			$switch = 2;
			next LINE;
		}
		elsif((defined($axt[0]))&&($switch == 2)){	
			$switch = 0;
			#print "papio ",$line1;
			@outgroup=split("",$line1);
			#print "yo ",$outgroup[0]," hey ",$#outgroup," ey ",$outgroup[$#outgroup-1],"\n";
			$coordinate=$begin-1;
			for ($y = 0 ; $y < $#outgroup ; $y++ ) {
				if($ingroup[$y] ne "-"){
					$coordinate+=1;
				}
				#print $chr."_".$coordinate,"\n";
				# check if this position is a gap in the ingroup
				if(defined($macaque_coordinates_key{$chr."_".$coordinate})){
					#print $chr,"_",$coordinate,"  ",$macaque_coordinates_key{$chr."_".$coordinate}," ";
					if($outgroup[$y] ne "-"){
						$macaque_coordinates_key{$chr."_".$coordinate} = $outgroup[$y];
					}
					else{
						$macaque_coordinates_key{$chr."_".$coordinate} = "N";
					}	
					#print $macaque_coordinates_key{$chr."_".$coordinate},"\n";
					if($y == 0){
						print "begin ",$chr,"_",$coordinate,"  ",$macaque_coordinates_key{$chr."_".$coordinate}," ";
						print $macaque_coordinates_key{$chr."_".$coordinate},"\n";
					}
					elsif($y == ($#outgroup)){
						print "end ",$chr,"_",$coordinate,"  ",$macaque_coordinates_key{$chr."_".$coordinate}," ";
						print $macaque_coordinates_key{$chr."_".$coordinate},"\n";
					}
				}
			}	
		}
	}
close DATAINPUT1;	
}

# now add the outgroup data to the output file in a new column
my @data; 
unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}
while ( my $line = <DATAINPUT>) {
	@data=split(/\s+/,$line);
	if ($data[0] ne "#CHROM"){
		for ($y = 0 ; $y <= 2 ; $y++ ) {
			print OUTFILE $data[$y],"\t";
		}
		if(defined($macaque_coordinates_key{$data[0]."_".$data[1]})){
			print OUTFILE $macaque_coordinates_key{$data[0]."_".$data[1]},"\t";	
		}
		else{
			print "Missed one; problem!\n";	
			print OUTFILE "N\t";	
		}	
		for ($y = 3 ; $y < $#data ; $y++ ) {
			print OUTFILE $data[$y],"\t";
		}
		print OUTFILE $data[$#data],"\n";
	}
	else{
		for ($y = 0 ; $y <= 2 ; $y++ ) {
			print OUTFILE $data[$y],"\t";
		}
		print OUTFILE "papAnu2\t";	
		for ($y = 3 ; $y < $#data ; $y++ ) {
			print OUTFILE $data[$y],"\t";
		}
		print OUTFILE $data[$#data],"\n";
	}
}


```
This works quite well and I have also used it to add the human outgroup sequence.  A small correction to the header is then required using sed:

`sed -i -e 's/papAnu2\tpapAnu2/hg19\tpapAnu2/g' final_round2_filtered.vcf.gz_with_baboon_and_human.tab`

# Automating addition of outgroups and calculation of popgen stats for each species

The script below adds the outgroup species can be executed for each tab file (there are 3 - one including genes plus 1000 bp on each end, one from 1001-51000 bp from genes, and one with the other bits). This needs to be done first for baboons, and then for humans and then a sed command to fix the header.  And for each of the three recal and non recal files (6 total). This can be done by pipine moving the tab delimited files into a folder and then globing them into the script I wrote for adding the outgroup.

Here is the script (`17_adds_outgroup_to_lots_of_tab_files.pl`):

``` perl
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
	$commandline = "16_Gets_outgroup_sequence_from_axt_files_NEW2015.pl ".$_." ".$_."_with_baboon.tab"
;
	$status = system($commandline);

}
```


## Popgen stats

I've polished up a script that calculates the important population statistics.  This script is quite flexible and it also does bootstrapping to provide confidence intervals (Boot_from_tab_diverge_poly_2015.pl):

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
my $number_of_bootstraps=1000;
my $lower=int($number_of_bootstraps*0.025);
my $upper=int($number_of_bootstraps*0.975);
my $JC_divergence_aDNA;
my $JC_divergence_xDNA;
my $JC_divergence_yDNA;
my $RAD_tag_count_aDNA;
my $RAD_tag_count_xDNA;
my $RAD_tag_count_yDNA;
my $distance_between_RAD_tags=500;
my $previousone_aDNA=0-$distance_between_RAD_tags;
my $previousone_xDNA=0-$distance_between_RAD_tags;
my $previousone_yDNA=0-$distance_between_RAD_tags;

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

print "Num automomal alleles: ",$n_aDNA,"\n";
print "Num chrX alleles: ",$n_xDNA,"\n";
print "Num chrY alleles: ",$n_yDNA,"\n";
my $expected_number_of_adna_pairwise_comparisons=0;
my $expected_number_of_xdna_pairwise_comparisons=0;
my $expected_number_of_ydna_pairwise_comparisons=0;

for ($y = 1 ; $y < $n_aDNA ; $y++ ) {
	$expected_number_of_adna_pairwise_comparisons+=$y;	
}

my $singleton_pi_aDNA = (($n_aDNA-1)/$expected_number_of_adna_pairwise_comparisons);
my $aDNA_singleton_sites;
my @aDNA_singleton_sites;

for ($y = 1 ; $y < $n_xDNA ; $y++ ) {
	$expected_number_of_xdna_pairwise_comparisons+=$y;	
}

my $singleton_pi_xDNA = (($n_xDNA-1)/$expected_number_of_xdna_pairwise_comparisons);
my $xDNA_singleton_sites;
my @xDNA_singleton_sites;

for ($y = 1 ; $y < $n_yDNA ; $y++ ) {
	$expected_number_of_ydna_pairwise_comparisons+=$y;	
}

my $singleton_pi_yDNA;
my $yDNA_singleton_sites;
my @yDNA_singleton_sites;

if($expected_number_of_ydna_pairwise_comparisons > 1){
	$singleton_pi_yDNA = (($n_yDNA-1)/$expected_number_of_ydna_pairwise_comparisons);
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
		if(($temp[0] ne "chrX")&&($temp[0] ne "chrY")&&($temp[0] ne "chrM")){
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
								if(($diff/$pi_counter) == $singleton_pi_aDNA){
									$aDNA_singleton_sites+=1;
									push(@aDNA_singleton_sites,1);
								}
								else{
									push(@aDNA_singleton_sites,0);
								}
								$diff=0;
								$aDNA_segregating_sites+=1;	
								push(@S_aDNA,1);					
								$pi_counter=0;	

						}
						if(($x_uniq == 1)||($x_uniq == 2)){
							$aDNA_sites+=1;
							if(uc $temp[$whotoinclude[0]-1] ne uc $temp1[0]){
								$aDNA_divergence+=1;
								push(@d_aDNA,1);
							}
							else{
								push(@d_aDNA,0);
							}
							if($temp[1] > ($previousone_aDNA+$distance_between_RAD_tags)){
								$RAD_tag_count_aDNA+=1;
							}
							$previousone_aDNA=$temp[1];
						}	
						if($aDNA_sites != ($#S_aDNA+1)){
							print $aDNA_sites,"a\t",$#S_aDNA+1,"a\t",$line,"a\n";
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
						if(($diff/$pi_counter) == $singleton_pi_xDNA){
							$xDNA_singleton_sites+=1;
							push(@xDNA_singleton_sites,1);
						}
						else{
							push(@xDNA_singleton_sites,0);
						}
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
						if($temp[1] > ($previousone_xDNA+$distance_between_RAD_tags)){
							$RAD_tag_count_xDNA+=1;
						}
						$previousone_xDNA=$temp[1];

					}	
					if($xDNA_sites != ($#S_xDNA+1)){
						print $xDNA_sites,"x\t",$#S_xDNA+1,"x\t",$line,"\n";
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
						if(($diff/$pi_counter) == $singleton_pi_yDNA){
							$yDNA_singleton_sites+=1;
							push(@yDNA_singleton_sites,1);
						}
						else{
							push(@yDNA_singleton_sites,0);
						}
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
						if($temp[1] > ($previousone_yDNA+$distance_between_RAD_tags)){
							$RAD_tag_count_yDNA+=1;
						}
						$previousone_yDNA=$temp[1];
					}	
					if($yDNA_sites != ($#S_yDNA+1)){
						print $yDNA_sites,"y\t",$#S_yDNA+1,"y\t",$line,"\n";
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

my @aDNA_singleton_sites_boot;
my @xDNA_singleton_sites_boot;
my @aDNA_singleton_sites_indexes;
my @xDNA_singleton_sites_indexes;

for ($m = 0 ; $m < $number_of_bootstraps ; $m++ ) {
	# generate an array with bootstrapped indexes
	# first aDNA
	print "bootstrap ",$m,"\n";
	
	for ($n = 0 ; $n < $aDNA_sites ; $n++ ) {
		push(@aDNA_bootstrapped_indexes,int(rand($aDNA_sites)));
	} # end $n

	# same for segregating sites
	for ($n = 0 ; $n < $aDNA_segregating_sites ; $n++ ) {
		push(@aDNA_singleton_sites_indexes,int(rand($aDNA_segregating_sites)));
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

	if((($e1_obs_aDNA*$S_aDNA_boot[$m])+($e2_obs_aDNA*$S_aDNA_boot[$m]*($S_aDNA_boot[$m]-1))) > 0){
		push(@aDNA_TajD_boot,($pi_aDNA_boot[$m] - ($S_aDNA_boot[$m]/$a1_obs_aDNA))/((($e1_obs_aDNA*$S_aDNA_boot[$m])+($e2_obs_aDNA*$S_aDNA_boot[$m]*($S_aDNA_boot[$m]-1)))**0.5));
 	}

	# singleton bootstrap
 	for ($n = 0 ; $n < $aDNA_segregating_sites ; $n++ ) {
		#calculate number of singletons for each replicate
		$aDNA_singleton_sites_boot[$m]+=$aDNA_singleton_sites[$aDNA_singleton_sites_indexes[$n]]
	} # end $n
	
	# now convert to a proportion
	$aDNA_singleton_sites_boot[$m]=($aDNA_singleton_sites_boot[$m]/$aDNA_segregating_sites);
  

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

	# same for segregating sites X
	for ($n = 0 ; $n < $xDNA_segregating_sites ; $n++ ) {
		push(@xDNA_singleton_sites_indexes,int(rand($xDNA_segregating_sites)));
	} # end $n

	push(@pi_xDNA_persite_boot,$pi_xDNA_boot[$m]/$xDNA_sites);

	# apply JC correction to divergence
	$JC_AP_divergence_xDNA_boot[$m]  = (-3/4)*log(1-(4/3)*($d_xDNA_boot[$m]/$xDNA_sites));
	# apply AP correction to divergence
	$JC_AP_divergence_xDNA_boot[$m]  = $JC_AP_divergence_xDNA_boot[$m] - ($pi_xDNA_boot[$m]/$xDNA_sites);

	# calculate pi/divergence for bootstrapped data
	push(@thetapi_over_divergence_xDNA,($pi_xDNA_boot[$m]/$xDNA_sites)/$JC_AP_divergence_xDNA_boot[$m]);
	# calculate TajD for bootstrapped data
	if((($e1_obs_xDNA*$S_xDNA_boot[$m])+($e2_obs_xDNA*$S_xDNA_boot[$m]*($S_xDNA_boot[$m]-1))) > 0){
		push(@xDNA_TajD_boot,($pi_xDNA_boot[$m] - ($S_xDNA_boot[$m]/$a1_obs_xDNA))/((($e1_obs_xDNA*$S_xDNA_boot[$m])+($e2_obs_xDNA*$S_xDNA_boot[$m]*($S_xDNA_boot[$m]-1)))**0.5));
	}

	# singleton bootstrap
 	for ($n = 0 ; $n < $xDNA_segregating_sites ; $n++ ) {
		#calculate number of singletons for each replicate
		$xDNA_singleton_sites_boot[$m]+=$xDNA_singleton_sites[$xDNA_singleton_sites_indexes[$n]]
	} # end $n
	
	# now convert to a proportion
	$xDNA_singleton_sites_boot[$m]=($xDNA_singleton_sites_boot[$m]/$xDNA_segregating_sites);


	# reset the indexes for next bootstrap
	@aDNA_bootstrapped_indexes=();
	@xDNA_bootstrapped_indexes=();
	@aDNA_singleton_sites_indexes=();
	@xDNA_singleton_sites_indexes=();



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


@aDNA_singleton_sites_boot = sort { $a <=> $b } @aDNA_singleton_sites_boot;
@xDNA_singleton_sites_boot = sort { $a <=> $b } @xDNA_singleton_sites_boot;



print OUTFILE "aDNA\n";
print OUTFILE "#_alleles\t",$n_aDNA,"\n";
print OUTFILE "#_Sites\t",$aDNA_sites,"\n";
print OUTFILE "RAD_tag_count\t",$RAD_tag_count_aDNA,"\n";
print OUTFILE "S\t",$aDNA_segregating_sites," (",$S_aDNA_boot[$lower]," - ",$S_aDNA_boot[$upper],")\n";
print OUTFILE "thetaW\t",sprintf("%.5f",$aDNA_segregating_sites/$a1_obs_aDNA/$aDNA_sites)," (",sprintf("%.5f",$S_aDNA_boot[$lower]/$a1_obs_aDNA/$aDNA_sites)," - ",sprintf("%.5f",$S_aDNA_boot[$upper]/$a1_obs_aDNA/$aDNA_sites),")\n";
print OUTFILE "pi\t",sprintf("%.5f",$asum/($#pi_aDNA+1))," (",sprintf("%.5f",$pi_aDNA_persite_boot[$lower])," - ",sprintf("%.5f",$pi_aDNA_persite_boot[$upper]),")\n";
my $pi_a = $asum/($#pi_aDNA+1);
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
#singleton proportion
print OUTFILE "Se/S\t",sprintf("%.3f",($aDNA_singleton_sites/$aDNA_segregating_sites))," (",sprintf("%.3f",($aDNA_singleton_sites_boot[$lower]))," - ",sprintf("%.3f",($aDNA_singleton_sites_boot[$upper])),")\n";

print OUTFILE "\n";
print OUTFILE "xDNA\n";
print OUTFILE "#_alleles\t",$n_xDNA,"\n";
print OUTFILE "#_Sites\t",$xDNA_sites,"\n";
print OUTFILE "RAD_tag_count\t",$RAD_tag_count_xDNA,"\n";
print OUTFILE "S\t",$xDNA_segregating_sites," (",$S_xDNA_boot[$lower]," - ",$S_xDNA_boot[$upper],")\n";
print OUTFILE "thetaW\t",sprintf("%.5f",$xDNA_segregating_sites/$a1_obs_xDNA/$xDNA_sites)," (",sprintf("%.5f",$S_xDNA_boot[$lower]/$a1_obs_xDNA/$xDNA_sites)," - ",sprintf("%.5f",$S_xDNA_boot[$upper]/$a1_obs_xDNA/$xDNA_sites),")\n";
print OUTFILE "pi\t",sprintf("%.5f",$xsum/($#pi_xDNA+1))," (",sprintf("%.5f",$pi_xDNA_persite_boot[$lower])," - ",sprintf("%.5f",$pi_xDNA_persite_boot[$upper]),")\n";
my $pi_x = $xsum/($#pi_xDNA+1);
print OUTFILE "d\t",sprintf("%.5f",$xDNA_divergence/$xDNA_sites)," (",sprintf("%.5f",$d_xDNA_boot[$lower]/$xDNA_sites)," - ",sprintf("%.5f",$d_xDNA_boot[$upper]/$xDNA_sites),")\n";
# apply JC correction
$JC_divergence_xDNA = (-3/4)*log(1-(4/3)*($xDNA_divergence/$xDNA_sites));
#print OUTFILE "d_JC\t",sprintf("%.5f",$JC_divergence_xDNA),"\n";
# apply correction for ancestral polymorphism
$JC_divergence_xDNA = $JC_divergence_xDNA- $xsum/($#pi_xDNA+1);
print OUTFILE "d_JC_AP\t",sprintf("%.5f",$JC_divergence_xDNA)," (",sprintf("%.5f",$JC_AP_divergence_xDNA_boot[$lower])," - ",sprintf("%.5f",$JC_AP_divergence_xDNA_boot[$upper]),")\n";
print OUTFILE "pi/d_JC_AP\t",sprintf("%.3f",($xsum/($#pi_xDNA+1))/$JC_divergence_xDNA)," (",sprintf("%.5f",$thetapi_over_divergence_xDNA[$lower])," - ",sprintf("%.5f",$thetapi_over_divergence_xDNA[$upper]),")\n";
# TajD xDNA
if(($n_xDNA > 2)&&(((($e1_obs_xDNA*$xDNA_segregating_sites)+($e2_obs_xDNA*$xDNA_segregating_sites*($xDNA_segregating_sites-1))))>0)){
	$TajD_xDNA = ($xsum - ($xDNA_segregating_sites/$a1_obs_xDNA))/((($e1_obs_xDNA*$xDNA_segregating_sites)+($e2_obs_xDNA*$xDNA_segregating_sites*($xDNA_segregating_sites-1)))**0.5);
	print OUTFILE "TajD\t",sprintf("%.3f",$TajD_xDNA)," (",sprintf("%.3f",$xDNA_TajD_boot[$lower])," - ",sprintf("%.3f",$xDNA_TajD_boot[$upper]),")\n";
}
else{
	$TajD_xDNA = "undefined";
 	print OUTFILE "TajD\tundefined\n";
}
#singleton proportion
print OUTFILE "Se/S\t",sprintf("%.3f",($xDNA_singleton_sites/$xDNA_segregating_sites))," (",sprintf("%.3f",($xDNA_singleton_sites_boot[$lower]))," - ",sprintf("%.3f",($xDNA_singleton_sites_boot[$upper])),")\n";

print OUTFILE "\n";
if((($number_of_individuals_genotyped-$number_of_female_individuals_genotyped)>1)&&($yDNA_sites>0)){
	print OUTFILE "yDNA\n";
	print OUTFILE "#_alleles\t",$n_yDNA,"\n";
	print OUTFILE "#_Sites\t",$yDNA_sites,"\n";
	print OUTFILE "RAD_tag_count\t",$RAD_tag_count_yDNA,"\n";
	print OUTFILE "S\t",$yDNA_segregating_sites,"\n";
	print OUTFILE "thetaW\t",sprintf("%.5f",$yDNA_segregating_sites/$a1_obs_yDNA/$yDNA_sites),"\n";
	print OUTFILE "pi\t",sprintf("%.5f",$ysum/($#pi_yDNA+1)),"\n";
	my $pi_y = $ysum/($#pi_yDNA+1);
	print OUTFILE "d\t",sprintf("%.5f",$yDNA_divergence/$yDNA_sites),"\n";
	# apply JC correction
	$JC_divergence_yDNA = (-3/4)*log(1-(4/3)*($yDNA_divergence/$yDNA_sites));
	#print OUTFILE "d_JC\t",sprintf("%.5f",$JC_divergence_yDNA),"\n";
	# apply correction for ancestral polymorphism
	$JC_divergence_yDNA = $JC_divergence_yDNA- $ysum/($#pi_yDNA+1);
	print OUTFILE "d_JC_AP\t",sprintf("%.5f",$JC_divergence_yDNA),"\n";
	print OUTFILE "pi/d_JC_AP\t",sprintf("%.3f",($ysum/($#pi_yDNA+1))/$JC_divergence_yDNA),"\n";
	# TajD yDNA
	if(($n_yDNA>2)&&($a1_obs_yDNA>0)&&(($e1_obs_yDNA*$yDNA_segregating_sites)+($e2_obs_yDNA*$yDNA_segregating_sites*($yDNA_segregating_sites-1))>0)){
		$TajD_yDNA = ($ysum - ($yDNA_segregating_sites/$a1_obs_yDNA))/((($e1_obs_yDNA*$yDNA_segregating_sites)+($e2_obs_yDNA*$yDNA_segregating_sites*($yDNA_segregating_sites-1)))**0.5);
		print OUTFILE "TajD\t",sprintf("%.3f",$TajD_yDNA),"\n";
	}
	else{
		$TajD_yDNA = "undefined";
		print OUTFILE "TajD\t",$TajD_yDNA,"\n\n";
	}
	#singleton proportion
	if($yDNA_segregating_sites>0){
		print OUTFILE "Se/S\t",sprintf("%.3f",($yDNA_singleton_sites/$yDNA_segregating_sites)),"\n\n";
	}
}

print OUTFILE "ratio_of_X_to_A_pi_over_d\n";
print OUTFILE "pi_X/d_jc_ad_X/pi_a/d_jc_ad_a\t",sprintf("%.3f",(($xsum/($#pi_xDNA+1))/$JC_divergence_xDNA)/(($asum/($#pi_aDNA+1))/$JC_divergence_aDNA))," ";

# now calculate and print the confidence intervals for this statistic using the delta method
# as in Evans et al. 2014
my $scaled_X_over_A = (($xsum/($#pi_xDNA+1))/$JC_divergence_xDNA)/(($asum/($#pi_aDNA+1))/$JC_divergence_aDNA);
my $var_dA_over_dX = (($JC_divergence_aDNA**2)/($JC_divergence_xDNA**4))*(sqrt($JC_divergence_xDNA)/sqrt($xDNA_sites))**2+(1/($JC_divergence_xDNA**2))*(sqrt($JC_divergence_aDNA)/sqrt ($aDNA_sites))**2;
my $dA_over_dX = ($JC_divergence_aDNA/$JC_divergence_xDNA);
my $var_piA = ((($n_aDNA+1)*$pi_a)/((3*$n_aDNA-1)*$aDNA_sites))+((2*($n_aDNA**2+$n_aDNA+3)*($pi_a**2))/(9*$n_aDNA*($n_aDNA-1)))/$RAD_tag_count_aDNA;
my $var_piX = ((($n_xDNA+1)*$pi_x)/((3*$n_xDNA-1)*$xDNA_sites))+((2*($n_xDNA**2+$n_xDNA+3)*($pi_x**2))/(9*$n_xDNA*($n_xDNA-1)))/$RAD_tag_count_xDNA;
my $var_piX_over_piA = ($pi_x**2/$pi_a**4)*$var_piA+(1/($pi_a**2))*$var_piX;
my $stdev_scaled_X_over_A = sqrt ($scaled_X_over_A**2*$var_dA_over_dX+$dA_over_dX**2*$var_piX_over_piA);
print OUTFILE " (",sprintf("%.3f",($scaled_X_over_A-1.96*$stdev_scaled_X_over_A))," - ",sprintf("%.3f",($scaled_X_over_A+1.96*$stdev_scaled_X_over_A)),")\n";

if(($number_of_individuals_genotyped-$number_of_female_individuals_genotyped)>1){
	print OUTFILE "pi_Y/d_jc_ad_Y/pi_a/d_jc_ad_a is\t",(($ysum/($#pi_yDNA+1))/$JC_divergence_yDNA)/(($asum/($#pi_aDNA+1))/$JC_divergence_aDNA),"\n";
}

close DATAINPUT;
close OUTFILE;

```

