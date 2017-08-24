# Making a phylogeny for each chromosome

Below is a script to convert the tab data into a nexus file for each chromosome (21_tab_to_nexus.pl).  This script is somewhat clunky because it is hardcoded to have three outgroup sequences:

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;



#  This program reads in a tab delimited genotype file generated
#  by the perl program '17_adds_outgroup_to_lots_of_tab_files.pl'
#  or from vcftools vcf_to_tab
#  and generates a nexus file that includes degenerate bases and gaps

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

unless (open DATAINPUT, $inputfile) {
	print 'Can not find the input file.\n';
	exit;
}


my @temp;
my @temp1;
my $current_chromosome;
my $previous_chromosome='blah';
my @names;
my %datahash;
my $y;
my $watisitnow;
my $count=0;
my $outputfile;

# Read in datainput file
while ( my $line = <DATAINPUT>) {
	chomp($line);
	#@temp=split(/[\/'\t']+/,$line);
	@temp=split('\t',$line);
	if($temp[0] eq '#CHROM'){
		@names=@temp;
		for ($y = 2; $y <= $#names; $y++ ) {
			$datahash{$names[$y]}=' ';
		}	
	}
	else{	
		# only print ones that are not microsats or indels
		if((length($temp[2]) == 1)&&(length($temp[3]) == 1)&&(length($temp[4]) == 1)){
			$current_chromosome=$temp[0];
			$count=$count+1;
			for ($y = 2 ; $y <= 4; $y++ ) {	
				$watisitnow = $datahash{$names[$y]}.$temp[$y];
				$datahash{$names[$y]} = $watisitnow;
				#print $watisitnow,'\n';
			}
			for ($y = 5 ; $y <= $#temp; $y++ ) {
				# for these, we need to use IUPAC codes
				if(($temp[$y] eq 'G/G')||($temp[$y] eq 'C/C')||($temp[$y] eq 'T/T')||($temp[$y] eq 'A/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.substr($temp[$y],0,1);
				}
				elsif($temp[$y] eq './.'){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
				}
				elsif(($temp[$y] eq 'C/T')||($temp[$y] eq 'T/C')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'Y';
				}
				elsif(($temp[$y] eq 'A/G')||($temp[$y] eq 'G/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'R';
				}
				elsif(($temp[$y] eq 'A/C')||($temp[$y] eq 'C/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'M';
				}
				elsif(($temp[$y] eq 'A/T')||($temp[$y] eq 'T/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'W';
				}
				elsif(($temp[$y] eq 'C/G')||($temp[$y] eq 'G/C')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'S';
				}
				elsif(($temp[$y] eq 'G/T')||($temp[$y] eq 'T/G')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'K';
				}
				else{
					print 'unknown genotype ',$temp[$y],'\n';
				}
			}
			if($previous_chromosome eq 'blah'){
				$previous_chromosome = $current_chromosome;
			}	
			elsif($previous_chromosome ne $current_chromosome){
				$outputfile = $previous_chromosome.".phy";
				unless (open(OUTFILE, ">$outputfile"))  {
					print "I can\'t write to $outputfile\n";
					exit;
				}
				print "Creating output file: $outputfile\n";
				print OUTFILE "#NEXUS\n";
				print OUTFILE "BEGIN DATA\;\nDIMENSIONS NTAX=",$#names-1," NCHAR=",$count,"\n";
				print OUTFILE "FORMAT DATATYPE=DNA  MISSING=? GAP=- \;\n";
				print OUTFILE "MATRIX\n";
				for ($y = 2; $y <= $#names; $y++ ) {
					print OUTFILE $names[$y],"\t\t",$datahash{$names[$y]},"\n";
				}
				print OUTFILE "\;\nEND\;";
				close OUTFILE;
				$count=0;
				for ($y = 2; $y <= $#names; $y++ ) {
					$datahash{$names[$y]}=' ';
				}
				$previous_chromosome = $current_chromosome;
			}
		}
	}
}		
$outputfile = $previous_chromosome.".phy";
unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";
print OUTFILE "#NEXUS\n";
print OUTFILE "BEGIN DATA\;\nDIMENSIONS NTAX=",$#names-1," NCHAR=",$count,"\;\n";
print OUTFILE "FORMAT DATATYPE=DNA  MISSING=? GAP=- \;\n";
print OUTFILE "MATRIX\n";
for ($y = 2; $y <= $#names; $y++ ) {
	print OUTFILE $names[$y],"\t\t",$datahash{$names[$y]},"\n";
}
print OUTFILE "\;\nEND\;";
close OUTFILE;
$count=0;
```

# Making a phylogeny for the whole genome

Here is a cooler script that converts the same tab delimited file to an interleaved nexus file with the entire dataset in it (21_tab_to_interleave_nexus.pl):

```
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;



#  This program reads in a tab delimited genotype file generated
#  by the perl program '17_adds_outgroup_to_lots_of_tab_files.pl'
#  and generates an interleaved nexus file that includes degenerate bases and gaps

# run it like this
# 21_tab_to_interleave_nexus.pl input.tab output_interleave.nxs

# it is hardcoded to expect one base from two outgroup sequences after the reference (rhesus) seq

# Notes: 
# ##### nigra_PM1000 is actually nigrescens_PM1000 #####

# ##### these samples have very low coverage (<10x): ####
# ##### nigrescens_PF654_sorted (7.33X) ####
# ##### maura_PM613_sorted (8.65X) ####
# ##### ochreata_PM596_sorted (9.02X)####
# ##### nigra_660_sorted (9.61X) ####
# ##### togeanus_PF549 (9.63X) ####



my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

unless (open DATAINPUT, $inputfile) {
	print 'Can not find the input file.\n';
	exit;
}


my @temp;
my @temp1;
my @names;
my %datahash;
my $y;
my $x;
my $watisitnow;
my $count=0;
my $interleave=0;

# Read in datainput file

# Ideally, I'd like to print out the interleaved sections as the vcf file is read
# and keep track of the number of bases, and print this out to screen at the end
# then this number could be added to the nxs file after the script runs


while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] eq '#CHROM'){
		@names=@temp;
		for ($y = 2; $y <= $#names; $y++ ) {
			$names[$y] =~ s/-//gi; # get rid of dashes in names
			$datahash{$names[$y]}=''; # initialize the hash
		}	
		# print preamble to output file
		unless (open(OUTFILE, ">$outputfile"))  {
			print "I can\'t write to $outputfile\n";
			exit;
		}
		print "Creating output file: $outputfile\n";
		print OUTFILE "#NEXUS\n\n";
		print OUTFILE "BEGIN DATA\;\nDIMENSIONS NTAX=",$#names-1," NCHAR= XXXX\;\n";
		print OUTFILE "FORMAT DATATYPE=DNA  MISSING=? INTERLEAVE GAP=- \;\n";
		print OUTFILE "MATRIX\n";
		print OUTFILE "\n";
	}
	else{	
		# only print ones that are not microsats or indels in the outgroups
		if($interleave>79){  # print this section of the data
			for ($y = 2; $y <= $#names; $y++ ) {
				print OUTFILE $names[$y],"\t\t",$datahash{$names[$y]},"\n";
			}
			print OUTFILE "\n";
			$interleave=0;
			# clear the hash
			for ($y = 2 ; $y <= $#names; $y++ ) {
				$datahash{$names[$y]}='';
			}	
		}
		if((length($temp[2]) == 1)&&(length($temp[3]) == 1)&&(length($temp[4]) == 1)){ # all the outgroup seqs are single bp
			$count=$count+1; # this is the count of all positions
			$interleave=$interleave+1; # this is the count of the interleave length
			
			
			# now add data to the hash
			for ($y = 2 ; $y <= 4; $y++ ) {	# first the three outgroups which are haploid
				if($temp[$y] ne '*'){
					$datahash{$names[$y]} = $datahash{$names[$y]}.uc($temp[$y]);
				}
				else{
					$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
				}	
			}
			for ($y = 5 ; $y <= $#temp; $y++ ) { # now the ingroups, which are diploid, usually (except chrX and chrY)
				# for these, we need to use IUPAC codes
				if(($temp[$y] eq 'G/G')||($temp[$y] eq 'C/C')||($temp[$y] eq 'T/T')||($temp[$y] eq 'A/A')||($temp[$y] eq 'G/')||($temp[$y] eq 'C/')||($temp[$y] eq 'T/')||($temp[$y] eq 'A/')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.substr($temp[$y],0,1);
				}
				elsif(($temp[$y] eq './.')||($temp[$y] eq './')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
				}
				elsif(($temp[$y] eq 'C/T')||($temp[$y] eq 'T/C')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'Y';
				}
				elsif(($temp[$y] eq 'A/G')||($temp[$y] eq 'G/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'R';
				}
				elsif(($temp[$y] eq 'A/C')||($temp[$y] eq 'C/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'M';
				}
				elsif(($temp[$y] eq 'A/T')||($temp[$y] eq 'T/A')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'W';
				}
				elsif(($temp[$y] eq 'C/G')||($temp[$y] eq 'G/C')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'S';
				}
				elsif(($temp[$y] eq 'G/T')||($temp[$y] eq 'T/G')){
					$datahash{$names[$y]} = $datahash{$names[$y]}.'K';
				}
				else{ # this is a microsat, so substitute missing data
					$datahash{$names[$y]} = $datahash{$names[$y]}.'-';
				}
			}	
		}
	}
}		

# print the last line
for ($y = 2; $y <= $#names; $y++ ) {
	print OUTFILE $names[$y],"\t\t",$datahash{$names[$y]},"\n";
}
#print OUTFILE "\n";



print OUTFILE "\;\nEND\;";
print OUTFILE "\n";
print OUTFILE "\n";
#print OUTFILE "BEGIN Mrbayes\;\n";
#print OUTFILE "Prset statefreqpr=dirichlet(1,1,1,1)\;\n";
#print OUTFILE "Lset  nst=6  rates=invgamma\;\n";
#print OUTFILE "mcmc ngen=2000000 savebrlens=yes\;\n";
#print OUTFILE "sumt burnin=10000\;\n";
#print OUTFILE "quit\;\n";

close OUTFILE;
print "The number of sites is $count\n";

# now update the number of bases

my $status;
$status = system("perl -p -i -e 's/XXXX/$count/g' $outputfile");

```
# Phylogeny with iqtree

```
~/2015_SulaRADtag/good_merged_samples/iqtree-1.5.0a-Linux/bin/iqtree -s FINAL_RADseq_alldata_noBSQR_2016_haploiddepth_X_nomaleXhets.nxs -m TEST -nt 1 -pre FINAL_RADseq_alldata_noBSQR_2016_haploiddepth_X_nomaleXhets.nxs_
```
```
~/2015_SulaRADtag/good_merged_samples/iqtree-1.5.0a-Linux/bin/iqtree -s GenotypeVCFs_noBSQR_filtered_xDNA_only_genotyped_based_on_depth_with_baboon_and_human.nxs -m K3Pu+G4 -bb 1000
```

