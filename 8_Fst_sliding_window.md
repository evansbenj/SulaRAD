# Analysis of introgression across the Makassar Strait

I'd like to calculate FST in a sliging window across the genome for comparisons between Borneo nemestrina and Sulawesi macaques and for comparisons between Borneo nemestrina and individual Sulawesi macaques, such as tonkeana, hecki, and maura.  Vcftools looks like the most straightforward way to do this.  I first made a file containing names of all of the individual populations I wanted to compare.  These have a file name ending with "_names" in the `/home/ben/2015_SulaRADtag/good_merged_samples` directory on info.  Using the commands below, vcftools should generate a file with the suffix `weir.fst`.


```bash
~/tabix-0.2.6/bgzip final_round2_filtered.vcf
~/tabix-0.2.6/tabix -p vcf final_round2_filtered.vcf.gz
/usr/local/vcftools/src/cpp/vcftools --gzvcf final_round2_filtered.vcf.gz --weir-fst-pop borneo_names --weir-fst-pop tonkeana_names --fst-window-size 100000 --fst-window-step 100000
```

My initial analyses suggest that there is not much obvious signal from Fst.  Based on the ABBA_BABA analysis below, this is probably because only a subset of the hecki populations have a strong indication of geneflow with nemestrina, and this may be too subtle to pick up with a population based comparison.

#ANGSD analysis

The ABBA BABA test is much cooler and the results are not ambiguous.  ```

## Making a bed file to use to filter a bam file

First I used a perl script to make a bed file from my vcf files:
```
./vcf2bed.pl bad_sex_bad_het.vcf > bad_sex_bad_het.bed
./vcf2bed.pl round2_indels_only.vcf > round2_indels_only.bed
```
then I had to remove an extra column:
```
awk '{print $1 "\t" $2 "\t" $3}' bad_sex_bad_het.bed > bad_sex_bad_het_.bed
awk '{print $1 "\t" $2 "\t" $3}' round2_indels_only.bed > round2_indels_only_.vcf
mv bad_sex_bad_het_.bed bad_sex_bad_het.bed
mv round2_indels_only_.vcf round2_indels_only.bed
```
then I made bed files for flanking regions:
```
~/bedtools2/bin/bedtools flank -b 200 -i bad_sex_bad_het.bed -g rhesus_chromosome_lengths > bad_sex_bad_het_100bp_buffer.bed
~/bedtools2/bin/bedtools flank -b 3 -i round2_indels_only.vcf -g rhesus_chromosome_lengths > round2_indels_only_3bp_buffer.bed
```
I made this tab delimited `rhesus_chromosome_lengths` file from the `.dict` genome file:
```
chr1	228252215
chr2	189746636
chr3	196418989
chr4	167655696
chr5	182086969
chr6	178205221
chr7	169801366
chr8	147794981
chr9	133323859
chr10	94855758
chr11	134511895
chr12	106505843
chr13	138028943
chr14	133002572
chr15	110119387
chr16	78773432
chr17	94452569
chr18	73567989
chr19	64391591
chr20	88221753
chrX	153947521
chrY	11339300
chrM	16564
```
Then I did a bunch of steps to sort and merge these files into one file:

```
sort -k1,1 -k2,2n bad_sex_bad_het.bed > bad_sex_bad_het_sorted.bed
sort -k1,1 -k2,2n bad_sex_bad_het_200bp_buffer.bed > bad_sex_bad_het_200bp_buffer_sorted.bed
sort -k1,1 -k2,2n round2_indels_only.bed > round2_indels_only_sorted.bed 
sort -k1,1 -k2,2n round2_indels_only_3bp_buffer.bed > round2_indels_only_3bp_buffer_sorted.bed
cat bad_sex_bad_het_sorted.bed bad_sex_bad_het_200bp_buffer_sorted.bed round2_indels_only_sorted.bed round2_indels_only_3bp_buffer_sorted.bed > merge.bed
sort -k1,1 -k2,2n merge.bed > merge_sorted.bed
cat merge_sorted.bed | ~/bedtools2/bin/bedtools merge -i stdin > bad_sex_bad_het_and_200bp_buffer_and_round2_indels_and_3bp_buffer_sorted.bed
```

And then (finally) I filtered the bam file like this:

```
~/ngsutils-ngsutils-0.5.7/bin/bamutils filter recal_stampy_round2_all.bam recal_stampy_round2_all_filtered.bam -excludebed bad_sex_bad_het_and_200bp_buffer_and_round2_indels_and_3bp_buffer_sorted.bed
```
Then I divided up the filtered bam file into individual bam files like this:

```
samtools view -bhl brunescens_PF707_stampy_sorted recal_stampy_round2_all_filtered.bam > brunescens_PF707_stampy_sorted_filtered.bam
samtools view -bhl hecki_PF643_stampy_sorted recal_stampy_round2_all_filtered.bam > hecki_PF643_stampy_sorted_filtered.bam
samtools view -bhl tonk_PM582_stampy_sorted recal_stampy_round2_all_filtered.bam > tonk_PM582_stampy_sorted_filtered.bam
samtools view -bhl nem_pagensis_stampy_sorted recal_stampy_round2_all_filtered.bam > nem_pagensis_stampy_sorted_filtered.bam
samtools view -bhl tonk_PM584_stampy_sorted recal_stampy_round2_all_filtered.bam > tonk_PM584_stampy_sorted_filtered.bam
samtools view -bhl tonk_PM592_stampy_sorted recal_stampy_round2_all_filtered.bam > tonk_PM592_stampy_sorted_filtered.bam
samtools view -bhl tonk_PM602_stampy_sorted recal_stampy_round2_all_filtered.bam > tonk_PM602_stampy_sorted_filtered.bam
samtools view -bhl hecki_PF644_stampy_sorted recal_stampy_round2_all_filtered.bam > hecki_PF644_stampy_sorted_filtered.bam
samtools view -bhl hecki_PF648_stampy_sorted recal_stampy_round2_all_filtered.bam > hecki_PF648_stampy_sorted_filtered.bam
samtools view -bhl hecki_PF651_stampy_sorted recal_stampy_round2_all_filtered.bam > hecki_PF651_stampy_sorted_filtered.bam
samtools view -bhl hecki_PM639_stampy_sorted recal_stampy_round2_all_filtered.bam > hecki_PM639_stampy_sorted_filtered.bam
samtools view -bhl hecki_PM645_stampy_sorted recal_stampy_round2_all_filtered.bam > hecki_PM645_stampy_sorted_filtered.bam
samtools view -bhl maura_PF615_stampy_sorted recal_stampy_round2_all_filtered.bam > maura_PF615_stampy_sorted_filtered.bam
samtools view -bhl maura_PF713_stampy_sorted recal_stampy_round2_all_filtered.bam > maura_PF713_stampy_sorted_filtered.bam
samtools view -bhl maura_PM613_stampy_sorted recal_stampy_round2_all_filtered.bam > maura_PM613_stampy_sorted_filtered.bam
samtools view -bhl maura_PM614_stampy_sorted recal_stampy_round2_all_filtered.bam > maura_PM614_stampy_sorted_filtered.bam
samtools view -bhl maura_PM616_stampy_sorted recal_stampy_round2_all_filtered.bam > maura_PM616_stampy_sorted_filtered.bam
samtools view -bhl maura_PM618_stampy_sorted recal_stampy_round2_all_filtered.bam > maura_PM618_stampy_sorted_filtered.bam
samtools view -bhl nem_Kedurang_stampy_sorted recal_stampy_round2_all_filtered.bam > nem_Kedurang_stampy_sorted_filtered.bam
samtools view -bhl nem_Malay_stampy_sorted recal_stampy_round2_all_filtered.bam > nem_Malay_stampy_sorted_filtered.bam
samtools view -bhl nem_Ngasang_stampy_sorted recal_stampy_round2_all_filtered.bam > nem_Ngasang_stampy_sorted_filtered.bam
samtools view -bhl nem_Gumgum_stampy_sorted recal_stampy_round2_all_filtered.bam > nem_Gumgum_stampy_sorted_filtered.bam
samtools view -bhl nem_PM664_stampy_sorted recal_stampy_round2_all_filtered.bam > nem_PM664_stampy_sorted_filtered.bam
samtools view -bhl nem_PM665_stampy_sorted recal_stampy_round2_all_filtered.bam > nem_PM665_stampy_sorted_filtered.bam
samtools view -bhl nem_Sukai_male_stampy_sorted recal_stampy_round2_all_filtered.bam > nem_Sukai_male_stampy_sorted_filtered.bam
samtools view -bhl nigra_PF1001_stampy_sorted recal_stampy_round2_all_filtered.bam > nigra_PF1001_stampy_sorted_filtered.bam
samtools view -bhl nigra_PF660_stampy_sorted recal_stampy_round2_all_filtered.bam > nigra_PF660_stampy_sorted_filtered.bam
samtools view -bhl nigra_PM1000_stampy_sorted recal_stampy_round2_all_filtered.bam > nigra_PM1000_stampy_sorted_filtered.bam
samtools view -bhl nigra_PM1003_stampy_sorted recal_stampy_round2_all_filtered.bam > nigra_PM1003_stampy_sorted_filtered.bam
samtools view -bhl nigrescens_PF654_stampy_sorted recal_stampy_round2_all_filtered.bam > nigrescens_PF654_stampy_sorted_filtered.bam
samtools view -bhl ochreata_PF625_stampy_sorted recal_stampy_round2_all_filtered.bam > ochreata_PF625_stampy_sorted_filtered.bam
samtools view -bhl ochreata_PM571_stampy_sorted recal_stampy_round2_all_filtered.bam > ochreata_PM571_stampy_sorted_filtered.bam
samtools view -bhl ochreata_PM596_stampy_sorted recal_stampy_round2_all_filtered.bam > ochreata_PM596_stampy_sorted_filtered.bam
samtools view -bhl togeanus_PF549_stampy_sorted recal_stampy_round2_all_filtered.bam > togeanus_PF549_stampy_sorted_filtered.bam
samtools view -bhl togeanus_PM545_stampy_sorted recal_stampy_round2_all_filtered.bam > togeanus_PM545_stampy_sorted_filtered.bam
samtools view -bhl tonk_PF515_stampy_sorted recal_stampy_round2_all_filtered.bam > tonk_PF515_stampy_sorted_filtered.bam
samtools view -bhl tonk_PM561_stampy_sorted recal_stampy_round2_all_filtered.bam > tonk_PM561_stampy_sorted_filtered.bam
samtools view -bhl tonk_PM565_stampy_sorted recal_stampy_round2_all_filtered.bam > tonk_PM565_stampy_sorted_filtered.bam
samtools view -bhl tonk_PM566_stampy_sorted recal_stampy_round2_all_filtered.bam > tonk_PM566_stampy_sorted_filtered.bam
samtools view -bhl tonk_PM567_stampy_sorted recal_stampy_round2_all_filtered.bam > tonk_PM567_stampy_sorted_filtered.bam
```
## ABBA_BABA test

Now for Sumatra, Borneo, and pagensis:

```
angsd -out all_samples -doAbbababa 1 -doCounts 1 -b all_samples_bamfilez -anc /home/ben/2015_BIO720/rhesus_genome/macaque_masked_chromosomes_ym.fasta
Rscript R/jackKnife.R file=all_samples.abbababa indNames=all_samples_bamfilez outfile=all_samples_abbababa_out
```
where this is the `all_samples_bamfilez` file:

```
/home/ben/2015_SulaRADtag/good_merged_samples/nem_pagensis_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/nem_Kedurang_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/nem_Malay_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/nem_Ngasang_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/nem_Gumgum_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/nem_PM664_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/nem_PM665_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/nem_Sukai_male_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/brunescens_PF707_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/hecki_PF643_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/hecki_PF644_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/hecki_PF648_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/hecki_PF651_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/hecki_PM639_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/hecki_PM645_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/maura_PF615_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/maura_PF713_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/maura_PM613_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/brunescens_PF707_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/hecki_PF643_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/hecki_PF644_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/hecki_PF648_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/hecki_PF651_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/hecki_PM639_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/hecki_PM645_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/maura_PF615_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/maura_PF713_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/maura_PM613_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/maura_PM614_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/maura_PM616_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/maura_PM618_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/nigra_PF1001_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/nigra_PF660_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/nigra_PM1003_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/nigra_PM1000_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/nigrescens_PF654_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/ochreata_PF625_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/ochreata_PM571_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/ochreata_PM596_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/togeanus_PF549_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/togeanus_PM545_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/tonk_PF515_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/tonk_PM561_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/tonk_PM565_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/tonk_PM566_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/tonk_PM567_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/tonk_PM582_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/tonk_PM584_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/tonk_PM592_stampy_sorted_filtered.bam
/home/ben/2015_SulaRADtag/good_merged_samples/tonk_PM602_stampy_sorted_filtered.bam
```

Then I wrote a script to parse the output of the R script.  This script (Parse_abba_baba.pl) counts up the significant comparisons that support gene flow for each individual using either Borneo, Sumatra, or pagensis as outgroups.  I used a modified one to do the same thing for Sumatra and Borneo and with pagensis as the outgroup (I will add Sulawesi as the outgroup also).

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;

# This program parses the output of angsd's abba_baba analysis
# generated by the R script jackKnife.R

# to run type this "Parse_abba_babba.pl inputfilename outputfilename"

#my $inputfile = "sumatra_sulawesi_abbababa_out.txt";

my $inputfile = $ARGV[0];
unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

my $outputfile = $ARGV[1];

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

my @borneo = ("nem_PM664_stampy_sorted_filtered.bam","nem_PM665_stampy_sorted_filtered.bam","nem_Sukai_male_stampy_sorted_filtered.bam","nem_Gumgum_stampy_sorted_filtered.bam");
my @sumatra = ("nem_Kedurang_stampy_sorted_filtered.bam","nem_Malay_stampy_sorted_filtered.bam","nem_Ngasang_stampy_sorted_filtered.bam");
my $pagensis = "nem_pagensis_stampy_sorted_filtered.bam";
my @sulawesi = ("brunescens_PF707_stampy_sorted_filtered.bam","hecki_PF643_stampy_sorted_filtered.bam","hecki_PF644_stampy_sorted_filtered.bam","hecki_PF648_stampy_sorted_filtered.bam","hecki_PF651_stampy_sorted_filtered.bam","hecki_PM639_stampy_sorted_filtered.bam","hecki_PM645_stampy_sorted_filtered.bam","maura_PF615_stampy_sorted_filtered.bam","maura_PF713_stampy_sorted_filtered.bam","maura_PM613_stampy_sorted_filtered.bam","maura_PM614_stampy_sorted_filtered.bam","maura_PM616_stampy_sorted_filtered.bam","maura_PM618_stampy_sorted_filtered.bam","nigra_PF1001_stampy_sorted_filtered.bam","nigra_PF660_stampy_sorted_filtered.bam","nigra_PM1003_stampy_sorted_filtered.bam","nigra_PM1000_stampy_sorted_filtered.bam","nigrescens_PF654_stampy_sorted_filtered.bam","ochreata_PF625_stampy_sorted_filtered.bam","ochreata_PM571_stampy_sorted_filtered.bam","ochreata_PM596_stampy_sorted_filtered.bam","togeanus_PF549_stampy_sorted_filtered.bam","togeanus_PM545_stampy_sorted_filtered.bam","tonk_PF515_stampy_sorted_filtered.bam","tonk_PM561_stampy_sorted_filtered.bam","tonk_PM565_stampy_sorted_filtered.bam","tonk_PM566_stampy_sorted_filtered.bam","tonk_PM567_stampy_sorted_filtered.bam","tonk_PM582_stampy_sorted_filtered.bam","tonk_PM584_stampy_sorted_filtered.bam","tonk_PM592_stampy_sorted_filtered.bam","tonk_PM602_stampy_sorted_filtered.bam");
my %borneohash;
my %borneohash_numcomparisons;
my %sumatrahash;
my %sumatrahash_numcomparisons;
my %pagensishash;
my %pagensishash_numcomparisons;

my $critical_value=3.1; #this is the critical value for anpha = 0.001
my @temp;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne 'H1'){
		if((in_array(\@borneo,$temp[2]) == 1)&&(in_array(\@sulawesi,$temp[0]) == 1)&&(in_array(\@sulawesi,$temp[1]) == 1)){
			if($temp[8] > $critical_value){
				# this supports ABBA; geneflow into H2
				$borneohash{$temp[1]}+=1;
			}
			elsif($temp[8] < -$critical_value){
				# this supports BABA; geneflow into H1
				$borneohash{$temp[0]}+=1;
			}
			$borneohash_numcomparisons{$temp[0]}+=1;
			$borneohash_numcomparisons{$temp[1]}+=1;
		}
		elsif((in_array(\@sumatra,$temp[2]) == 1)&&(in_array(\@sulawesi,$temp[0]) == 1)&&(in_array(\@sulawesi,$temp[1]) == 1)){
			if($temp[8] > $critical_value){
				# this supports ABBA; geneflow into H2
				$sumatrahash{$temp[1]}+=1;
			}
			elsif($temp[8] < -$critical_value){
				# this supports BABA; geneflow into H1
				$sumatrahash{$temp[0]}+=1;
			}
			$sumatrahash_numcomparisons{$temp[0]}+=1;
			$sumatrahash_numcomparisons{$temp[1]}+=1;
		}
		elsif(($temp[2] eq $pagensis)&&(in_array(\@sulawesi,$temp[0]) == 1)&&(in_array(\@sulawesi,$temp[1]) == 1)){
			if($temp[8] > $critical_value){
				# this supports ABBA; geneflow into H2
				$pagensishash{$temp[1]}+=1;
			}
			elsif($temp[8] < -$critical_value){
				# this supports BABA; geneflow into H1
				$pagensishash{$temp[0]}+=1;
			}
			$pagensishash_numcomparisons{$temp[0]}+=1;
			$pagensishash_numcomparisons{$temp[1]}+=1;
		}
	}	
}			

print "Borneo\n";
foreach(@sulawesi){
        print $_,"\t",;
        if(defined($borneohash{$_})){
        	print $borneohash{$_},"\t";
        }
        else{
        	print "0\t";
        }	
       	if(defined($borneohash_numcomparisons{$_})){
        	print $borneohash_numcomparisons{$_},"\n";
        }
        else{
        	print "0\n";
        }
}

print "\nSumatra\n";
foreach(@sulawesi){
        print $_,"\t",;
        if(defined($sumatrahash{$_})){
        	print $sumatrahash{$_},"\t";
        }
        else{
        	print "0\t";
        }	
      	if(defined($sumatrahash_numcomparisons{$_})){
        	print $sumatrahash_numcomparisons{$_},"\n";
        }
        else{
        	print "0\n";
        }
}

print "\npagensis\n";
foreach(@sulawesi){
        print $_,"\t",;
        if(defined($pagensishash{$_})){
        	print $pagensishash{$_},"\t";
        }
        else{
        	print "0\t";
        }	
      	if(defined($pagensishash_numcomparisons{$_})){
        	print $pagensishash_numcomparisons{$_},"\n";
        }
        else{
        	print "0\n";
        }
}




close DATAINPUT;
close OUTFILE;


sub in_array {
my ($arr,$search_for) = @_;
foreach my $value (@$arr) {
	return 1 if $value eq $search_for;
}
 	return 0;
}

```

Because ANGSD works with bam files, I wanted to check if a custom script would recover similar results.  I also found it weird that the tonkeana samples, which had much higher coverage than the others, always seemed like they had no evidence of gene flow when compared to the other samples.  So I wrote this script (`Performs_ABBA_BABA.pl`):

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;


#  This program reads in a tab delimited genotype file generated
#  by the perl program "Gets_outgroup_sequence_from_axt_files.pl"
#  or from vcftools and performs the ABBA_BABA test using four individuals
#  that will be specified in the command line. 


# to execute type BPerforms_ABBA_BABA.pl inputfile.tab 1111100110000111100011100110010100000000 
# 4_5_14_22_32 out.abbababa  
# where 1111100110000111100011100110010100000000 refers to whether or not each individual in the ingroup 
# in the vcf file is (1) or is not (0) female ,and 4_5_14_22_32 refers to taxa and columns involved in the 
# expected (H4(H3,(H1,H2))) relationship as described next.

# To keep the format the same as ANGSD, an ABBA site has the same SNP in H2 and H3
# and a an BABA site has the same SNP in H1 and H3

# In "4_5_14_22_32", the first number is the column number of the of the outgroup 
# (4 in this case).  The second number is the column  where the ingroup begins. These begin with 1.
# The other columns are the ingroup individuals as enumerated below

# Example:
# Performs_ABBA_BABA.pl recal_1000_51000.vcf.gz.tab_with_baboon.tab_and_human.tab 1111100110000111100011100110010100000000 3_6_14_2_32 out.abbababa

# This should include gumgum, hecki_PF643, and tonk PF515.  If gene flow occurred from gumgum to PF643, 

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




my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $input3 = $ARGV[2];
my $outputfile = $ARGV[3];

my $string;

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

my $sliding_window=5000000;
my $current_window=0;
my $current_chromosome="blah";

my @temp;
my @temp1;


my $w;
my $y;
my $x;
my @unique;
my $x_uniq;

my $number_of_individuals_genotyped=($#whotoinclude - 1);

print "This number should be 3: ",$number_of_individuals_genotyped,"\n";

my $number_of_female_individuals_genotyped;
for ($y = 2 ; $y <= $#whotoinclude ; $y++ ) {
	if($sexes[$whotoinclude[$y]-1] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}	

print "This includes ",$number_of_female_individuals_genotyped," female(s)\n";

my %ABBA_hash;
my %BABA_hash;


while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line);
	if($temp[0] ne '#CHROM'){
		# base the genomic location on the outgroup	
		# this could be changed later to always rely on the most closely related ingroup
		if($temp[0] ne $current_chromosome){
			$current_chromosome = $temp[0];
			$current_window = 0;
		}
		until($temp[1] < ($current_window+$sliding_window)){
			$current_window = $current_window+$sliding_window;
		}
		if(($temp[0] ne "chrX")&&($temp[0] ne "chrY")&&($temp[0] ne "chrM")){
			$string=();
			if((uc $temp[$whotoinclude[0]-1] eq "A")||(uc $temp[$whotoinclude[0]-1] eq "C")||(uc $temp[$whotoinclude[0]-1] eq "T")||(uc $temp[$whotoinclude[0]-1] eq "G")){
				# the outgroup is defined
				$w = uc $temp[$whotoinclude[0]-1];
				$string=$string.$w;
				# now load the autosomal ingroup data
				for ($y = 0 ; $y < $number_of_individuals_genotyped; $y++ ) {
					# load a random allele from each individual
					if(rand() < 0.5){
						# pick the first allele
						if($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ] ne '.'){
							$w = uc $temp[$whotoinclude[1]-2 + 2*$whotoinclude[$y+2]-1 ];
							$string=$string.$w;
						}
					}
					else{
						# pick the second allele
						if($temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]] ne '.'){
							$w = uc $temp[ $whotoinclude[1]-2 + 2*$whotoinclude[$y+2]];
							$string=$string.$w;
						}	
					}
				}
				if(defined($string)){
					@temp1=split('',$string);
					if(($#temp1 == ($number_of_individuals_genotyped))){
						# there are no missing data
						$x_uniq = uniq @temp1;
						if($x_uniq == 2){
							# this is a polymorphic position 
							if(($temp1[0] eq $temp1[2])&&($temp1[1] eq $temp1[3])){
								# this is an ABBA site
								$ABBA_hash{$current_chromosome."_".$current_window}+=1;
								# To keep the format the same as ANGSD, an ABBA site has the same SNP in H2 and H3
							}	
							elsif(($temp1[0] eq $temp1[3])&&($temp1[1] eq $temp1[2])){
								# this is an BABA site	
								$BABA_hash{$current_chromosome."_".$current_window}+=1;
								# a BABA site has the same SNP in H1 and H3
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

@common_keys = map  { $_->[0] }
             sort { $a->[1] <=> $b->[1] }
             map  { [$_, $_=~/(\d+)/] }
                 @common_keys;
#@common_keys = sort @common_keys;


foreach (@common_keys) {
	@temp1=split('_',$_);
	print OUTFILE $temp1[0],"\t",$temp1[1]+1,"\t",$temp1[1]+$sliding_window,"\t";
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

close OUTFILE;

```

And the results can be analyzed using the ANGSD Rscript like this:

```
Rscript R/jackKnife.R file=/home/ben/2015_SulaRADtag/good_merged_samples/out.abbababa indNames=temp_bamfilez outfile=temp_abbababa_out
```

So far my preliminary results did recover support for migration with pagensis and the Ngasang sample from Sumatra (near the Mentawais).  

I also wrote a wrapper for `Performs_ABBA_BABA.pl` to automate each comparison.  Here it is (`Wrapper_for_Performs_ABBA_BABA.pl`):

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;

# This program is a wrapper for the script "Performs_ABBA_BABA.pl"

# It will feed in commandlines for all of the comparisons I want to do,
# and then analyze the results with the ANGSD R script 'jackKnife.R'
# and then compile the results

# to run type this: "Wrapper_for_Performs_abba_babba.pl outputfile"

# outputfile is the name of the summary file



#my $outputfile = $ARGV[0];

#unless (open(OUTFILE, ">$outputfile"))  {
#	print "I can\'t write to $outputfile\n";
#	exit;
#}
#print "Creating output file: $outputfile\n";

my @borneo = ("nem_Gumgum_stampy_sorted","nem_PM664_stampy_sorted","nem_PM665_stampy_sorted","nem_Sukai_male_stampy_sorted");
my @borneo_numbers = ("14","18","19","20");
my @sumatra = ("nem_Kedurang_stampy_sorted","nem_Malay_stampy_sorted","nem_Ngasang_stampy_sorted");
my @sumatra_numbers=("15","16","17");
my $pagensis = "nem_pagensis_stampy_sorted";
my $pagensis_number=21;
my @sulawesi = ("brunescens_PF707_stampy_sorted","hecki_PF643_stampy_sorted","hecki_PF644_stampy_sorted","hecki_PF648_stampy_sorted","hecki_PF651_stampy_sorted","hecki_PM639_stampy_sorted","hecki_PM645_stampy_sorted","maura_PF615_stampy_sorted","maura_PF713_stampy_sorted","maura_PM613_stampy_sorted","maura_PM614_stampy_sorted","maura_PM616_stampy_sorted","maura_PM618_stampy_sorted","nigra_PF1001_stampy_sorted","nigra_PF660_stampy_sorted","nigrescens_PM1000_stampy_sorted","nigra_PM1003_stampy_sorted","nigrescens_PF654_stampy_sorted","ochreata_PF625_stampy_sorted","ochreata_PM571_stampy_sorted","ochreata_PM596_stampy_sorted","togeanus_PF549_stampy_sorted","togeanus_PM545_stampy_sorted","tonk_PF515_stampy_sorted","tonk_PM561_stampy_sorted","tonk_PM565_stampy_sorted","tonk_PM566_stampy_sorted","tonk_PM567_stampy_sorted","tonk_PM582_stampy_sorted","tonk_PM584_stampy_sorted","tonk_PM592_stampy_sorted","tonk_PM602_stampy_sorted");
my @sulawesi_numbers=("1","2","3","4","5","6","7","8","9","10","11","12","13","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40");
#my @sulawesi = ("brunescens_PF707_stampy_sorted_filtered.bam","hecki_PF643_stampy_sorted_filtered.bam","hecki_PF644_stampy_sorted_filtered.bam","hecki_PF648_stampy_sorted_filtered.bam","hecki_PF651_stampy_sorted_filtered.bam","hecki_PM639_stampy_sorted_filtered.bam","hecki_PM645_stampy_sorted_filtered.bam","maura_PF615_stampy_sorted_filtered.bam","maura_PF713_stampy_sorted_filtered.bam","maura_PM613_stampy_sorted_filtered.bam","maura_PM614_stampy_sorted_filtered.bam","maura_PM616_stampy_sorted_filtered.bam","maura_PM618_stampy_sorted_filtered.bam","nigra_PF1001_stampy_sorted_filtered.bam","nigra_PF660_stampy_sorted_filtered.bam","nigra_PM1003_stampy_sorted_filtered.bam","nigra_PM1000_stampy_sorted_filtered.bam","nigrescens_PF654_stampy_sorted_filtered.bam","ochreata_PF625_stampy_sorted_filtered.bam","ochreata_PM571_stampy_sorted_filtered.bam","ochreata_PM596_stampy_sorted_filtered.bam","togeanus_PF549_stampy_sorted_filtered.bam","togeanus_PM545_stampy_sorted_filtered.bam");
my %borneohash;
my %borneohash_numcomparisons;
my %sumatrahash;
my %sumatrahash_numcomparisons;
my %pagensishash;
my %pagensishash_numcomparisons;
my $sliding_window=5000000;
my $outgroup_number=3; # for the tab file with human and baboon outgroup, 3=rhesus, 4=human, and 5=baboon
my $ingroup_column_begin_number=6; # for the tab file with human and baboon outgroup, this is 6
my $infile_tab = "final_round2_filtered.vcf.gz_with_baboon_and_human.tab";
my $commandline;
my $status;

my $critical_value=3.1; #this is the critical value for anpha = 0.001
my @temp;
my $x;
my $y;
my $z;

# first do all pairs of Sulawesi macaques with each of the Borneo samples
for ($x = 0 ; $x < $#sulawesi_numbers ; $x++ ) { # This is the first Sulawesi sample
	for ($y = 0 ; $y <= $#borneo_numbers ; $y++ ) { # This is each of the Borneo samples
		for ($z = ($x+1) ; $z <= $#sulawesi_numbers ; $z++ ) { # This is the second Sulawesi sample
			$commandline = "Performs_ABBA_BABA.pl ".$infile_tab." 1111100110000111100011100110010100000000 ".$outgroup_number."_".$ingroup_column_begin_number."_".$borneo_numbers[$y]."_".$sulawesi_numbers[$x]."_".$sulawesi_numbers[$z]." ".$sulawesi_numbers[$x]."_".$sulawesi_numbers[$z]."_".$borneo_numbers[$y].".abbababa";
			$status = system($commandline);
			# now make a file with the names of the taxa
			my $outputfile2 = $sulawesi_numbers[$x]."_".$sulawesi_numbers[$z]."_".$borneo_numbers[$y].".names";
			unless (open(OUTFILE2, ">$outputfile2"))  {
				print "I can\'t write to $outputfile2\n";
				exit;
			}
			print "Creating output file: $outputfile2\n";
			print OUTFILE2 $borneo[$y],"\n",$sulawesi[$x],"\n",$sulawesi[$z],"\n";
			close OUTFILE2;
		}
	}
}

```

And I wrote a wrapper for the jacknife.R script as well (Wrapper_for_Abbababa_jacknife.pl):

```perl
#!/usr/bin/env perl
use strict;
use warnings;

# This program is a wrapper for the script "jacknife.R"

# to run type this: "Wrapper_for_Abbababa_jacknife.pl"



my @borneo = ("nem_Gumgum_stampy_sorted","nem_PM664_stampy_sorted","nem_PM665_stampy_sorted","nem_Sukai_male_stampy_sorted");
my @borneo_numbers = ("14","18","19","20");
my @sumatra = ("nem_Kedurang_stampy_sorted","nem_Malay_stampy_sorted","nem_Ngasang_stampy_sorted");
my @sumatra_numbers=("15","16","17");
my $pagensis = "nem_pagensis_stampy_sorted";
my $pagensis_number=21;
my @sulawesi = ("brunescens_PF707_stampy_sorted","hecki_PF643_stampy_sorted","hecki_PF644_stampy_sorted","hecki_PF648_stampy_sorted","hecki_PF651_stampy_sorted","hecki_PM639_stampy_sorted","hecki_PM645_stampy_sorted","maura_PF615_stampy_sorted","maura_PF713_stampy_sorted","maura_PM613_stampy_sorted","maura_PM614_stampy_sorted","maura_PM616_stampy_sorted","maura_PM618_stampy_sorted","nigra_PF1001_stampy_sorted","nigra_PF660_stampy_sorted","nigrescens_PM1000_stampy_sorted","nigra_PM1003_stampy_sorted","nigrescens_PF654_stampy_sorted","ochreata_PF625_stampy_sorted","ochreata_PM571_stampy_sorted","ochreata_PM596_stampy_sorted","togeanus_PF549_stampy_sorted","togeanus_PM545_stampy_sorted","tonk_PF515_stampy_sorted","tonk_PM561_stampy_sorted","tonk_PM565_stampy_sorted","tonk_PM566_stampy_sorted","tonk_PM567_stampy_sorted","tonk_PM582_stampy_sorted","tonk_PM584_stampy_sorted","tonk_PM592_stampy_sorted","tonk_PM602_stampy_sorted");
my @sulawesi_numbers=("1","2","3","4","5","6","7","8","9","10","11","12","13","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40");
my %borneohash;
my %borneohash_numcomparisons;
my %sumatrahash;
my %sumatrahash_numcomparisons;
my %pagensishash;
my %pagensishash_numcomparisons;
my $sliding_window=5000000;
my $outgroup_number=3; # for the tab file with human and baboon outgroup, 3=rhesus, 4=human, and 5=baboon
my $ingroup_column_begin_number=6; # for the tab file with human and baboon outgroup, this is 6
#my $infile_tab = "final_round2_filtered.vcf.gz_with_baboon_and_human.tab";
my $commandline;
my $status;

my $critical_value=3.1; #this is the critical value for anpha = 0.001
my @temp;
my $x;
my $y;
my $z;

# first do all pairs of Sulawesi macaques with each of the Borneo samples
for ($x = 0 ; $x < $#sulawesi_numbers ; $x++ ) { # This is the first Sulawesi sample
	for ($y = 0 ; $y <= $#borneo_numbers ; $y++ ) { # This is each of the Borneo samples
		for ($z = ($x+1) ; $z <= $#sulawesi_numbers ; $z++ ) { # This is the second Sulawesi sample
			$commandline = "Rscript jackKnife.R file=".$sulawesi_numbers[$x]."_".$sulawesi_numbers[$z]."_".$borneo_numbers[$y].".abbababa indNames=".$sulawesi_numbers[$x]."_".$sulawesi_numbers[$z]."_".$borneo_numbers[$y].".names outfile=".$sulawesi_numbers[$x]."_".$sulawesi_numbers[$z]."_".$borneo_numbers[$y].".out";
			$status = system($commandline);
		}
	}
}

```

Then type `cat *out.txt > concat.out`

and the concatenated file can be parsed with this script (Parse_bens_abba_baba_sulawesi.pl):

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;

# This program parses the output of angsd's abba_baba analysis
# generated by the R script jackKnife.R

# to run type this "Parse_abba_babba.pl inputfilename outputfilename"

#my $inputfile = "sumatra_sulawesi_abbababa_out.txt";

my $inputfile = $ARGV[0];
unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

my $outputfile = $ARGV[1];

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

my @borneo = ("nem_Gumgum_stampy_sorted","nem_PM664_stampy_sorted","nem_PM665_stampy_sorted","nem_Sukai_male_stampy_sorted");
my @sumatra = ("nem_Kedurang_stampy_sorted","nem_Malay_stampy_sorted","nem_Ngasang_stampy_sorted");
my $pagensis = "nem_pagensis_stampy_sorted";
my @sulawesi = ("brunescens_PF707_stampy_sorted","hecki_PF643_stampy_sorted","hecki_PF644_stampy_sorted","hecki_PF648_stampy_sorted","hecki_PF651_stampy_sorted","hecki_PM639_stampy_sorted","hecki_PM645_stampy_sorted","maura_PF615_stampy_sorted","maura_PF713_stampy_sorted","maura_PM613_stampy_sorted","maura_PM614_stampy_sorted","maura_PM616_stampy_sorted","maura_PM618_stampy_sorted","nigra_PF1001_stampy_sorted","nigra_PF660_stampy_sorted","nigrescens_PM1000_stampy_sorted","nigra_PM1003_stampy_sorted","nigrescens_PF654_stampy_sorted","ochreata_PF625_stampy_sorted","ochreata_PM571_stampy_sorted","ochreata_PM596_stampy_sorted","togeanus_PF549_stampy_sorted","togeanus_PM545_stampy_sorted","tonk_PF515_stampy_sorted","tonk_PM561_stampy_sorted","tonk_PM565_stampy_sorted","tonk_PM566_stampy_sorted","tonk_PM567_stampy_sorted","tonk_PM582_stampy_sorted","tonk_PM584_stampy_sorted","tonk_PM592_stampy_sorted","tonk_PM602_stampy_sorted");
my %borneohash;
my %borneohash_numcomparisons;
my %sumatrahash;
my %sumatrahash_numcomparisons;
my %pagensishash;
my %pagensishash_numcomparisons;

my $critical_value=4; #this is the critical value for anpha = 0.001
my @temp;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne 'H1'){
		if((in_array(\@borneo,$temp[2]) == 1)&&(in_array(\@sulawesi,$temp[0]) == 1)&&(in_array(\@sulawesi,$temp[1]) == 1)){
			if($temp[8] > $critical_value){
				# this supports ABBA; geneflow into H2
				$borneohash{$temp[1]}+=1;
				if ($temp[0] eq "maura_PM618_stampy_sorted_filtered.bam"){
					print $temp[1],"\t",$temp[2],"\t",$temp[8],"\n";
				}
				if ($temp[1] eq "maura_PM618_stampy_sorted_filtered.bam"){
					print $temp[0],"\t",$temp[2],"\t",$temp[8],"\n";
				}
			}
			elsif($temp[8] < -$critical_value){
				# this supports BABA; geneflow into H1
				$borneohash{$temp[0]}+=1;
				if ($temp[0] eq "maura_PM618_stampy_sorted_filtered.bam"){
					print $temp[1],"\t",$temp[2],"\t",$temp[8],"\n";
				}
				if ($temp[1] eq "maura_PM618_stampy_sorted_filtered.bam"){
					print $temp[0],"\t",$temp[2],"\t",$temp[8],"\n";
				}
			}
			$borneohash_numcomparisons{$temp[0]}+=1;
			$borneohash_numcomparisons{$temp[1]}+=1;
		}
		elsif((in_array(\@sumatra,$temp[2]) == 1)&&(in_array(\@sulawesi,$temp[0]) == 1)&&(in_array(\@sulawesi,$temp[1]) == 1)){
			if($temp[8] > $critical_value){
				# this supports ABBA; geneflow into H2
				$sumatrahash{$temp[1]}+=1;
			}
			elsif($temp[8] < -$critical_value){
				# this supports BABA; geneflow into H1
				$sumatrahash{$temp[0]}+=1;
			}
			$sumatrahash_numcomparisons{$temp[0]}+=1;
			$sumatrahash_numcomparisons{$temp[1]}+=1;
		}
		elsif(($temp[2] eq $pagensis)&&(in_array(\@sulawesi,$temp[0]) == 1)&&(in_array(\@sulawesi,$temp[1]) == 1)){
			if($temp[8] > $critical_value){
				# this supports ABBA; geneflow into H2
				$pagensishash{$temp[1]}+=1;
			}
			elsif($temp[8] < -$critical_value){
				# this supports BABA; geneflow into H1
				$pagensishash{$temp[0]}+=1;
			}
			$pagensishash_numcomparisons{$temp[0]}+=1;
			$pagensishash_numcomparisons{$temp[1]}+=1;
		}
	}	
}			

print "Borneo\n";
foreach(@sulawesi){
        print $_,"\t",;
        if(defined($borneohash{$_})){
        	print $borneohash{$_},"\t";
        }
        else{
        	print "0\t";
        }	
       	if(defined($borneohash_numcomparisons{$_})){
        	print $borneohash_numcomparisons{$_},"\n";
        }
        else{
        	print "0\n";
        }
}

print "\nSumatra\n";
foreach(@sulawesi){
        print $_,"\t",;
        if(defined($sumatrahash{$_})){
        	print $sumatrahash{$_},"\t";
        }
        else{
        	print "0\t";
        }	
      	if(defined($sumatrahash_numcomparisons{$_})){
        	print $sumatrahash_numcomparisons{$_},"\n";
        }
        else{
        	print "0\n";
        }
}

print "\npagensis\n";
foreach(@sulawesi){
        print $_,"\t",;
        if(defined($pagensishash{$_})){
        	print $pagensishash{$_},"\t";
        }
        else{
        	print "0\t";
        }	
      	if(defined($pagensishash_numcomparisons{$_})){
        	print $pagensishash_numcomparisons{$_},"\n";
        }
        else{
        	print "0\n";
        }
}




close DATAINPUT;
close OUTFILE;


sub in_array {
my ($arr,$search_for) = @_;
foreach my $value (@$arr) {
	return 1 if $value eq $search_for;
}
 	return 0;
}

```
