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

Then I wrote a script to parse the output of the R script.  This script (Parse_abba_baba.pl) counts up the significant comparisons that support gene flow for each individual.  I modified this to make it use different outgroups, such as M. pagensis instead of Borneo and Sumatra.  (Actually I think I will just add this)

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
my $critical_value=4;
my @temp;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne 'H1'){
		if((in_array(\@borneo,$temp[2]) == 1)&&(in_array(\@borneo,$temp[0]) == 0)&&(in_array(\@borneo,$temp[1]) == 0)&&(in_array(\@sumatra,$temp[0]) == 0)&&(in_array(\@sumatra,$temp[1]) == 0)){
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
		elsif((in_array(\@sumatra,$temp[2]) == 1)&&(in_array(\@sumatra,$temp[0]) == 0)&&(in_array(\@sumatra,$temp[1]) == 0)&&(in_array(\@borneo,$temp[0]) == 0)&&(in_array(\@borneo,$temp[1]) == 0)){
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

